subroutine WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
    use file_operations, only: out_i, outfile
    ! T, P and Density of the calculated envelope
    DOUBLE PRECISION, dimension(800) :: Tv
    DOUBLE PRECISION, dimension(800) :: Pv
    DOUBLE PRECISION, dimension(800) :: Dv

    ! number of valid elements in To, Po and Do arrays
    integer :: n_points

    ! positions of the last saturation points before each critical point
    integer, dimension(4) :: icri
    ! T, P and Density of critical points
    DOUBLE PRECISION, dimension(4) :: Tcri
    DOUBLE PRECISION, dimension(4) :: Pcri
    DOUBLE PRECISION, dimension(4) :: Dcri

    ! number of valid elements in icri, Tcri, Pcri and Dcri arrays
    integer :: ncri

    character(len=200) :: filename

    out_i = out_i + 1
    filename = "envelout"
    filename = outfile(filename, out_i)

    open(unit=out_i, file=filename)

    !WRITE (out_i, *)
    WRITE (out_i, *) '   T(K)        P(bar)        D(mol/L)'
    do i = 1, n_points
        WRITE (out_i, 1) Tv(i), Pv(i), Dv(i)
    end do
1   FORMAT(F12.4, 2E14.4, x, I4)
    WRITE (out_i, *)
    WRITE (out_i, *) ' Number of critical points found: ', ncri
    WRITE (out_i, *) '   T(K)        P(bar)        D(mol/L)'
    do i = 1, ncri
        WRITE (out_i, 1) Tcri(i), Pcri(i), Dcri(i), icri(i)
    end do
end subroutine WriteEnvel

subroutine CheckCross(XpairA, YpairA, XpairB, YpairB, Cross, Xcr, Ycr)
    ! For given consecutive pairs of points along two different curves A and B, this subroutine:
    ! 1) Finds the straight line equation connecting each pair of points ( Y = c*X + d )
    ! 2) Determines whether there is a crossing between the two segments
    implicit DOUBLE PRECISION(A - H, O - Z)
    DOUBLE PRECISION, dimension(2) :: XpairA, YpairA, XpairB, YpairB
    logical :: run, passingcri, Cross
    Cross = .false.
    cA = (YpairA(2) - YpairA(1))/(XpairA(2) - XpairA(1))
    dA = YpairA(1) - cA*XpairA(1)
    cB = (YpairB(2) - YpairB(1))/(XpairB(2) - XpairB(1))
    dB = YpairB(1) - cB*XpairB(1)
    Xcr = (dB - dA)/(cA - cB)
    if ((Xcr - XpairA(1))*(Xcr - XpairA(2)) < 0) then
        if ((Xcr - XpairB(1))*(Xcr - XpairB(2)) < 0) then
            Cross = .true.
            Ycr = cA*Xcr + dA
        end if
    end if
end subroutine CheckCross

subroutine envelope2(ichoice, model, n, z, T, P, KFACT, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)

    implicit DOUBLE PRECISION(A - H, O - Z)
    PARAMETER(nco=64)

    ! M&M means the book by Michelsen and Mollerup, 2nd Edition (2007)

    ! eos id, number of compounds in the system and starting point type
    integer, intent(in) :: model, n, ichoice

    DOUBLE PRECISION Kinf

    ! estimated T and P for first point (then used for every point)
    DOUBLE PRECISION :: T, P

    ! estimated K factors for first point (then used for every point)
    DOUBLE PRECISION, dimension(n) :: KFACT
    DOUBLE PRECISION, dimension(nco) :: KFcr1, Kscr1, KFcr2, Kscr2, KlowT, PHILOGxlowT ! go in COMMONS (cannot have n dimension)

    ! composition of the system
    DOUBLE PRECISION, dimension(n), intent(in) :: z

    ! pure compound physical constants
    DOUBLE PRECISION, dimension(n), intent(in) :: tcn
    DOUBLE PRECISION, dimension(n), intent(in) :: pcn
    DOUBLE PRECISION, dimension(n), intent(in) :: omgn

    ! eos parameters
    DOUBLE PRECISION, dimension(n), intent(in) :: acn  ! in bar*(L/mol)**2
    DOUBLE PRECISION, dimension(n), intent(in) :: bn   ! in L/mol
    DOUBLE PRECISION, dimension(n), intent(in) :: delta1n  !only required for RKPR
    DOUBLE PRECISION, dimension(n), intent(in) :: k_or_mn  ! k for RKPR ; m for SRK/PR

    ! interaction parameters matrices
    DOUBLE PRECISION, dimension(n, n), intent(in) :: Kij_or_K0n
    DOUBLE PRECISION, dimension(n, n), intent(in) :: Tstarn
    DOUBLE PRECISION, dimension(n, n), intent(in) :: Lijn

    ! T, P and Density of the calculated envelope
    DOUBLE PRECISION, dimension(800), intent(out) :: Tv
    DOUBLE PRECISION, dimension(800), intent(out) :: Pv
    DOUBLE PRECISION, dimension(800), intent(out) :: Dv

    ! number of valid elements in Tv, Pv and Dv arrays
    integer, intent(out) :: n_points

    ! positions of the last saturation points before each critical point
    integer, dimension(4), intent(out) :: icri
    ! T, P and Density of critical points
    DOUBLE PRECISION, dimension(4), intent(out) :: Tcri
    DOUBLE PRECISION, dimension(4), intent(out) :: Pcri
    DOUBLE PRECISION, dimension(4), intent(out) :: Dcri

    ! number of valid elements in icri, Tcri, Pcri and Dcri arrays
    integer, intent(out) :: ncri

    ! Intermediate variables during calculation process
    DOUBLE PRECISION, dimension(n) :: y, PHILOGy, PHILOGx
    DOUBLE PRECISION, dimension(n) :: DLPHITx, DLPHIPx, DLPHITy, DLPHIPy
    DOUBLE PRECISION, dimension(n, n) :: FUGNx, FUGNy
    integer, dimension(n + 2) :: ipiv
    DOUBLE PRECISION, dimension(n + 2) :: X, Xold, Xold2, delX, bd, F, dFdS, dXdS
    DOUBLE PRECISION, dimension(n + 2, n + 2) :: JAC, AJ
    DOUBLE PRECISION :: Vy, Vx
    DOUBLE PRECISION, dimension(2) :: TpairA, PpairA, TpairB, PpairB
    logical :: run, passingcri, Cross, minT, minmaxT

    DOUBLE PRECISION, dimension(nco, nco) :: Kij_or_K0, Tstar
    COMMON/CRIT/TC(nco), PC(nco), DCeos(nco), omg(nco)
    COMMON/COMPONENTS/ac(nco), b(nco), delta1(nco), rk_or_m(nco), Kij_or_K0, NTDEP
    COMMON/MODEL/NMODEL
    COMMON/rule/ncomb
    COMMON/bcross/bij(nco, nco)
    COMMON/Tdep/Kinf, Tstar
    COMMON/DewCurve/ilastDewC, TdewC(800), PdewC(800), dewK(800, nco)
    COMMON/CrossingPoints/Tcr1, Pcr1, Tcr2, Pcr2, KFcr1, Kscr1, KFcr2, Kscr2
    COMMON/lowTbub/TlowT, PlowT, KlowT, PHILOGxlowT

    minT = .false.
    minmaxT = .false.
    Told2 = 0.0
    Told = 10.0
    maxP = 0.d0
    ! Charging the commons(nco) from input arguments (n)
    NMODEL = model
    TC(:n) = tcn
    PC(:n) = pcn
    OMG(:n) = omgn
    ac(:n) = acn
    b(:n) = bn
    delta1(:n) = delta1n
    rk_or_m(:n) = k_or_mn
    Kij_or_K0(:n, :n) = Kij_or_K0n
    Kinf = 0.0d0
    ncomb = 0  ! only  vdW combining rules and quadratic mixing rules by  the moment
    Tstar(:n, :n) = Tstarn
    ! b matrix for Classical or van der Waals combining rules:
    do i = 1, n
        do j = i, n
            bij(i, j) = (1 - lijn(i, j))*(b(i) + b(j))/2
            bij(j, i) = bij(i, j)
        end do
    end do
    !
    !-----------------------------------------------------------

    ! Continuation method for tracing the envelope starts here
    run = .true.
    i = 0
    ncri = 0
    JAC(n + 1, :) = 0.d0
    lda = n + 2
    ldb = n + 2
    X(:n) = log(KFACT)
    X(n + 1) = log(T)
    X(n + 2) = log(P)
    iy = 1
    ix = 1
    if (ichoice <= 2) then  ! low T bub (1) or dew (2)
        if (ichoice == 1) iy = -1
        if (ichoice == 2) ix = -1  ! x will be vapor phase during the first part, and liquid after a critical point is crossed
        ns = n + 1
        S = log(T)
        delS = 0.005
        y = KFACT*z ! Wilson estimate for vapor (or liquid) composition
    else    ! (ichoice==3) high P L-L sat
        PmaxDewC = maxval(PdewC(1:ilastDewC))
        ns = n + 2
        S = log(P)
        delS = -0.005
        y = 0.d0
        y(n) = 1.d0
    end if
    Xold = 0.d0
    dFdS = 0.d0
    dFdS(n + 2) = -1.d0
    do while (run)
        i = i + 1
        ! Newton starts here
        delX = 1.0
        iter = 0
        do while (maxval(abs(delX)) > 1.d-5 .and. iter <= 70)
            iter = iter + 1
            !      nc,MTYP,INDIC,T,P,rn,V,PHILOG,DLPHI,DLPHIP,DLPHIT,FUGN
            call TERMO(n, iy, 4, T, P, y, Vy, PHILOGy, DLPHIPy, DLPHITy, FUGNy)
            call TERMO(n, ix, 2, T, P, z, Vx, PHILOGx, DLPHIPx, DLPHITx, FUGNx)
            F(:n) = X(:n) + PHILOGy - PHILOGx  ! X(:n) are LOG_K
            F(n + 1) = sum(y - z)
            F(n + 2) = X(ns) - S
            do j = 1, n
                JAC(1:n, j) = FUGNy(:, j)*y(j)  ! z*K=y
                JAC(j, j) = JAC(j, j) + 1.d0
            end do
            JAC(1:n, n + 1) = T*(DLPHITy - DLPHITx)
            JAC(1:n, n + 2) = P*(DLPHIPy - DLPHIPx)
            JAC(n + 1, 1:n) = y
            JAC(n + 2, :) = 0.d0
            JAC(n + 2, ns) = 1.d0
            ! call dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
            bd = -F
            AJ = JAC
            call dgesv(n + 2, 1, AJ, lda, ipiv, bd, ldb, info)
            if (info .ne. 0) write (6, *) "error with dgesv in parameter ", info
            delX = bd
            if (i == 1) then
                do while (maxval(abs(delX)) > 5.0)   ! Too large Newton step --> Reduce it
                    delX = delX/2
                end do
            else
                do while (maxval(abs(delX)) > 0.08)   ! Too large Newton step --> Reduce it
                    delX = delX/2
                end do
                if (iter > 8) delX = delX/2  ! too many iterations (sometimes due to oscillatory behavior near crit) --> Reduce it
            end if
            X = X + delX
            if (.not. passingcri .and. i /= 1 .and. iter == 8 .and. maxval(abs(delX)) > 0.001) then ! Too many iterations-->Reduce step to new point
                delS = delS/2
                S = S - delS
                X = Xold + dXdS*delS
                iter = 0
            end if
            KFACT = exp(X(:n))
            y = z*KFACT
            T = exp(X(n + 1))
            P = exp(X(n + 2))
        end do
        ! Point converged (unless it jumped out because of high number of iterations)

        if (ichoice == 2) dewK(i, :n) = X(:n)
        if (ichoice == 1 .and. i == 1) then ! save for starting later a 3-envel
            TlowT = T
            PlowT = P
            KlowT(1:n) = KFACT
            PHILOGxlowT(1:n) = PHILOGx
        end if
        if (iter > 70) run = .false.
        if (P > maxP) maxP = P
        if (ichoice == 2 .and. i > 1) then
            if (P < Pv(i - 1) .and. P < maxP/5 .and. T > 300) run = .false.  ! to finish envelope going to low T bubble
            if (Tcr1 == 0.d0) then
                if (minmaxT) then ! check for self-crossing of the envelope
                    iD = imin
                    do while (Tv(iD) < Tv(i - 1))
                        iD = iD - 1
                    end do
                    if (abs(Pv(iD) - P) < 10) then
                        TpairA = [Tv(iD + 1), Tv(iD)]
                        PpairA = [Pv(iD + 1), Pv(iD)]
                        TpairB = [Tv(i - 1), T]
                        PpairB = [Pv(i - 1), P]
                        call CheckCross(TpairA, PpairA, TpairB, PpairB, Cross, Tcr, Pcr)
                        if (Cross) then
                            Tcr1 = Tcr
                            Pcr1 = Pcr
                            KFcr1(1:n) = XOLD(1:n) + (X(1:n) - XOLD(1:n))*(Tcr - Told)/(T - Told)
                            Kscr1(1:n) = dewK(iD, :n) + (dewK(iD + 1, :n) - dewK(iD, :n))*(Tcr - TpairA(2))/(TpairA(1) - TpairA(2))
                        end if
                        if (.not. cross) then
                            if (Tcr > Tv(i - 1)) then  ! try with the previous step
                                TpairB = [Tv(i - 2), Tv(i - 1)]
                                PpairB = [Pv(i - 2), Pv(i - 1)]
                                call CheckCross(TpairA, PpairA, TpairB, PpairB, Cross, Tcr, Pcr)
                                if (Cross) then
                                    Tcr1 = Tcr
                                    Pcr1 = Pcr
                                    KFcr1(1:n) = XOLD2(1:n) + (XOLD(1:n) - XOLD2(1:n))*(Tcr - Told2)/(Told - Told2)
                             Kscr1(1:n) = dewK(iD, :n) + (dewK(iD + 1, :n) - dewK(iD, :n))*(Tcr - TpairA(2))/(TpairA(1) - TpairA(2))
                                end if
                            end if
                        end if
                    end if
                else if (minT) then  ! check for max after minT, which could indicate that a crossing is coming
                    if (Told > T .and. Told > Told2) minmaxT = .true.
                else    ! check for minT
                    if (Told < T .and. Told < Told2) minT = .true.
                    if (minT) imin = i - 1
                end if
            end if
        else if (ichoice == 1 .and. i > 1) then ! along construction of "bubble" curve from low T...

            if (P < Pv(i - 1) .and. Tcr1 == 0.d0) then ! Max Pressure passed without crossing -->
                !  Start 3envel from low T bub (added June 2017)
                !            run = .false.
            else if (Tcr1 == 0.d0) then ! search for 1st crossing with dew curve coming from high P
                iD = ilastDewC
                ! do while (PdewC(iD)>Pv(i-1))
                do while (PdewC(iD) > P)
                    iD = iD - 1
                end do
                if (abs(TdewC(iD) - T) < 20) then
                    TpairA = [TdewC(iD + 1), TdewC(iD)]
                    PpairA = [PdewC(iD + 1), PdewC(iD)]
                    TpairB = [Tv(i - 1), T]
                    PpairB = [Pv(i - 1), P]
                    call CheckCross(TpairA, PpairA, TpairB, PpairB, Cross, Tcr, Pcr)
                    if (Cross) then
                        Tcr1 = Tcr
                        Pcr1 = Pcr
                        KFcr1(1:n) = XOLD(1:n) + (X(1:n) - XOLD(1:n))*(Tcr - Told)/(T - Told)
                        Kscr1(1:n) = dewK(iD, :n) + (dewK(iD + 1, :n) - dewK(iD, :n))*(Tcr - TpairA(2))/(TpairA(1) - TpairA(2))
                    end if
                end if
                ! do while (Tcr1==0.d0.and.PdewC(iD)>P) !  in case the previous check was not yet for the right segment
                do while (Tcr1 == 0.d0 .and. Tcr > TdewC(iD)) !  in case the previous check was not yet for the right segment
                    iD = iD - 1
                    if (abs(TdewC(iD) - T) < 20) then
                        TpairA = [TdewC(iD + 1), TdewC(iD)]
                        PpairA = [PdewC(iD + 1), PdewC(iD)]
                        TpairB = [Tv(i - 1), T]
                        PpairB = [Pv(i - 1), P]
                        call CheckCross(TpairA, PpairA, TpairB, PpairB, Cross, Tcr, Pcr)
                        if (Cross) then
                            Tcr1 = Tcr
                            Pcr1 = Pcr
                            KFcr1(1:n) = XOLD(1:n) + (X(1:n) - XOLD(1:n))*(Tcr - Told)/(T - Told)
                            Kscr1(1:n) = dewK(iD, :n) + (dewK(iD + 1, :n) - dewK(iD, :n))*(Tcr - TpairA(2))/(TpairA(1) - TpairA(2))
                        end if
                    end if
                end do
            else ! search for 2nd crossing with dew curve coming from high P
                do while (TdewC(iD) < Tv(i - 1))
                    iD = iD - 1
                end do
                if (abs(PdewC(iD) - P) < 20) then
                    TpairA = [TdewC(iD + 1), TdewC(iD)]
                    PpairA = [PdewC(iD + 1), PdewC(iD)]
                    TpairB = [Tv(i - 1), T]
                    PpairB = [Pv(i - 1), P]
                    call CheckCross(TpairA, PpairA, TpairB, PpairB, Cross, Tcr, Pcr)
                    if (Cross) then
                        Tcr2 = Tcr
                        Pcr2 = Pcr
                        KFcr2(1:n) = XOLD(1:n) + (X(1:n) - XOLD(1:n))*(Tcr - Told)/(T - Told)
                        Kscr2(1:n) = dewK(iD, :n) + (dewK(iD + 1, :n) - dewK(iD, :n))*(Tcr - TpairA(2))/(TpairA(1) - TpairA(2))
                        T = Tcr2
                        P = Pcr2
                        run = .false.
                    end if
                end if
                do while (run .and. TdewC(iD) < T)  !  in case the previous check was not yet for the right segment
                    iD = iD - 1
                    if (abs(PdewC(iD) - P) < 20) then
                        TpairA = [TdewC(iD + 1), TdewC(iD)]
                        PpairA = [PdewC(iD + 1), PdewC(iD)]
                        TpairB = [Tv(i - 1), T]
                        PpairB = [Pv(i - 1), P]
                        call CheckCross(TpairA, PpairA, TpairB, PpairB, Cross, Tcr, Pcr)
                        if (Cross) then
                            Tcr2 = Tcr
                            Pcr2 = Pcr
                            KFcr2(1:n) = XOLD(1:n) + (X(1:n) - XOLD(1:n))*(Tcr - Told)/(T - Told)
                            Kscr2(1:n) = dewK(iD, :n) + (dewK(iD + 1, :n) - dewK(iD, :n))*(Tcr - TpairA(2))/(TpairA(1) - TpairA(2))
                            T = Tcr2
                            P = Pcr2
                            run = .false.
                        end if
                    end if
                end do
            end if
        else if (ichoice == 3 .and. P < PmaxDewC) then ! search for crossing with dew curve, coming from high P L-L sat
            iD = ilastDewC
            do while (TdewC(iD) < Tv(i))
                iD = iD - 1
            end do
            if (abs(PdewC(iD) - P) < 30) then
                TpairA = [TdewC(iD + 1), TdewC(iD)]
                PpairA = [PdewC(iD + 1), PdewC(iD)]
                TpairB = [Tv(i - 1), T]
                PpairB = [Pv(i - 1), P]
                call CheckCross(TpairA, PpairA, TpairB, PpairB, Cross, Tcr, Pcr)
                if (Cross) then
                    Tcr1 = Tcr
                    Pcr1 = Pcr
                    Kscr1(1:n) = XOLD(1:n) + (X(1:n) - XOLD(1:n))*(Tcr - Told)/(T - Told)
                    KFcr1(1:n) = dewK(iD, :n) + (dewK(iD + 1, :n) - dewK(iD, :n))*(Tcr - TpairA(2))/(TpairA(1) - TpairA(2))
                    T = Tcr
                    P = Pcr
                    run = .false.
                end if
            end if
            do while (run .and. TdewC(iD) < Tv(i - 1))  !  in case the previous check was not yet for the right segment
                iD = iD - 1
                if (abs(PdewC(iD) - P) < 30) then
                    TpairA = [TdewC(iD + 1), TdewC(iD)]
                    PpairA = [PdewC(iD + 1), PdewC(iD)]
                    TpairB = [Tv(i - 1), T]
                    PpairB = [Pv(i - 1), P]
                    call CheckCross(TpairA, PpairA, TpairB, PpairB, Cross, Tcr, Pcr)
                    if (Cross) then
                        Tcr1 = Tcr
                        Pcr1 = Pcr
                        Kscr1(1:n) = XOLD(1:n) + (X(1:n) - XOLD(1:n))*(Tcr - Told)/(T - Told)
                        KFcr1(1:n) = dewK(iD, :n) + (dewK(iD + 1, :n) - dewK(iD, :n))*(Tcr - TpairA(2))/(TpairA(1) - TpairA(2))
                        T = Tcr
                        P = Pcr
                        run = .false.
                    end if
                end if
            end do
        end if
        print *, T, P, ns, iter
        Tv(i) = T
        Pv(i) = P
        Dv(i) = 1/Vx    ! saturated phase density
        ! rho_y = 1/Vy     incipient phase density
        if (sum(X(:n)*Xold(:n)) < 0) then  ! critical point detected
            ncri = ncri + 1
            icri(ncri) = i - 1
            frac = -Xold(ns)/(X(ns) - Xold(ns))
            Tcri(ncri) = Tv(i - 1) + frac*(T - Tv(i - 1))
            Pcri(ncri) = Pv(i - 1) + frac*(P - Pv(i - 1))
            Dcri(ncri) = Dv(i - 1) + frac*(Dv(i) - Dv(i - 1))
            iy = -iy
            ix = -yx
        end if
        if (run) then
        ! Calculation of sensitivities (dXdS)
        ! dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
            bd = -dFdS
            AJ = JAC
            call dgesv(n + 2, 1, AJ, lda, ipiv, bd, ldb, info)
            if (info .ne. 0) write (6, *) "error with dgesv in parameter ", info
            dXdS = bd
            ! Selection of (the most changing) variable to be specified for the next point
            nsold = ns
            ns = maxloc(abs(dXdS), DIM=1)
            if (maxval(abs(X(:n))) < 0.2) then
                ns = maxloc(abs(dXdS(1:n)), DIM=1)  ! T and P not allowed to be chosen close to a critical point
            end if
            if (ns /= nsold) then
                delS = dXdS(ns)*delS    ! translation of delS to the  new specification variable
                dXdS = dXdS/dXdS(ns)  ! translation of sensitivities
                S = X(ns)               ! update of S
            end if
            ! Setting step in S for the next point to be calculated
            delmax = max(sqrt(abs(X(ns)))/10, 0.1)
            updel = delS*3/iter
            if (passingcri) updel = delS
            if (delS > 0) then
                delS = min(updel, delmax)
            else
                delS = max(updel, -delmax)
            end if
            S = S + delS
            ! Generation of estimates for the next point
            Told2 = Told
            Told = T
            Xold2 = Xold
            Xold = X
            X = Xold + dXdS*delS
            if (passingcri) passingcri = .false.
            do while (maxval(abs(X(:n))) < 0.03)  ! approaching the black hole... get out of there! (0.03)
                stepX = maxval(abs(X(:n) - Xold(:n))) ! the step given by the most changing logK to fall into the black hole
                passingcri = .true.
                if (stepX > 0.07) then
                    S = S - delS/2
                    X = X - dXdS*delS/2   !  half step back
                else
                    S = S + delS
                    X = X + dXdS*delS   ! one more step to jump over the critical point
                end if
            end do
            T = exp(X(n + 1))
            if (.not. passingcri .and. abs(T - Told) > 7) then ! Delta T estimations > 7K are not allowed
                delS = delS/2
                S = S - delS
                X = Xold + dXdS*delS
                T = exp(X(n + 1))
            end if
            P = exp(X(n + 2))
            KFACT = exp(X(:n))
            y = z*KFACT
            if ((dXdS(n + 1)*delS < 0 .and. P < 0.1 .or. T < 120.0) &  ! dew line stops when P<0.1 bar or T<150K
                .or. (P > 1.0 .and. T < 150.0) &   ! bubble line stops when T<150K
                .or. (P > 1500) &
                .or. (abs(dels) < 1.d-4)) then
                run = .false.
            end if
        end if
    end do
    n_points = i
    !-----------------------------------------------------------

    print *, y
    print *, rho_x
    print *, rho_y
    print *, beta
end subroutine envelope2

subroutine envelope3(ichoice, model, n, z, T, P, beta, KFACT, KFsep, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)

    ! Routine for tracing the boundaries of a three-phase region, starting from a known previously detected point.
    ! Developed as an extension and modification of "envelope2". March 2017.

    !   In the main case (ichoice=1) y, xx and w denote molar fractions in the following phases:
    !   y: vapor (incipient phase)
    !  xx: hydrocarbon liquid
    !   w: the other liquid (e.g. aqueous of asphaltenic phase)
    !   beta is the fraction of w phase. Then we have (1-beta) of xx phase and 0 of y (incipient phase)
    !   KFACT = y/xx
    !   KFsep = w/xx

    !   ichoice=2 is used for cases where the incipient phase "y" is the 2nd liquid.
    !   w will correspond to vapor phase, with fraction beta

    !   ichoice=3 is used for cases where the initial saturated phase "xx" is the vapor (e.g. OilB with water from Lindeloff-Michelsen).
    !   y and w will correspond to the two liquid phases, with a beta fraction for w.

    implicit DOUBLE PRECISION(A - H, O - Z)
    PARAMETER(nco=64)

    ! M&M means the book by Michelsen and Mollerup, 2nd Edition (2007)

    ! eos id, number of compounds in the system and starting point type
    integer, intent(in) :: model, n, ichoice

    DOUBLE PRECISION Kinf

    ! estimated T, P and "w" phase fraction for first point (then used for every point)
    DOUBLE PRECISION :: T, P, beta

    ! estimated K factors for first point (then used for every point)
    DOUBLE PRECISION, dimension(n) :: KFACT, KFsep

    ! composition of the system
    DOUBLE PRECISION, dimension(n), intent(in) :: z

    ! pure compound physical constants
    DOUBLE PRECISION, dimension(n), intent(in) :: tcn
    DOUBLE PRECISION, dimension(n), intent(in) :: pcn
    DOUBLE PRECISION, dimension(n), intent(in) :: omgn

    ! eos parameters
    DOUBLE PRECISION, dimension(n), intent(in) :: acn  ! in bar*(L/mol)**2
    DOUBLE PRECISION, dimension(n), intent(in) :: bn   ! in L/mol
    DOUBLE PRECISION, dimension(n), intent(in) :: delta1n  !only required for RKPR
    DOUBLE PRECISION, dimension(n), intent(in) :: k_or_mn  ! k for RKPR ; m for SRK/PR

    ! interaction parameters matrices
    DOUBLE PRECISION, dimension(n, n), intent(in) :: Kij_or_K0n
    DOUBLE PRECISION, dimension(n, n), intent(in) :: Tstarn
    DOUBLE PRECISION, dimension(n, n), intent(in) :: Lijn

    ! T, P and Density of the calculated envelope
    DOUBLE PRECISION, dimension(800), intent(out) :: Tv
    DOUBLE PRECISION, dimension(800), intent(out) :: Pv
    DOUBLE PRECISION, dimension(800), intent(out) :: Dv

    ! number of valid elements in Tv, Pv and Dv arrays
    integer, intent(out) :: n_points

    ! positions of the last saturation points before each critical point
    integer, dimension(4), intent(out) :: icri
    ! T, P and Density of critical points
    DOUBLE PRECISION, dimension(4), intent(out) :: Tcri
    DOUBLE PRECISION, dimension(4), intent(out) :: Pcri
    DOUBLE PRECISION, dimension(4), intent(out) :: Dcri

    ! number of valid elements in icri, Tcri, Pcri and Dcri arrays
    integer, intent(out) :: ncri

    ! Intermediate variables during calculation process
    DOUBLE PRECISION, dimension(n) :: y, xx, w, PHILOGy, PHILOGx, PHILOGw
    DOUBLE PRECISION, dimension(n) :: dxdB, dydB, dwdB, dxdKs, dydKs, dwdKs, aux
    DOUBLE PRECISION, dimension(n) :: DLPHITx, DLPHIPx, DLPHITy, DLPHIPy, DLPHITw, DLPHIPw
    DOUBLE PRECISION, dimension(n, n) :: FUGNx, FUGNy, FUGNw
    integer, dimension(2*n + 3) :: ipiv
    DOUBLE PRECISION, dimension(2*n + 3) :: X, Xold, delX, bd, F, dFdS, dXdS
    DOUBLE PRECISION, dimension(2*n + 3) :: Fp, JACnumK2, JACnumKn, JACnumKs2, JACnumKsn, JACnumlT, JACnumlP, JACnumB
    DOUBLE PRECISION, dimension(2*n + 3, 2*n + 3) :: JAC, AJ
    DOUBLE PRECISION :: Vy, Vx, Vw
    ! DOUBLE PRECISION, dimension(2) :: TpairA, PpairA, TpairB, PpairB    ! crossing vars
    logical :: run, passingcri, Comp3ph  !, Cross                       ! crossing var (Cross)

    DOUBLE PRECISION, dimension(nco, nco) :: Kij_or_K0, Tstar
    DOUBLE PRECISION, dimension(nco) :: KFsep1
    COMMON/CRIT/TC(nco), PC(nco), DCeos(nco), omg(nco)
    COMMON/COMPONENTS/ac(nco), b(nco), delta1(nco), rk_or_m(nco), Kij_or_K0, NTDEP
    COMMON/MODEL/NMODEL
    COMMON/rule/ncomb
    COMMON/bcross/bij(nco, nco)
    COMMON/Tdep/Kinf, Tstar
    COMMON/lowTKsep/KFsep1
    COMMON/writeComp/Comp3ph, i1, i2
    ! COMMON /DewCurve/ ilastDewC, TdewC(500), PdewC(500)     ! crossing vars
    ! COMMON /CrossingPoints/ Tcr1,Pcr1,Tcr2,Pcr2,KFcr1,Kscr1,KFcr2,Kscr2

    ! Charging the commons(nco) from input arguments (n)
    NMODEL = model
    TC(:n) = tcn
    PC(:n) = pcn
    OMG(:n) = omgn
    ac(:n) = acn
    b(:n) = bn
    delta1(:n) = delta1n
    rk_or_m(:n) = k_or_mn
    Kij_or_K0(:n, :n) = Kij_or_K0n
    Kinf = 0.0d0
    ncomb = 0  ! only  vdW combining rules and quadratic mixing rules by  the moment
    Tstar(:n, :n) = Tstarn
    ! b matrix for Classical or van der Waals combining rules:
    do i = 1, n
        do j = i, n
            bij(i, j) = (1 - lijn(i, j))*(b(i) + b(j))/2
            bij(j, i) = bij(i, j)
        end do
    end do
    !
    !-----------------------------------------------------------

    ! Continuation method for tracing the envelope starts here
    run = .true.
    i = 0
    ncri = 0
    JAC = 0.d0
    lda = 2*n + 3
    ldb = 2*n + 3
    X(:n) = log(KFACT)
    X(n + 1:2*n) = log(KFsep)
    X(2*n + 1) = log(T)
    X(2*n + 2) = log(P)
    X(2*n + 3) = beta
    !
    ! To test Jacobian insert here block for numerical derivatives at the end of this code
    !
    iy = 1
    ix = 1
    iw = 1
    if (ichoice == 1) iy = -1
    if (ichoice == 2) iw = -1  ! w will be vapor phase during the first part
    if (ichoice == 3) ix = -1  ! x will be vapor phase, saturated in the first point
    if (beta == 0.0d0) then
        ns = 2*n + 3
        S = 0.00   ! for beta
        delS = 0.001
    else ! start from low T "bubble" point or low T "Lower AOP"
        ns = 2*n + 1
        S = X(2*n + 1)   ! for log(T)
        delS = 0.001
    end if
    xx = z/(1 - beta + beta*KFsep)
    y = KFACT*xx
    w = KFsep*xx
    Xold = 0.d0
    dFdS = 0.d0
    dFdS(2*n + 3) = -1.d0
    write (2, *) '     T       P      beta     X(1)     X(n)     X(n+1)     X(2*n)    ns  iter'
    if (Comp3ph) write (3, *) '     T       P      beta     xa     xb     ya     yb    wa     wb'
    do while (run)
        i = i + 1  ! number of point to be calculated along the line
        ! Newton starts here
        delX = 1.0
        iter = 0
        do while (maxval(abs(delX)) > 1.d-5 .and. iter <= 70)
            iter = iter + 1
            !      nc,MTYP,INDIC,T,P,rn,V,PHILOG,DLPHI,DLPHIP,DLPHIT,FUGN
            call TERMO(n, iy, 4, T, P, y, Vy, PHILOGy, DLPHIPy, DLPHITy, FUGNy)
            call TERMO(n, ix, 4, T, P, xx, Vx, PHILOGx, DLPHIPx, DLPHITx, FUGNx)
            call TERMO(n, iw, 4, T, P, w, Vw, PHILOGw, DLPHIPw, DLPHITw, FUGNw)
            F(:n) = X(:n) + PHILOGy - PHILOGx  ! X(:n) are LOG_K
            F(n + 1:2*n) = X(n + 1:2*n) + PHILOGw - PHILOGx  ! X(:n) are LOG_K
            F(2*n + 1) = sum(y - xx)
            F(2*n + 2) = sum(w - xx)
            F(2*n + 3) = X(ns) - S
            dxdB = -(KFsep - 1.d0)*xx*xx/z
            dydB = KFACT*dxdB
            dwdB = KFsep*dxdB
            aux = -beta*xx/z
            dxdKs = aux*xx
            dydKs = aux*y
            dwdKs = xx*(1 + aux*KFsep)
            do j = 1, n
                JAC(1:n, j) = FUGNy(:, j)*y(j)  ! y=K*xx
                JAC(j, j) = JAC(j, j) + 1.d0
            end do
            do j = n + 1, 2*n    ! wrt Ks
                JAC(1:n, j) = KFsep(j - n)*(FUGNy(:, j - n)*dydKs(j - n) - FUGNx(:, j - n)*dxdKs(j - n))
            end do
            JAC(1:n, 2*n + 1) = T*(DLPHITy - DLPHITx)  ! wrt T
            JAC(1:n, 2*n + 2) = P*(DLPHIPy - DLPHIPx)  ! wrt P
            do l = 1, n    ! wrt beta
                JAC(l, 2*n + 3) = sum(FUGNy(l, :)*dydB - FUGNx(l, :)*dxdB)
            end do
            ! ders of F(n+1:2*n) wrt K = 0
            do j = n + 1, 2*n  ! wrt Ks
                JAC(n + 1:2*n, j) = KFsep(j - n)*(FUGNw(:, j - n)*dwdKs(j - n) - FUGNx(:, j - n)*dxdKs(j - n))
                JAC(j, j) = JAC(j, j) + 1.d0
            end do
            JAC(n + 1:2*n, 2*n + 1) = T*(DLPHITw - DLPHITx)  ! wrt T
            JAC(n + 1:2*n, 2*n + 2) = P*(DLPHIPw - DLPHIPx)  ! wrt P
            do l = n + 1, 2*n    ! wrt beta
                JAC(l, 2*n + 3) = sum(FUGNw(l - n, :)*dwdB - FUGNx(l - n, :)*dxdB)
            end do
            JAC(2*n + 1, 1:n) = y             ! sum(y-x) wrt K
            JAC(2*n + 1, n + 1:2*n) = KFsep*(dydKs - dxdKs)    ! sum(y-x) wrt Ks
            ! sum(y-x) wrt T or P = 0
            JAC(2*n + 1, 2*n + 3) = sum(dydB - dxdB)  ! wrt beta
            ! sum(w-x) wrt K = 0
            JAC(2*n + 2, n + 1:2*n) = KFsep*(dwdKs - dxdKs)    ! sum(w-x) wrt Ks
            ! sum(w-x) wrt T or P = 0
            JAC(2*n + 2, 2*n + 3) = sum(dwdB - dxdB)  ! wrt beta
            JAC(2*n + 3, :) = 0.d0
            JAC(2*n + 3, ns) = 1.d0
            !After JAC calculated, print to compare with numerical values:
            ! write(2,*)'JACK2:  ',JAC(:,2)
            ! write(2,*)'JACKn:  ',JAC(:,n)
            ! write(2,*)'JACKs2: ',JAC(:,n+2)
            ! write(2,*)'JACKsn: ',JAC(:,2*n)
            ! write(2,*)'JAClT:  ',JAC(:,2*n+1)
            ! write(2,*)'JAClP:  ',JAC(:,2*n+2)
            ! write(2,*)'JACB:   ',JAC(:,2*n+3)

            !call dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
            bd = -F
            AJ = JAC
            call dgesv(2*n + 3, 1, AJ, lda, ipiv, bd, ldb, info)
            if (info .ne. 0) write (6, *) "error with dgesv in parameter ", info
            delX = bd
            if (i == 1) then
                do while (maxval(abs(delX)) > 1.0)   ! Too large Newton step --> Reduce it
                    delX = delX/2
                end do
            else
                do while (maxval(abs(delX)) > 0.08)   ! Too large Newton step --> Reduce it
                    delX = delX/2
                end do
                if (iter > 10) delX = delX/2  ! too many iterations (sometimes due to oscillatory behavior near crit) --> Reduce it
            end if
            X = X + delX
            if (.not. passingcri .and. i /= 1 .and. iter == 20 .and. maxval(abs(delX)) > 0.001) then ! Too many iterations-->Reduce step to new point
                delS = delS/2
                S = S - delS
                X = Xold + dXdS*delS
                iter = 0
            end if
            KFACT = exp(X(:n))
            KFsep = exp(X(n + 1:2*n))
            T = exp(X(2*n + 1))
            P = exp(X(2*n + 2))
            beta = X(2*n + 3)
            xx = z/(1 - beta + beta*KFsep)
            y = KFACT*xx
            w = KFsep*xx
        end do
        ! Point converged (unless it jumped out because of high number of iterations)
        if (iter > 70) run = .false.
        if (beta < 0) run = .false.
        if (ichoice == 1 .and. i == 1) KFsep1(1:n) = KFsep
        print *, T, P, ns, iter
        write (2, 1) T, P, beta, X(1), X(n), X(n + 1), X(2*n), ns, iter
        if (Comp3ph) write (3, 3) T, P, beta, xx(i1), xx(i2), y(i1), y(i2), w(i1), w(i2)
        Tv(i) = T
        Pv(i) = P
        Dv(i) = 1/Vx    ! saturated phase density
        ! rho_y = 1/Vy     incipient phase density
        if (sum(X(:n)*Xold(:n)) < 0) then  ! critical point detected between x and y phases
            ncri = ncri + 1
            icri(ncri) = i - 1
            frac = -Xold(ns)/(X(ns) - Xold(ns))
            Tcri(ncri) = Tv(i - 1) + frac*(T - Tv(i - 1))
            Pcri(ncri) = Pv(i - 1) + frac*(P - Pv(i - 1))
            Dcri(ncri) = Dv(i - 1) + frac*(Dv(i) - Dv(i - 1))
            iy = -iy
            ix = -yx
        end if
        if (sum(X(n + 1:2*n)*Xold(n + 1:2*n)) < 0) then  ! critical point detected between x and w phases
            ncri = ncri + 1
            icri(ncri) = i - 1
            frac = -Xold(ns)/(X(ns) - Xold(ns))
            Tcri(ncri) = Tv(i - 1) + frac*(T - Tv(i - 1))
            Pcri(ncri) = Pv(i - 1) + frac*(P - Pv(i - 1))
            Dcri(ncri) = Dv(i - 1) + frac*(Dv(i) - Dv(i - 1))
            iw = -iw
            ix = -yx
        end if
        if (run) then
            ! Calculation of sensitivities (dXdS)
            ! dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
            bd = -dFdS
            AJ = JAC
            call dgesv(2*n + 3, 1, AJ, lda, ipiv, bd, ldb, info)
            if (info .ne. 0) write (6, *) "error with dgesv in parameter ", info
            dXdS = bd
            ! Selection of (the most changing) variable to be specified for the next point
            nsold = ns
            ns = maxloc(abs(dXdS), DIM=1)
            if (maxval(abs(X(:n))) < 0.2) then
                ns = maxloc(abs(dXdS(1:n)), DIM=1)  ! vars other than logK not allowed to be chosen close to a y-x critical point
            end if
            if (maxval(abs(X(n + 1:2*n))) < 0.2) then
                ns = maxloc(abs(dXdS(n + 1:2*n)), DIM=1)  ! vars other than logKs not allowed to be chosen close to a w-x critical point
            end if
            if (ns /= nsold) then
                delS = dXdS(ns)*delS    ! translation of delS to the  new specification variable
                dXdS = dXdS/dXdS(ns)  ! translation of sensitivities
                S = X(ns)               ! update of S
            end if
            ! Setting step in S for the next point to be calculated
            delmax = max(sqrt(abs(X(ns)))/10, 0.1)
            updel = delS*3/iter
            if (passingcri) updel = delS
            if (delS > 0) then
                delS = min(updel, delmax)
            else
                delS = max(updel, -delmax)
            end if
            S = S + delS
            ! Generation of estimates for the next point
            Told = T
            Xold = X
            X = Xold + dXdS*delS
            if (passingcri) passingcri = .false.
            ! do while (maxval(abs(X(:n)))<0.03)  ! approaching the black hole... get out of there! (0.03)
            !   stepX = maxval(abs(X(:n)-Xold(:n))) ! the step given by the most changing logK to fall into the black hole
            !   passingcri = .true.
            !   if (stepX > 0.07) then
            !       S = S - delS/2
            !       X = X - dXdS * delS/2   !  half step back
            !   else
            !       S = S + delS
            !       X = X + dXdS * delS   ! one more step to jump over the critical point
            !   end if
            ! end do
            if (maxval(abs(X(:n))) < 0.03) then  ! approaching the black hole... get out of there!
                if (delS > 0) delS = 0.04 - S
                if (delS < 0) delS = -0.04 - S
                S = S + delS
                X = X + dXdS*delS   ! one more step to jump over the critical point
                passingcri = .true.
            end if
            do while (maxval(abs(X(n + 1:2*n))) < 0.03)  ! approaching the black hole... get out of there! (0.03)
                stepX = maxval(abs(X(n + 1:2*n) - Xold(n + 1:2*n))) ! the step given by the most changing logKs to fall into the black hole
                passingcri = .true.
                if (stepX > 0.07) then
                    S = S - delS/2
                    X = X - dXdS*delS/2   !  half step back
                else
                    S = S + delS
                    X = X + dXdS*delS   ! one more step to jump over the critical point
                end if
            end do
            T = exp(X(2*n + 1))
            if (.not. passingcri .and. abs(T - Told) > 7) then ! Delta T estimations > 7K are not allowed
                delS = delS/2
                S = S - delS
                X = Xold + dXdS*delS
                T = exp(X(2*n + 1))
            end if
            P = exp(X(2*n + 2))
            KFACT = exp(X(:n))
            KFsep = exp(X(n + 1:2*n))
            beta = X(2*n + 3)
            xx = z/(1 - beta + beta*KFsep)
            y = KFACT*xx
            w = KFsep*xx
            if ((dXdS(2*n + 1)*delS < 0 .and. P < 0.1 .or. T < 120.0) &  ! this may need adjustment
                .or. (P > 1.0 .and. T < 150.0) &
                .or. (P > 1500) &
                .or. (abs(dels) < 1.d-4)) then
                run = .false.
            end if
        end if
    end do
    n_points = i
    !-----------------------------------------------------------

    print *, y
    print *, rho_x
    print *, rho_y
    print *, beta
1   FORMAT(7F10.4, 2I4)
3   FORMAT(3F10.4, 6E12.3)
end subroutine envelope3

subroutine EvalFEnvel3(n, z, X, F)
    implicit DOUBLE PRECISION(A - H, O - Z)
    DOUBLE PRECISION, dimension(n) :: KFACT, KFsep
    DOUBLE PRECISION, dimension(n) :: z, y, xx, w, PHILOGy, PHILOGx, PHILOGw
    DOUBLE PRECISION, dimension(n) :: DLPHITx, DLPHIPx, DLPHITy, DLPHIPy, DLPHITw, DLPHIPw
    DOUBLE PRECISION, dimension(n, n) :: FUGNx, FUGNy, FUGNw
    DOUBLE PRECISION, dimension(2*n + 3) :: X, F

    S = 0.001
    KFACT = exp(X(:n))
    KFsep = exp(X(n + 1:2*n))
    T = exp(X(2*n + 1))
    P = exp(X(2*n + 2))
    beta = X(2*n + 3)
    xx = z/(1 - beta + beta*KFsep)
    y = KFACT*xx
    w = KFsep*xx
    call TERMO(n, 1, 4, T, P, y, Vy, PHILOGy, DLPHIPy, DLPHITy, FUGNy)
    call TERMO(n, -1, 4, T, P, xx, Vx, PHILOGx, DLPHIPx, DLPHITx, FUGNx)
    call TERMO(n, 1, 4, T, P, w, Vw, PHILOGw, DLPHIPw, DLPHITw, FUGNw)
    F(:n) = X(:n) + PHILOGy - PHILOGx  ! X(:n) are LOG_K
    F(n + 1:2*n) = X(n + 1:2*n) + PHILOGw - PHILOGx  ! X(:n) are LOG_K
    F(2*n + 1) = sum(y - xx)
    F(2*n + 2) = sum(w - xx)
    F(2*n + 3) = X(25) - S
end subroutine EvalFEnvel3
