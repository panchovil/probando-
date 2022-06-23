subroutine readcase(n)
   implicit DOUBLE PRECISION(A - H, O - Z)
   PARAMETER(nco=64, Pmax=700.0)
   dimension z(n), xx(n), w(n)

   DOUBLE PRECISION Kinf
   DOUBLE PRECISION, dimension(n) :: Kfact, KFsep
   DOUBLE PRECISION, dimension(nco) :: KFcr1, Kscr1, KFcr2, Kscr2, KlowT, PHILOGxlowT, KFsep1 ! go in COMMONS (cannot have n dimension)
   ! pure compound physical constants
   DOUBLE PRECISION, dimension(n) :: tcn
   DOUBLE PRECISION, dimension(n) :: pcn
   DOUBLE PRECISION, dimension(n) :: omgn

   ! eos parameters
   DOUBLE PRECISION, dimension(n) :: acn  ! in bar*(L/mol)**2
   DOUBLE PRECISION, dimension(n) :: bn   ! in L/mol
   DOUBLE PRECISION, dimension(n) :: delta1n  !only required for RKPR
   DOUBLE PRECISION, dimension(n) :: k_or_mn  ! k for RKPR ; m for SRK/PR

   ! interaction parameters matrices
   DOUBLE PRECISION, dimension(n, n) :: Kij_or_K0n, Lijn
   DOUBLE PRECISION, dimension(n, n) :: Tstarn

   ! interaction parameters matrices
   DOUBLE PRECISION, dimension(nco, nco) :: Kij_or_K0, Lij, Tstar

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

   DOUBLE PRECISION, dimension(n) :: y, PHILOGy, PHILOGx
   DOUBLE PRECISION, dimension(n) :: DLPHITx, DLPHIPx, DLPHITy, DLPHIPy
   DOUBLE PRECISION, dimension(n, n) :: FUGNx, FUGNy

   CHARACTER*4 :: spec
   LOGICAL FIRST

   COMMON/MODEL/NMODEL
   COMMON/CRIT/TC(nco), PC(nco), DCeos(nco), omg(nco)
   COMMON/COMPONENTS/ac(nco), b(nco), delta1(nco), rk_or_m(nco), Kij_or_K0, NTDEP
   COMMON/rule/ncomb
   COMMON/bcross/bij(nco, nco)
   COMMON/Tdep/Kinf, Tstar
   COMMON/lforin/lij
   COMMON/DewCurve/ilastDewC, TdewC(800), PdewC(800), dewK(800, nco)
   COMMON/CrossingPoints/Tcr1, Pcr1, Tcr2, Pcr2, KFcr1, Kscr1, KFcr2, Kscr2
   COMMON/lowTbub/TlowT, PlowT, KlowT, PHILOGxlowT !shared with envelope2
   COMMON/lowTKsep/KFsep1     !shared with envelope3

   Tcr1 = 0.d0 ! T of 1st crossing point detected between different envelope segments
   Tcr2 = 0.d0
   READ (1, *) (z(j), j=1, N)
   read (1, *) nmodel
   if (nmodel < 3) then
      call read2PcubicNC(N, 1, 2)
   else if (nmodel == 3) then
      call readRKPRNC(N, 1, 2)
   end if
   WRITE (2, *)
   write (2, 4) (z(i), i=1, n)
4  FORMAT('Molar fractions: ', 20F7.4)

   ! Passing values from commons(nco) to input arguments (n)
   TCn = tc(:n)
   PCn = pc(:n)
   OMGn = omg(:n)
   acn = ac(:n)
   bn = b(:n)
   delta1n = delta1(:n)
   k_or_mn = rk_or_m(:n)
   Kij_or_K0n = Kij_or_K0(:n, :n)
   Tstarn = Tstar(:n, :n)
   lijn = lij(:n, :n)

   ! Start from dew curve at low T and P
   ! Although the algorithm is written to start from a bubble point, here we invert the Wilson K factors
   ! and y will actually mean x (ichoice = 2)

   ichoice = 2
   P = 1.0
   T = 315.0
   do while (P > 0.1)
      T = T - 5.D0
      P = 1.d0/sum(z/(PCn*EXP(5.373*(1 + omgn)*(1 - TCn/T))))
   end do
   KFACT = PCn*EXP(5.373*(1 + omgn)*(1 - TCn/T))/P  ! standard Wilson K factors
   KFACT = 1.d0/KFACT  ! inversion

   call envelope2(ichoice, nmodel, n, z, T, P, KFACT, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                  Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
   call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
   TdewC = Tv
   PdewC = Pv
   ilastDewC = n_points

   IF (Tcr1 > 0.d0) then  ! self crossing detected (usual in some asymmetric hc mixtures)

   ELSE IF (P > Pmax) then  ! now run from Low T Bubble point
      ichoice = 1
      P = 11.0   ! 11.0
      T = 205.0
      do while (P > 10) ! > 10
         T = T - 5.D0
         P = sum(z*PCn*EXP(5.373*(1 + omgn)*(1 - TCn/T)))
      end do
      KFACT = PCn*EXP(5.373*(1 + omgn)*(1 - TCn/T))/P
      call envelope2(ichoice, nmodel, n, z, T, P, KFACT, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
   ELSE  ! now run from High P L-L saturation (incipient phase rich in last comp.)
      ichoice = 3
      P = Pmax
      T = 710.0
      iy = 1
      ix = 1
      y = 0.d0
      y(n) = 1.d0
      difw = -1.d0
      !        WRITE (2,*) '  P = ', P
      !        WRITE (2,*) '     T    dif LogFug'
      do while (difW < 0.d0)
         T = T - 10.D0
         call TERMO(n, ix, 1, T, P, z, Vx, PHILOGx, DLPHIPx, DLPHITx, FUGNx) ! for fluid
         call TERMO(n, iy, 1, T, P, y, Vy, PHILOGy, DLPHIPy, DLPHITy, FUGNy) ! for pure Water/Asphaltene
         difW = log(z(n)) + PHILOGx(n) - log(y(n)) - PHILOGy(n)
         !            WRITE (2,*) T, difW
      end do
      KFACT = 1.D-3
      KFACT(n) = 1/z(n)
      call envelope2(ichoice, nmodel, n, z, T, P, KFACT, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
   END IF

   print *, Tcr1, Pcr1, "cross1"
   print *, Tcr2, Pcr2, "cross2"

   if (Tcr1 == 0.0d0) then ! no crossings
      print *, "no cross"
      ichoice = 1   ! first inner curve: 3-phase bubble curve (incipient: V)
      T = TlowT
      P = PlowT
      beta = z(n) ! asphaltenes (or last, separating compound) molar fraction
      KFACT = KlowT(1:n)
      y = 0.d0
      y(n) = 1.d0
      call TERMO(n, 1, 1, T, P, y, Vy, PHILOGy, DLPHIPy, DLPHITy, FUGNy)
      KFsep = exp(PHILOGxlowT(1:n) - PHILOGy)
      call envelope3(ichoice, nmodel, n, z, T, P, beta, KFACT, KFsep, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)

      ichoice = 2   ! second inner curve: Lower AOP curve (incipient: Asph. or 2nd L)
      T = TlowT
      ! beta = sum(z(1:9)) ! volatile compounds until C5  (GENERALIZE LATER as sum of all compounds with Tc<470)
      ! s = 2
      ! do while (s>1.2)
      !     P = P/2
      !     beta = 0.8 * beta
      !     KFsep = PCn * EXP(5.373*(1+omgn)*(1-TCn/T)) / P  ! standard Wilson K factors
      !     xx = z/(1-beta+beta*KFsep)
      !     s = sum(xx)
      ! end do
      ! KFACT(n) = 1/xx(n)
      ! P = P1 + (P2-P1) * (1-S1) / (S2-S1)
      ! KFACT = KFsep1(1:n) ! from converged first point for 3-phase bubble curve
      P = PlowT/2
      FIRST = .TRUE.
      spec = 'TP'
      call flash(spec, FIRST, nmodel, n, z, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                 Kij_or_K0n, Tstarn, Lijn, t, p, v, xx, w, rho_x, rho_y, beta, iter)
      call TERMO(n, 1, 1, T, P, xx, Vx, PHILOGx, DLPHIPy, DLPHITy, FUGNy)
      y = 0.d0
      y(n) = 1.d0
      call TERMO(n, 1, 1, T, P, y, Vy, PHILOGy, DLPHIPy, DLPHITy, FUGNy)
      dif = PHILOGy(n) - (PHILOGx(n) + log(xx(n)))
      Pold = P
      if (dif < 0) then
         P = 0.9*P
      else
         P = 1.1*P
      end if
      ! yn = xx(n)*KFACT(n)
      ! do while (yn>2.0)
      do while (abs(dif) > 0.1 .and. P > 0.9)
         call flash(spec, FIRST, nmodel, n, z, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                    Kij_or_K0n, Tstarn, Lijn, t, p, v, xx, w, rho_x, rho_y, beta, iter)
         call TERMO(n, 1, 1, T, P, xx, Vx, PHILOGx, DLPHIPy, DLPHITy, FUGNy)
         call TERMO(n, 1, 1, T, P, y, Vy, PHILOGy, DLPHIPy, DLPHITy, FUGNy)
         dold = dif
         dif = PHILOGy(n) - (PHILOGx(n) + log(xx(n)))
         aux = P
         P = max(P/10, Pold - dold*(P - Pold)/(dif - dold))
         Pold = aux
      end do
      if (abs(dif) > 0.1) then
         FIRST = .TRUE.
         Told = T
         T = T + 10.0
      end if
      do while (abs(dif) > 0.1)
         call flash(spec, FIRST, nmodel, n, z, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                    Kij_or_K0n, Tstarn, Lijn, t, p, v, xx, w, rho_x, rho_y, beta, iter)
         call TERMO(n, 1, 1, T, P, xx, Vx, PHILOGx, DLPHIPy, DLPHITy, FUGNy)
         call TERMO(n, 1, 1, T, P, y, Vy, PHILOGy, DLPHIPy, DLPHITy, FUGNy)
         dold = dif
         dif = PHILOGy(n) - (PHILOGx(n) + log(xx(n)))
         aux = T
         T = min(Told - dold*(T - Told)/(dif - dold), T + 20.0)
         Told = aux
      end do
      KFACT = exp(PHILOGx - PHILOGy)
      KFsep = w/xx
      call envelope3(ichoice, nmodel, n, z, T, P, beta, KFACT, KFsep, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)

   else  ! at least one crossing
      if (ichoice /= 3) ichoice = 1
      print *, "at least one cross", Tcr1, Pcr1
      T = Tcr1
      P = Pcr1
      beta = 0.0d0
      KFACT = exp(KFcr1)
      KFsep = exp(Kscr1)
      call envelope3(ichoice, nmodel, n, z, T, P, beta, KFACT, KFsep, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)

      if (ichoice /= 3) ichoice = 2
      T = Tcr1
      P = Pcr1
      beta = 0.0d0
      KFACT = exp(Kscr1) ! now y (incipient phase in envelope3) will be the second liquid
      KFsep = exp(KFcr1) ! w will be vapor, with phase fraction beta
      call envelope3(ichoice, nmodel, n, z, T, P, beta, KFACT, KFsep, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
   end if

   if (Tcr2 > 0) then
      print *, "second cross?", Tcr2, Pcr2
      ichoice = 1
      T = Tcr2
      P = Pcr2
      beta = 0.0d0
      KFACT = exp(KFcr2)
      KFsep = exp(Kscr2)
      call envelope3(ichoice, nmodel, n, z, T, P, beta, KFACT, KFsep, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)

      ichoice = 2
      T = Tcr2
      P = Pcr2
      beta = 0.0d0
      KFACT = exp(Kscr2) ! now y (incipient phase in envelope3) will be the second liquid
      KFsep = exp(KFcr2) ! w will be vapor, with phase fraction beta
      call envelope3(ichoice, nmodel, n, z, T, P, beta, KFACT, KFsep, tcn, pcn, omgn, acn, bn, k_or_mn, delta1n, &
                     Kij_or_K0n, Tstarn, Lijn, n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
      call WriteEnvel(n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri)
   end if

end subroutine readcase
