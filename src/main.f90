! Listado de commons utilizados, que podrían pasar a un módulo:
!   COMMON /CRIT/TC(nco),PC(nco),DCeos(nco),omg(nco)
!        COMMON /COMPONENTS/ ac(nco),b(nco),delta1(nco),rk(nco),Kij_or_K0,NTDEP
!        COMMON /MODEL/ NMODEL
!        COMMON /rule/ncomb
!        COMMON /bcross/bij(nco,nco)
!        COMMON /Tdep/ Kinf,Tstar

program calc_envelope2and3
    implicit DOUBLE PRECISION(A - H, O - Z)
    LOGICAL Comp3ph
    COMMON/writeComp/Comp3ph, i1, i2
    OPEN (1, FILE='envelIN.txt')
    OPEN (2, FILE='envelOUT.txt')
    read (1, *) N
    !write (6, *) 'write extra output with compositions for 2 compounds along 3-phase lines?'
    !write (6, *) 'Enter 1 for YES. Otherwise, any other number.'
    ! read (5, *) i
    i = 0
    if (i == 1) Comp3ph = .true.
    if (Comp3ph) then
        OPEN (3, FILE='Comp3phOUT.txt')
        write (6, *) 'Enter order numbers for two selected compounds, separated by space.'
        read (5, *) i1, i2
        write (3, *) 'Molar fractions along three-phase boundaries are printed below for compounds with order:', i1, i2
    end if

    call readcase(n)
end program
