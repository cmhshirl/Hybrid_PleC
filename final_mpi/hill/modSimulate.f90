module simulation
  use helper
  use stochastic
  use solver
  implicit none

  logical:: flag
  double precision, parameter :: TD = 1.0
  
contains

  subroutine simulate(T1)
    double precision, intent(in) :: T1
    double precision :: ts, te, random, propensity
    integer :: i

    ! --- Init State ----
    flag = .true.

    !/* Reinitialize */
    a0 = 0.0d0;
    do i=1, RNUM
      call cal_propensity_ssa(i, propensity)
      a(i) = propensity
!print *, 'i: ', i, 'a(i): ', a(i), a0
      a0 = a0 + a(i)
    end do

    call random_number(random)
    y(neq) = log(random)
    ts = 0.0d0
    te = 1.0d0

    ! ---- Simulate ----
    do while (.true.)  
       istate = 1
       jroot  = 0        
       call integrate(ts, te, istate, jroot)

       if (istate == 3) then
!print *, 'te: ', te, 'a0:', a0 !, 'a(68):',a(68),'p(14):',p(14),'rate(68):',RATE(68)
             call propensity_update_ode()
             call ssa_alg()
             call random_number(random)
             y(neq) = log(random)

!stop -1
       elseif (istate == 2) then
             call write_file(te)
             te = te + TD
             if (te >= T1) then
                exit
             end if
!print *,  'CckA_phos:',p(17), 'CckA_kin:',p(18), 'divl:',p(15)
!print *, 'hill:',a(31),'rev_hill:',a(65)

!print *,'te:',te,'a0:', a0!, p(6) !,'a(7,8):', a(7),a(8)
!do i=1,RNUM
!    print *, 'a(', i, '):', a(i)
!end do

        else
            flag = .false.
            print *, 'istate = ', istate
            stop -1
        end if

    end do ! finish a cell simulation
    
    print *,'te:',te,'a0:', a0
    call write_file(te)
    call write_fire()

    print*, 'flag = ', flag

  end subroutine simulate


end module simulation

