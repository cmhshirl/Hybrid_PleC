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
      call cal_propensity(i, propensity)
      a(i) = propensity
!print *, 'i: ', i, 'a(i): ', a(i)
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
             call propensity_update_ode()
             call ssa_alg()
             call random_number(random)
             y(neq) = log(random)
!print *, 'a0: ', te
!stop -1
       elseif (istate == 2) then
             call write_file(te)
             te = te + TD
             if (te >= T1) then
                exit
             end if
       else
          flag = .false.
          print *, 'istate = ', istate
          stop -1
       end if

    end do ! finish a cell simulation

    print*, 'flag = ', flag

  end subroutine simulate


end module simulation

