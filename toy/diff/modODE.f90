# define IJth(y,i,j) y((i-1)*50+j)

module deterministic
  use parameters
  use rand
  integer, parameter:: sys_size = SNUM_ODE*MM + 2
  integer, parameter:: neq = sys_size
  double precision, dimension(neq) :: y

contains

  !****************************
  !*                          *
  !* Initialize ode           *
  !*                          *
  !****************************
  subroutine init_ode()
    implicit none
    integer :: i, j
    double precision:: random

    do i = 1, SNUM_ODE
        do j = 1, MM
            IJth(y,i,j) = 0.0
        end do
    end do

    y(neq-1) = 1.3/MM
    h = y(neq-1)
    call random_number(random)
    y(neq) = log(random)

  end subroutine init_ode


  !****************************
  !*                          *
  !* ODE System of Cell Cycle *
  !*                          *
  !****************************

  subroutine ode_func( t, y, dy)
    implicit none
    integer :: i, j, ileft, iright
    double precision :: c, clt, crt, diff, TWO=2.0
    double precision, intent(in) ::t, y(neq)
    double precision, intent(out)::dy(neq)


    call propensity_update_ode()

    !print *, 't:', t, 'a0:', a0, 'y(neq):', y(neq)

    do j = 1, MM
      if (j==1) then
        ileft = 0
      else
        ileft = -1
      end if

      if (j==MM) then
        iright = 0
      else
        iright = 1
      end if

      do i= 1, 8
        c = IJth(y,i,j)
        clt = IJth(y,i,j+ileft)
        crt = IJth(y,i,j+iright)
        if (i .LT. 8) then
	    diff = D_SP1*(crt - TWO*c + clt)
        else
	    diff = D_SP2*(crt - TWO*c + clt)
	end if
	IJth(dy, i, j) = diff/(y(neq-1)*y(neq-1))
      end do

    end do

    dy(neq) = a0
    dy(neq-1) = mu*y(neq-1)

!do j=1, RNUM
!print *, 'react1_pp: ', pp(j,react1_ita)
!print *, 'aux_pp:    ', pp(j,aux_ita)
!if (a(j) .GT. 0.01) then
!print *, 'j: ', j, ', a(j):', a(j), ', rate: ', RATE(j), ', pop: ', p(NETWORK(REACT1, j))
!end if
!end do


  end subroutine ode_func



  subroutine propensity_update_ode()
    implicit none
    integer :: i, j, ita
    double precision :: propensity

    do i=1, SNUM_ODE
        do j=1, MM
            pp(j,i) = IJth(y, i, j)
        end do
    end do

    h = y(neq-1)

    do i = 1, DPNUM_ODE
        ita = ODEDEPEND(i)
        a0 = a0 - a(ita)
        call cal_propensity_ode(ita, propensity)
        a(ita) = propensity
        a0 = a0 + a(ita)
    end do


  end subroutine propensity_update_ode



  subroutine cal_propensity_ode(ita, propensity)

    implicit none
    integer, intent(in) :: ita
    double precision, intent(out) :: propensity
    integer :: j
    integer :: react_type, react1_ita, react2_ita, aux_ita
    double precision :: rt, temp, Kmdl

    react_type = NETWORK(RTYPE, ita)
    react1_ita = NETWORK(REACT1, ita)
    react2_ita = NETWORK(REACT2, ita)
    aux_ita = NETWORK(AUX, ita)
    rt = RATE(ita)
    temp = 0.0

    select case (react_type)

      case (sticky)
        do j=1, MM
            prePro(j,ita) = pp(j,react1_ita)*pp(j,aux_ita)
            temp = temp + prePro(j,ita)
        end do
        totPro(ita) = temp
        propensity = rt*temp

     case (catalytic)
        do j=1, MM
            prePro(j,ita) = pp(j,react1_ita)*pp(j,aux_ita)
            temp = temp + prePro(j,ita)
        end do
        totPro(ita) = temp
        propensity = rt*temp/h

      case (displacement : combination)
        do j=1, MM
            prePro(j,ita) = pp(j,react1_ita)*pp(j,react2_ita)
            temp = temp + prePro(j,ita)
        end do
        totPro(ita) = temp
        propensity = rt*temp/h

      case (hill)
        Kmdl = 0.5*h*SCALE
        Kmdl = Kmdl*Kmdl*Kmdl*Kmdl
        do j=1,MM
            prePro(j,ita) = pp(j,react1_ita)*(auxt(j)/(auxt(j)+Kmdl))
            temp = temp + prePro(j,ita)
        end do
        totPro(ita) = temp
        propensity = rt*temp

      case (diff)
        propensity= 2*rt/(h*h)*p(react1_ita)

      case default
        print *, "Reaction_type in ita: ", ita, " not found !"
        print *, "Reaction_type in cal-propensity: ", react_type, " not found !"
        stop -1
    
    end select


!if (propensity .GT. 500.0)  then
!    print *, 'react_type: ', react_type, ', react1_ita:', react1_ita, ', react2_ita: ', react2_ita, ', aux_ita: ', aux_ita, ', rate: ', rt, ', h: ', h, ', temp: ', temp, ', react1_pop: ', p(react1_ita), ', aux_pop: ', p(aux_ita)
!do j=1, MM
!print *, 'react1_pp: ', pp(j,react1_ita)
!print *, 'aux_pp:    ', pp(j,aux_ita)
!print *, 'prePro:    ', prePro(j,ita)
!end do
!stop -1
!end if

  end subroutine cal_propensity_ode



end module deterministic
