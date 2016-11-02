# define IJth(y,i,j) y((i-1)*5+j)

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
    integer :: i, j, ileft, iright, TWO=2
    double precision :: c, clt, crt, avg, Kmdl, diff
    double precision, intent(in) ::t, y(neq)
    double precision, intent(out)::dy(neq)


    Kmdl = 0.25*h*SCALE
    Kmdl = Kmdl*Kmdl*Kmdl*Kmdl
    call propensity_update_ode()

    !print *, 't     ', t
    !print *, 'a0:', a0
    !print *, 'y(neq)', y(neq)

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
 
      do i= 10, 14
        c = IJth(y,i,j)
        clt = IJth(y,i,j+ileft)
        crt = IJth(y,i,j+iright)
        diff = D_mRNA*(crt - TWO*c + clt)
        IJth(dy, i, j) = diff/(y(neq-1)*y(neq-1)) + K_syn_mRNA*pp(j,gene_ori) + K_syn_mRNA*pp(j,gene_dup) - K_deg_mRNA*IJth(y, i, j)
      end do

      i = mRNA_DivJ
      c = IJth(y,i,j)
      clt = IJth(y,i,j+ileft)
      crt = IJth(y,i,j+iright)
      diff = D_mRNA*(crt - TWO*c + clt)
      IJth(dy, i, j) = diff/(y(neq-1)*y(neq-1)) + K_syn_mRNADivJ*pp(j,gene_ori) + K_syn_mRNADivJ*pp(j,gene_dup) - K_deg_mRNA*IJth(y, i, j)

      !DivKp 
      i = DivKp
      c = IJth(y,i,j)
      clt = IJth(y,i,j+ileft)
      crt = IJth(y,i,j+iright)  
      diff = D_SP*(crt - TWO*c + clt)    
      IJth(dy, i, j) = diff/(y(neq-1)*y(neq-1)) - K_bndKp_DivL*IJth(y, DivL, j)*IJth(y, i, j) + K_ubdKp_DivLKp*IJth(y, DivLKp, j) 

      !DivK 
      i = DivK
      c = IJth(y,i,j)
      clt = IJth(y,i,j+ileft)
      crt = IJth(y,i,j+iright)  
      diff = D_SP*(crt - TWO*c + clt)    
      IJth(dy, i, j) = diff/(y(neq-1)*y(neq-1)) + K_syn_DivK*IJth(y, mRNA_DivK, j)

      !DivJ_f 
      i = DivJ_f
      c = IJth(y,i,j)
      clt = IJth(y,i,j+ileft)
      crt = IJth(y,i,j+iright)  
      diff = D_SP*(crt - TWO*c + clt)    
      IJth(dy, i, j) = diff/(y(neq-1)*y(neq-1)) + K_syn_DivJf*IJth(y, mRNA_DivJ, j)

      !CckA 
      i = CckA
      c = IJth(y,i,j)
      clt = IJth(y,i,j+ileft)
      crt = IJth(y,i,j+iright)  
      diff = D_SP*(crt - TWO*c + clt)
      IJth(dy, i, j) = diff/(y(neq-1)*y(neq-1)) + K_syn_CckA*IJth(y, mRNA_CckA, j) - K_sticky*IJth(y, CckA, j)*pp(j,CckA_st)

      !CtrA 
      i = CtrA
      c = IJth(y,i,j)
      clt = IJth(y,i,j+ileft)
      crt = IJth(y,i,j+iright)  
      diff = D_SP*(crt - TWO*c + clt)    
      IJth(dy, i, j) = diff/(y(neq-1)*y(neq-1)) - K_phos_CtrA*IJth(y, i, j)*IJth(y, CckA_kin, j) + K_deph_CtrAp*IJth(y, CtrAp, j)*IJth(y, CckA_phos, j) + K_syn_CtrA*IJth(y, mRNA_CtrA, j) - K_deg_CtrA*IJth(y, CtrA, j)

      !CtrAp 
      i = CtrAp
      c = IJth(y,i,j)
      clt = IJth(y,i,j+ileft)
      crt = IJth(y,i,j+iright)  
      diff = D_SP*(crt - TWO*c + clt)    
      IJth(dy, i, j) = diff/(y(neq-1)*y(neq-1)) + K_phos_CtrA*IJth(y, CtrA, j)*IJth(y, CckA_kin, j) - K_deph_CtrAp*IJth(y, i, j)*IJth(y, CckA_phos, j) - K_deg_CtrAp*IJth(y, CtrAp, j)
            
      !DivL_f 
      i = DivL_f
      c = IJth(y,i,j)
      clt = IJth(y,i,j+ileft)
      crt = IJth(y,i,j+iright)  
      diff = D_SP*(crt - TWO*c + clt)    
      IJth(dy, i, j) = diff/(y(neq-1)*y(neq-1)) + K_syn_DivLf*IJth(y, mRNA_DivL, j) - K_deg_DivLf*IJth(y, i, j) - K_bnd_DivL*IJth(y, i, j)*pp(j,DivL_st) + K_ubd_DivL*IJth(y, DivL, j)

      !PleC_f 
      i = PleC_f
      c = IJth(y,i,j)
      clt = IJth(y,i,j+ileft)
      crt = IJth(y,i,j+iright)  
      diff = D_SP*(crt - TWO*c + clt)    
      IJth(dy, i, j) = diff/(y(neq-1)*y(neq-1)) + K_syn_PleCf*IJth(y, mRNA_PleC, j) - K_deg_PleCf*IJth(y, PleC_f, j)

      !/******************* no diffusion *********************/      
      !DivL 
      i = DivL
      IJth(dy, i, j) = K_bnd_DivL*IJth(y, DivL_f, j)*pp(j,DivL_st) - K_ubd_DivL*IJth(y, i, j) - K_deg_DivL*IJth(y, i, j)

      !DivLKp 
      i = DivLKp  
      IJth(dy, i, j) = K_bndKp_DivL*IJth(y, DivL, j)*IJth(y, DivKp, j) - K_ubdKp_DivLKp*IJth(y, i, j)*IJth(y, DivKp, j) - K_deg_DivLKp*IJth(y, i, j)

      !CckA_phos 
      i = CckA_phos
      avg = IJth(y, DivL, j)
      avg = avg*avg*avg*avg
      IJth(dy, i, j) = - K_ph2k_CckA*IJth(y, i, j)*avg/(Kmdl+avg) + K_sticky*IJth(y, CckA, j)*pp(j,CckA_st) + K_k2ph_CckA*IJth(y, CckA_kin, j)

      !CckA_kin 
      i = CckA_kin
      IJth(dy, i, j) = K_ph2k_CckA*IJth(y, CckA_phos, j)*avg/(Kmdl+avg) - K_k2ph_CckA*IJth(y, CckA_kin, j)

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
        p(i) = 0.0
        do j=1, MM
            pp(j,i) = IJth(y, i, j)
            p(i) = p(i) + IJth(y, i, j)
        end do
    end do

    h = y(neq-1)

    do i = 1, DPNUM_ODE
        ita = ODEDEPEND(i)
        a0 = a0 - a(ita)
        call cal_propensity(ita, propensity)
        a(ita) = propensity
        a0 = a0 + a(ita)
    end do


  end subroutine propensity_update_ode



  subroutine cal_propensity(ita, propensity)

    implicit none
    integer, intent(in) :: ita
    double precision, intent(out) :: propensity
    integer :: j
    integer :: react_type, react1_ita, react2_ita, aux_ita
    double precision :: rt, temp

    react_type = NETWORK(RTYPE, ita)
    react1_ita = NETWORK(REACT1, ita)
    react2_ita = NETWORK(REACT2, ita)
    aux_ita = NETWORK(AUX, ita)
    rt = RATE(ita)
    temp = 0.0

    select case (react_type)
    
      case (syn)
        propensity = rt * p(aux_ita)

      case (deg : decomposition)
        propensity = rt * p(react1_ita)
      
      case (sticky)
        do j=1, MM
            prePro(j,ita) = pp(j,react1_ita)*pp(j,aux_ita)
            temp = temp + prePro(j,ita)
        end do
        totPro(ita) = temp
        propensity = rt*temp
      
      case (displacement : combination)
        do j=1, MM
            prePro(j,ita) = pp(j,react1_ita)*pp(j,react2_ita)
            temp = temp + prePro(j,ita)
        end do
        totPro(ita) = temp
        propensity = rt*temp/h
      
      case default
        print *, "Reaction_type in ita: ", ita, " not found !"
        print *, "Reaction_type in cal-propensity: ", react_type, " not found !"
        stop -1
    
    end select


if (propensity .GT. 500.0)  then
    print *, 'react_type: ', react_type, ', react1_ita:', react1_ita, ', react2_ita: ', react2_ita, ', aux_ita: ', aux_ita, ', rate: ', rt, ', h: ', h, ', temp: ', temp, ', react1_pop: ', p(react1_ita), ', aux_pop: ', p(aux_ita)
do j=1, MM
!print *, 'react1_pp: ', pp(j,react1_ita)
!print *, 'aux_pp:    ', pp(j,aux_ita)
print *, 'prePro:    ', prePro(j,ita)
end do
stop -1
end if

  end subroutine cal_propensity



end module deterministic
