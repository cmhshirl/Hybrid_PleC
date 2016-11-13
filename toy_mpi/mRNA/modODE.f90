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


    if(mu .NE. 0.0) then
        call propensity_update_ode()
    end if
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

      do i= DivJ_s, SNUM_ODE
        IJth(dy, i, j) = 0.0
      end do


      ! 37 ODE species
      ! 8 diffusive species

      ! DivKp
      c = IJth(y,DivKp,j)
      clt = IJth(y,DivKp,j+ileft)
      crt = IJth(y,DivKp,j+iright)
      diff = D_SP1*(crt - TWO*c + clt)/(h*h)

IJth(dy,DivKp,j) = diff - kdeg_dkp*IJth(y,DivKp,j) &
         + kph1_pc*IJth(y,PleCph1,j) - kpc_ph1*IJth(y,DivKp,j)*IJth(y,PleCphos,j)/h &
         + kpk0_pk1*IJth(y,PleCkin0,j) -  kpk1_pk0*IJth(y,PleCkin1,j)*IJth(y,DivKp,j)/h &
         + kpk1_pk2p*IJth(y,PleCkin1,j) - kpk2p_pk1*IJth(y,PleC2p,j)*IJth(y,DivKp,j)/h &
         + kpt2_pk1h*IJth(y,PleCpt2,j) - kpk1h_pt2*IJth(y,PleCkin01p,j)*IJth(y,DivKp,j)/h &
         + kjkp_dj*IJth(y,DivJKp,j) - kdj_jkp*IJth(y,DivKp,j)*IJth(y,DivJ,j)/h &
         + kpk2_pk3*IJth(y,PleCkin2,j) - kpk3_pk2*IJth(y,PleCkin3,j)*IJth(y,DivKp,j)/h &
         + kpt3_pk1p*IJth(y,PleCkin10p,j) - kpk1p_pt3*IJth(y,PleC1p,j)*IJth(y,DivKp,j)/h &
         - kpk1p_1h*IJth(y,PleC1p,j)*IJth(y,DivKp,j)/h + k1h_pk1p*IJth(y,PleCkin01p,j) &
         - kph2_p12*IJth(y,PleCph2,j)*IJth(y,DivKp,j)/h + kp12_ph2*IJth(y,PleCkin12,j) &
         + kpt4_pk3h*IJth(y,PleCpt4,j) - kpk3h_pt4*IJth(y,PleCkin02p,j)*IJth(y,DivKp,j)/h &
         + kh1_h2*IJth(y,PleCph1,j)*IJth(y,DivK,j)/h - kh2_h1*IJth(y,PleCph2,j)*IJth(y,DivKp,j)/h &
         - kph1_p11*IJth(y,PleCph1,j)*IJth(y,DivKp,j)/h + kp11_ph1*IJth(y,PleCkin11,j) &
         - kbdl_dldk*IJth(y,DivL,j)*IJth(y,DivKp,j)/h + kudl_dldk*IJth(y,DivLKp,j) &
         + kdeg_plc*(IJth(y,PleCph1,j)+2*IJth(y,PleCkin11,j)+2*IJth(y,PleCkin0,j) &
         + IJth(y,PleCkin1,j)+IJth(y,PleCkin12,j)+IJth(y,PleCkin2,j)&
         + 2*IJth(y,PleCpt2,j)+IJth(y,PleCkin01p,j)+IJth(y,PleCkin10p,j)+IJth(y,PleCpt4,j)) &
         - Kdeph_DivKp*IJth(y,DivKp,j) + Kphos_DivK*IJth(y,DivK,j)

      ! DivK
      c = IJth(y,DivK,j)
      clt = IJth(y,DivK,j+ileft)
      crt = IJth(y,DivK,j+iright)
      diff = D_SP1*(crt - TWO*c + clt)/(h*h)

IJth(dy,DivK,j) = diff + ksyn_dk*pp(j,mRNA_DivK)- kdeg_dk*IJth(y,DivK,j) &
         - kpk1_pk2*IJth(y,PleCkin1,j)*IJth(y,DivK,j)/h + kpk2_pk1*IJth(y,PleCkin2,j) &
         - kj_jk*IJth(y,DivJ,j)*IJth(y,DivK,j)/h + kjk_j*IJth(y,DivJK,j) &
         + kpk3_pk2p*IJth(y,PleCkin3,j) - kpk2p_pk3*IJth(y,PleC2p,j)*IJth(y,DivK,j)/h &
         - kpk3_pk4*IJth(y,PleCkin3,j)*IJth(y,DivK,j)/h + kpk4_pk3*IJth(y,PleCkin4,j) &
         - kpc_ph2*IJth(y,PleCphos,j)*IJth(y,DivK,j)/h + kph2_pc*IJth(y,PleCph2,j) &
         - kpk1p_p3h*IJth(y,PleC1p,j)*IJth(y,DivK,j)/h + kp3h_pk1p*IJth(y,PleCkin02p,j) &
         - kph1_p12*IJth(y,PleCph1,j)*IJth(y,DivK,j)/h + kp12_ph1*IJth(y,PleCkin12,j) &
         + kh2_h1*IJth(y,PleCph2,j)*IJth(y,DivKp,j)/h - kh1_h2*IJth(y,PleCph1,j)*IJth(y,DivK,j)/h &
         - kph2_p22*IJth(y,PleCph2,j)*IJth(y,DivK,j)/h + kp22_ph2*IJth(y,PleCkin22,j) &
         + Kdeph_DivKp*IJth(y,DivKp,j) - Kphos_DivK*IJth(y,DivK,j)


! DivJ_f
      c = IJth(y,DivJ_f,j)
      clt = IJth(y,DivJ_f,j+ileft)
      crt = IJth(y,DivJ_f,j+iright)
      diff = D_SP1*(crt - TWO*c + clt)/(h*h)

IJth(dy,DivJ_f,j) = diff + ksyn_divj*pp(j,mRNA_DivJ) - kdeg_divj*IJth(y,DivJ_f,j) &
         - kdj_djp*IJth(y,DivJ_f,j)*IJth(y,DivJ_s,j) + kdjp_dj*IJth(y,DivJ,j)

! CckA free
      c = IJth(y,CckA,j)
      clt = IJth(y,CckA,j+ileft)
      crt = IJth(y,CckA,j+iright)
      diff = D_SP1*(crt - TWO*c + clt)/(h*h)

IJth(dy,CckA,j) = diff + ksyn_ccka*pp(j,mRNA_CckA) - kdeg_ccka*IJth(y,CckA,j) &
          - kcf_cb*IJth(y,CckA,j)*IJth(y,CckA_s,j) + kcb_cf*IJth(y,CckAphos,j)

! CTRA
      c = IJth(y,CtrA,j)
      clt = IJth(y,CtrA,j+ileft)
      crt = IJth(y,CtrA,j+iright)
      diff = D_SP1*(crt - TWO*c + clt)/(h*h)

IJth(dy,CtrA,j) = diff + ksyn_ctr*pp(j,mRNA_CtrA) - kdeg_ctr*IJth(y,CtrA,j) &
          + kctr_phos*IJth(y,CtrAp,j)*IJth(y,CckAphos,j)/h - kctr_kin*IJth(y,CtrA,j)*IJth(y,CckAkin,j)/h

! CTRAP
      c = IJth(y,CtrAp,j)
      clt = IJth(y,CtrAp,j+ileft)
      crt = IJth(y,CtrAp,j+iright)
      diff = D_SP1*(crt - TWO*c + clt)/(h*h)

IJth(dy,CtrAp,j) = diff -kdeg_ctr*IJth(y,CtrAp,j) &
          - kctr_phos*IJth(y,CtrAp,j)*IJth(y,CckAphos,j)/h + kctr_kin*IJth(y,CtrA,j)*IJth(y,CckAkin,j)/h

! DivLfree
      c = IJth(y,DivL_f,j)
      clt = IJth(y,DivL_f,j+ileft)
      crt = IJth(y,DivL_f,j+iright)
      diff = D_SP1*(crt - TWO*c + clt)/(h*h)

IJth(dy,DivL_f,j) = diff + ksyn_dl*pp(j,mRNA_DivL) - kdeg_dl*IJth(y,DivL_f,j) &
          - kdlf_dlb*IJth(y,DivL_s,j)*IJth(y,DivL_f,j) + kdlb_dlf*IJth(y,DivL,j)

! PleC_f
      c = IJth(y,PleC_f,j)
      clt = IJth(y,PleC_f,j+ileft)
      crt = IJth(y,PleC_f,j+iright)
      diff = D_SP2*(crt - TWO*c + clt)/(h*h)

IJth(dy,PleC_f,j) = diff + ksyn_plc2*pp(j,mRNA_PleC) - kdeg_plc2*IJth(y,PleC_f,j) &
         - kplc_plcb*IJth(y,PleC_f,j)*IJth(y,PleC_s,j) + kplcb_plc*IJth(y,PleCphos,j)

! non-diffuse speciese     
! DivL bound
IJth(dy,DivL,j) = - kdeg_dl*IJth(y,DivL,j) &
          + kdlf_dlb*IJth(y,DivL_s,j)*IJth(y,DivL_f,j) - kdlb_dlf*IJth(y,DivL,j) &
          - kbdl_dldk*IJth(y,DivL,j)*IJth(y,DivKp,j)/h + kudl_dldk*IJth(y,DivLKp,j)

! DivL:DivK~P
IJth(dy,DivLKp,j) = - kdeg_dl*IJth(y,DivLKp,j) &
          + kbdl_dldk*IJth(y,DivL,j)*IJth(y,DivKp,j)/h - kudl_dldk*IJth(y,DivLKp,j)

! CCKAphos
IJth(dy,CckAphos,j) = kcf_cb*IJth(y,CckA,j)*IJth(y,CckA_s,j) - kcb_cf*IJth(y,CckAphos,j) &
          - kdeg_ccka*IJth(y,CckAphos,j) &
          - kcp_ck*IJth(y,CckAphos,j)*((IJth(y,DivL,j))**4/((IJth(y,DivL,j))**4+(Kmdl*h*SCALE)**4)) &
          + kck_cp*IJth(y,CckAkin,j) 

! Cckakin
IJth(dy,CckAkin,j) = kcp_ck*IJth(y,CckAphos,j)*((IJth(y,DivL,j))**4/((IJth(y,DivL,j))**4+(Kmdl*h*SCALE)**4)) &
          - kck_cp*IJth(y,CckAkin,j) - kdeg_ccka*IJth(y,CckAkin,j) 

!DIVJ
IJth(dy,DivJ,j) = kdj_djp*IJth(y,DivJ_f,j)*IJth(y,DivJ_s,j) - kdjp_dj*IJth(y,DivJ,j) &
         - kj_jk*IJth(y,DivJ,j)*IJth(y,DivK,j)/h + kjk_j*IJth(y,DivJK,j) &
         + kjkp_dj*IJth(y,DivJKp,j) - kdj_jkp*IJth(y,DivKp,j)*IJth(y,DivJ,j)/h &
         - kdegpp2*IJth(y,DivJ,j)

!DJk
IJth(dy,DivJK,j) = kj_jk*IJth(y,DivJ,j)*IJth(y,DivK,j)/h - kjk_j*IJth(y,DivJK,j) &
         - kjk_jkp*IJth(y,DivJK,j) + kjkp_jk*IJth(y,DivJKp,j) &
         - kdegpp2*IJth(y,DivJK,j)

! DIvJ:DivK-P
IJth(dy,DivJKp,j) = kjk_jkp*IJth(y,DivJK,j) - kjkp_jk*IJth(y,DivJKp,j) &
          - kjkp_dj*IJth(y,DivJKp,j) + kdj_jkp*IJth(y,DivKp,j)*IJth(y,DivJ,j)/h &
          - kdegpp2*IJth(y,DivJKp,j)

! PLEC
IJth(dy,PleCphos,j) = - kdeg_plc*IJth(y,PleCphos,j) &
         - kpc_ph1*IJth(y,PleCphos,j)*IJth(y,DivKp,j)/h + kph1_pc*IJth(y,PleCph1,j) &
         + kpk2p_pc*IJth(y,PleC2p,j) + kpk1p_pc*IJth(y,PleC1p,j) &
         - kpc_ph2*IJth(y,PleCphos,j)*IJth(y,DivK,j)/h + kph2_pc*IJth(y,PleCph2,j) &
         + kplc_plcb*IJth(y,PleC_f,j)*IJth(y,PleC_s,j) - kplcb_plc*IJth(y,PleCphos,j)    

! PLECH1
IJth(dy,PleCph1,j) = kpc_ph1*IJth(y,PleCphos,j)*IJth(y,DivKp,j)/h - kph1_pc*IJth(y,PleCph1,j) &
         - kph1_p11*IJth(y,PleCph1,j)*IJth(y,DivKp,j)/h + kp11_ph1*IJth(y,PleCkin11,j) &
         - kph1_p12*IJth(y,PleCph1,j)*IJth(y,DivK,j)/h + kp12_ph1*IJth(y,PleCkin12,j) &
         + kh2_h1*IJth(y,PleCph2,j)*IJth(y,DivKp,j)/h - kh1_h2*IJth(y,PleCph1,j)*IJth(y,DivK,j)/h &
         - kph1_ph2*IJth(y,PleCph1,j) + kph2_ph1*IJth(y,PleCph2,j) &
         - kdeg_plc*IJth(y,PleCph1,j)

! PLECH2
IJth(dy,PleCph2,j) = kpc_ph2*IJth(y,PleCphos,j)*IJth(y,DivK,j)/h - kph2_pc*IJth(y,PleCph2,j) &
         - kph2_p22*IJth(y,PleCph2,j)*IJth(y,DivK,j)/h + kp22_ph2*IJth(y,PleCkin22,j)&
         - kph2_p12*IJth(y,PleCph2,j)*IJth(y,DivKp,j)/h + kp12_ph1*IJth(y,PleCkin12,j) &
         + kh1_h2*IJth(y,PleCph1,j)*IJth(y,DivK,j)/h - kh2_h1*IJth(y,PleCph2,j)*IJth(y,DivKp,j)/h &
         + kph1_ph2*IJth(y,PleCph1,j) - kph2_ph1*IJth(y,PleCph2,j) &
         - kdeg_plc*IJth(y,PleCph2,j)

! Pk11
IJth(dy,PleCkin11,j) = kph1_p11*IJth(y,PleCph1,j)*IJth(y,DivKp,j)/h- kp11_ph1*IJth(y,PleCkin11,j) &
         + kpt4_p11*IJth(y,PleCpt4,j) - kp11_pt4*IJth(y,PleCkin11,j) &
         - kp11_pk0*IJth(y,PleCkin11,j) + kpk0_p11*IJth(y,PleCkin0,j) &
         - kdeg_plc*IJth(y,PleCkin11,j)

!Pk12
IJth(dy,PleCkin12,j) = kph1_p12*IJth(y,PleCph1,j)*IJth(y,DivK,j)/h- kp12_ph1*IJth(y,PleCkin12,j) &
         + kph2_p12*IJth(y,PleCph2,j)*IJth(y,DivKp,j)/h- kp12_ph2*IJth(y,PleCkin12,j) &
         - kp12_pk2*IJth(y,PleCkin12,j) + kpk2_p12*IJth(y,PleCkin2,j) &
         - kdeg_plc*IJth(y,PleCkin12,j)

! Pk22
IJth(dy,PleCkin22,j) = kph2_p22*IJth(y,PleCph2,j)*IJth(y,DivK,j)/h - kp22_ph2*IJth(y,PleCkin22,j) &
         - kp22_pk4*IJth(y,PleCkin22,j) + kpk4_p22*IJth(y,PleCkin4,j) &
         - kdeg_plc*IJth(y,PleCkin22,j)

! Pk0
IJth(dy,PleCkin0,j) = kp11_pk0*IJth(y,PleCkin11,j) - kpk0_p11*IJth(y,PleCkin0,j) &
         - kpk0_pk1*IJth(y,PleCkin0,j) + kpk1_pk0*IJth(y,PleCkin1,j)*IJth(y,DivKp,j)/h &
         - kdeg_plc*IJth(y,PleCkin0,j)

!Pk1
IJth(dy,PleCkin1,j) = kpk0_pk1*IJth(y,PleCkin0,j) -  kpk1_pk0*IJth(y,PleCkin1,j)*IJth(y,DivKp,j)/h  &
         - kpk1_pk2*IJth(y,PleCkin1,j)*IJth(y,DivK,j)/h + kpk2_pk1*IJth(y,PleCkin2,j) &
         - kpk1_pk2p*IJth(y,PleCkin1,j) + kpk2p_pk1*IJth(y,PleC2p,j)*IJth(y,DivKp,j)/h &
         - kpk1_pk1h*IJth(y,PleCkin1,j) + kpk1h_pk1*IJth(y,PleCkin01p,j) &
         - kdeg_plc*IJth(y,PleCkin1,j)

! Pk2
IJth(dy,PleCkin2,j) = kpk1_pk2*IJth(y,PleCkin1,j)*IJth(y,DivK,j)/h - kpk2_pk1*IJth(y,PleCkin2,j) &
         + kp12_pk2*IJth(y,PleCkin12,j) - kpk2_p12*IJth(y,PleCkin2,j) &
         - kpk2_pt2*IJth(y,PleCkin2,j) + kpt2_pk2*IJth(y,PleCpt2,j) &
         - kpk2_pk3*IJth(y,PleCkin2,j) + kpk3_pk2*IJth(y,PleCkin3,j)*IJth(y,DivKp,j)/h &
         - kdeg_plc*IJth(y,PleCkin2,j) 

! Pk3
IJth(dy,PleCkin3,j) = kpk2_pk3*IJth(y,PleCkin2,j) - kpk3_pk2*IJth(y,PleCkin3,j)*IJth(y,DivKp,j)/h &
         - kpk3_pt3*IJth(y,PleCkin3,j) + kpt3_pk3*IJth(y,PleCkin10p,j) &
         - kpk3_pk2p*IJth(y,PleCkin3,j) + kpk2p_pk3*IJth(y,PleC2p,j)*IJth(y,DivK,j)/h &
         - kpk3_pk4*IJth(y,PleCkin3,j)*IJth(y,DivK,j)/h + kpk4_pk3*IJth(y,PleCkin4,j) &
         - kpk3_pk3h*IJth(y,PleCkin3,j) + kpk3h_pk3*IJth(y,PleCkin02p,j) &
         - kdeg_plc*IJth(y,PleCkin3,j)

! Pk4
IJth(dy,PleCkin4,j) = kpk3_pk4*IJth(y,PleCkin3,j)*IJth(y,DivK,j)/h - kpk4_pk3*IJth(y,PleCkin4,j)&
         - kpk4_pt4*IJth(y,PleCkin4,j) + kpt4_pk4*IJth(y,PleCpt4,j) &
         + kp22_pk4*IJth(y,PleCkin22,j) - kpk4_p22*IJth(y,PleCkin4,j) &
         - kdeg_plc*IJth(y,PleCkin4,j)

! Pk1P
IJth(dy,PleC1p,j) = kpt3_pk1p*IJth(y,PleCkin10p,j) - kpk1p_pt3*IJth(y,PleC1p,j)*IJth(y,DivKp,j)/h &
         - kpk1p_1h*IJth(y,PleC1p,j)*IJth(y,DivKp,j)/h + k1h_pk1p*IJth(y,PleCkin01p,j) &
         - kpk1p_pc*IJth(y,PleC1p,j) &
         - kpk1p_p3h*IJth(y,PleC1p,j)*IJth(y,DivK,j)/h + kp3h_pk1p*IJth(y,PleCkin02p,j) &
         - kdeg_plc*IJth(y,PleC1p,j)

! Pk2P
IJth(dy,PleC2p,j) = kpk1_pk2p*IJth(y,PleCkin1,j) - kpk2p_pk1*IJth(y,PleC2p,j)*IJth(y,DivKp,j)/h &
         - kpk2p_pc*IJth(y,PleC2p,j)  + kpk3_pk2p*IJth(y,PleCkin3,j) &
         - kpk2p_pk3*IJth(y,PleC2p,j)*IJth(y,DivK,j)/h &
         - kdeg_plc*IJth(y,PleC2p,j)

! Pk1H
IJth(dy,PleCkin01p,j) =  kpt2_pk1h*IJth(y,PleCpt2,j) - kpk1h_pt2*IJth(y,PleCkin01p,j)*IJth(y,DivKp,j)/h &
         + kpk1_pk1h*IJth(y,PleCkin1,j) - kpk1h_pk1*IJth(y,PleCkin01p,j) &
         + kpk1p_1h*IJth(y,PleC1p,j)*IJth(y,DivKp,j)/h - k1h_pk1p*IJth(y,PleCkin01p,j) &
         - kdeg_plc*IJth(y,PleCkin01p,j)

! PT3
IJth(dy,PleCkin10p,j) = kpk3_pt3*IJth(y,PleCkin3,j) - kpt3_pk3*IJth(y,PleCkin10p,j) &
         + kpk1p_pt3*IJth(y,PleC1p,j)*IJth(y,DivKp,j)/h - kpt3_pk1p*IJth(y,PleCkin10p,j) &
         - kdeg_plc*IJth(y,PleCkin10p,j)

! Pk3H<=>PleCkin02p
IJth(dy,PleCkin02p,j) = kpt4_pk3h*IJth(y,PleCpt4,j) - kpk3h_pt4*IJth(y,PleCkin02p,j)*IJth(y,DivKp,j)/h &
+ kpk3_pk3h*IJth(y,PleCkin3,j) - kpk3h_pk3*IJth(y,PleCkin02p,j) &
+ kpk1p_p3h*IJth(y,PleC1p,j)*IJth(y,DivK,j)/h - kp3h_pk1p*IJth(y,PleCkin02p,j) &
- kdeg_plc*IJth(y,PleCkin02p,j)

! PT2 <=> PleCpt2
IJth(dy,PleCpt2,j) = kpk2_pt2*IJth(y,PleCkin2,j) - kpt2_pk2*IJth(y,PleCpt2,j) &
         - kpt2_pk1h*IJth(y,PleCpt2,j) + kpk1h_pt2*IJth(y,PleCkin01p,j)*IJth(y,DivKp,j)/h &
         - kdeg_plc*IJth(y,PleCpt2,j)

! PT4 <=> PleCpt4
IJth(dy,PleCpt4,j) = kpk4_pt4*IJth(y,PleCkin4,j) - kpt4_pk4*IJth(y,PleCpt4,j) &
         - kpt4_p11*IJth(y,PleCpt4,j) + kp11_pt4*IJth(y,PleCkin11,j)&
         - kpt4_pk3h*IJth(y,PleCpt4,j) + kpk3h_pt4*IJth(y,PleCkin02p,j)*IJth(y,DivKp,j)/h &
         - kdeg_plc*IJth(y,PleCpt4,j)


    end do

    dy(neq) = a0
    dy(neq-1) = mu*y(neq-1)


  end subroutine ode_func

subroutine propensity_update_ode()
    implicit none
    integer :: i, j, ita
    double precision :: propensity

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
