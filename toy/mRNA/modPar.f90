module parameters
  implicit none

  !************************
  !*                      *
  !* Parameters for system*
  !*                      *
  !************************
  integer, parameter:: MM = 50
  integer, parameter:: SNUM_SSA = 8
  integer, parameter:: SNUM_ODE = 37
  integer, parameter:: RNUM = 24
  integer, parameter:: DPNUM_ODE = 6
  integer, parameter:: DPNUM_SSA = 4

  integer, parameter:: HMEAN = 10
  integer, parameter:: syn = 1
  integer, parameter:: deg = 2
  integer, parameter:: first_order= 3
  integer, parameter:: decomposition = 4
  integer, parameter:: sticky = 5
  integer, parameter:: catalytic = 6
  integer, parameter:: displacement = 7
  integer, parameter:: combination = 8
  integer, parameter:: hill = 9
  integer, parameter:: diff = 10

  !************************
  !*                      *
  !* Parameters for SPECIE*
  !*                      *
  !************************
  integer, parameter:: RSPACE = 6
  integer, parameter:: RTYPE  = 1
  integer, parameter:: REACT1 = 2
  integer, parameter:: REACT2 = 3
  integer, parameter:: PROD1  = 4
  integer, parameter:: PROD2  = 5
  integer, parameter:: AUX    = 6
  
  integer, parameter :: DivKp   = 1
  integer, parameter :: DivK    = 2
  integer, parameter :: DivJ_f  = 3
  integer, parameter :: CckA    = 4
  integer, parameter :: CtrA    = 5
  integer, parameter :: CtrAp   = 6
  integer, parameter :: DivL_f  = 7
  integer, parameter :: PleC_f	= 8

  integer, parameter :: DivL = 9
  integer, parameter :: DivLKp = 10
  integer, parameter :: CckAphos = 11
  integer, parameter :: CckAkin = 12
  integer, parameter :: DivJ = 13
  integer, parameter :: DivJK = 14
  integer, parameter :: DivJKp = 15
  integer, parameter :: PleCphos = 16
  integer, parameter :: PleCph1 = 17
  integer, parameter :: PleCph2 = 18
  integer, parameter :: PleCkin11 = 19
  integer, parameter :: PleCkin12 = 20
  integer, parameter :: PleCkin22 = 21
  integer, parameter :: PleCkin0 = 22
  integer, parameter :: PleCkin1 = 23
  integer, parameter :: PleCkin2 = 24
  integer, parameter :: PleCkin3 = 25
  integer, parameter :: PleCkin4 = 26
  integer, parameter :: PleC1p = 27
  integer, parameter :: PleC2p = 28
  integer, parameter :: PleCkin01p = 29
  integer, parameter :: PleCkin10p = 30
  integer, parameter :: PleCkin02p = 31
  integer, parameter :: PleCpt2 = 32
  integer, parameter :: PleCpt4 = 33
  integer, parameter :: DivJ_s = 34
  integer, parameter :: DivL_s = 35
  integer, parameter :: CckA_s = 36
  integer, parameter :: PleC_s = 37

  integer, parameter :: mRNA_DivJ	= 1
  integer, parameter :: mRNA_DivK	= 2
  integer, parameter :: mRNA_PleC	= 3
  integer, parameter :: mRNA_DivL	= 4
  integer, parameter :: mRNA_CckA	= 5
  integer, parameter :: mRNA_CtrA	= 6
  
  integer, parameter :: gene_ori	= 7
  integer, parameter :: gene_dup	= 8

  integer, parameter :: syn_OrimrDivJ = 7
  integer, parameter :: syn_DupmrDivJ = 8


  !************************
  !*                      *
  !*Para for ODE RATE     *
  !*                      *
  !************************
  double precision, parameter :: SCALE = 1.0e+3
  double precision, parameter :: D_SP1           = 1.0e+0
  double precision, parameter :: D_SP2           = 1.0e+0  
  !add in matlab version reactions
  
  ! Synthesis rate constants [units --> 1/min]
  double precision, parameter :: ksyn_dk = 25.0      !null -> DivK
  double precision, parameter :: kdeg_dk = 0.005     !DivK -> null
  double precision, parameter :: kdeg_dkp = 0.005    !DivKp -> null

  double precision, parameter :: ksyn_plc2 = 50.0     ! 0 for delta_pleC
  double precision, parameter :: kdeg_plc = 0.05
  double precision, parameter :: kdeg_plc2 = 0.05

  double precision, parameter :: ksyn_divj = 25.0
  double precision, parameter :: kdeg_divj = 0.05
  double precision, parameter :: kdegpp2 = 0.05
  
  double precision, parameter :: ksyn_ccka = 12.5
  double precision, parameter :: kdeg_ccka= 0.05
  
  double precision, parameter :: ksyn_ctr = 25.0
  double precision, parameter :: kdeg_ctr = 0.05
  
  double precision, parameter :: ksyn_dl = 12.5
  double precision, parameter :: kdeg_dl = 0.05

  ! Binding constants
  double precision, parameter :: kdj_djp = 1.0          ! DivJ binding
  double precision, parameter :: kdjp_dj = 0
  
  double precision, parameter :: kplc_plcb = 1.0        ! pleC binding
  double precision, parameter :: kplcb_plc = 0.5
  
  double precision, parameter :: kdlf_dlb = 1.0         ! DivL binding
  double precision, parameter :: kdlb_dlf = 0.1

  double precision, parameter :: kbdl_dldk = 2.5/SCALE      ! DivL + DivKp -> DivLKp
  double precision, parameter :: kudl_dldk = 0.5      !              <-

  double precision, parameter :: kcf_cb = 1.0           ! CckA binding
  double precision, parameter :: kcb_cf = 0.1

  double precision, parameter :: kcp_ck=10.0            ! CckA phos-to-kin
  double precision, parameter :: kck_cp=1.0
  double precision, parameter :: Kmdl=0.5

  ! DivJ phosphorylation paramters
  double precision, parameter :: kj_jk = 5.0/SCALE            ! DIvJ + DivK -> DivJK
  double precision, parameter :: kjk_j = 0.0016	!             <-
  double precision, parameter :: kjk_jkp = 5.0
  double precision, parameter :: kjkp_jk = 0.16
  double precision, parameter :: kjkp_dj = 1.0
  double precision, parameter :: kdj_jkp = 5.0/SCALE      !DivJ + DivKp -> DivJKp

  ! PleC phosphorylation paramters
  double precision, parameter :: kpc_ph1 = 5.0/SCALE	    !       <-
  double precision, parameter :: kph1_pc = 5.0	    !PleCph1 -> PleCphos + DivKp

  double precision, parameter :: kph1_p11 = 5.0/SCALE	   ! PleCph1 + DivKp -> PleCkin11
  double precision, parameter :: kp11_ph1 = 2.5  !                 <-

  double precision, parameter :: kp11_pk0 = 2.5    !          <-
  double precision, parameter :: kpk0_p11 = 5.0      ! PleCkin11 -> PleCkin0

  double precision, parameter :: kpk0_pk1 = 0.16   ! PleCkin0 -> PleCkin1 + DivKp
  double precision, parameter :: kpk1_pk0 = 5.0/SCALE      !          <-

  double precision, parameter :: kpk1_pk2 = 5.0/SCALE      !PleCkin1 + DivK -> PleCkin2
  double precision, parameter :: kpk2_pk1 = 0.0016 !                <-

  double precision, parameter :: kpk1_pk1h = 5.0     ! PleCkin1 -> PleCkin01p
  double precision, parameter :: kpk1h_pk1 = 5.0     !          <-

  double precision, parameter :: kpk1_pk2p = 0.16  !PleCkin1 -> PleC2p + DivKp
  double precision, parameter :: kpk2p_pk1 = 5.0/SCALE     !         <-

  double precision, parameter :: kpk2_pt2 = 5.0	     ! PleCkin2 -> PleCpt2
  double precision, parameter :: kpt2_pk2 = 0.16   !          <-

  double precision, parameter :: kpt2_pk1h = 0.16  ! PleCpt2 -> PleCkin01p + DivKp
  double precision, parameter :: kpk1h_pt2 = 5.0/SCALE     !         <-

  double precision, parameter :: kpk2p_pc = 5.0	     ! PleC2p -> PleCphos
  double precision, parameter :: kpk1p_pc = 5.0     ! PleC1p -> PleCphos

  double precision, parameter :: kpk2_pk3 = 0.16   !PleCkin2 -> PleCkin3 + DivKp
  double precision, parameter :: kpk3_pk2 = 5.0/SCALE      !         <-

  double precision, parameter :: kpk3_pt3 = 5.0	     !PleCkin3 -> pleCpt3
  double precision, parameter :: kpt3_pk3 = 0.16   !         <-

  double precision, parameter :: kpt3_pk1p = 0.16 !PleCkin10p -> PleC1p + DivKp
  double precision, parameter :: kpk1p_pt3 = 5.0/SCALE    !        <-

  double precision, parameter :: kpk1p_1h = 5.0/SCALE     !PleC1p + DivKp ->PleCkin01p
  double precision, parameter :: k1h_pk1p = 0.16  !              <-

  double precision, parameter :: kpk3_pk2p = 0.0016! PleCkin3  -> PleC2p + DivK
  double precision, parameter :: kpk2p_pk3 = 5.0/SCALE     !           <-

  double precision, parameter :: kpk3_pk4 = 5.0/SCALE	     ! PleCkin3 + DivK -> PleCkin4
  double precision, parameter :: kpk4_pk3 = 0.0016 !                 <-

  double precision, parameter :: kpk4_pt4 = 5.0	     ! PleCkin4 -> PleCpt4
  double precision, parameter :: kpt4_pk4 = 0.16   !          <-

  double precision, parameter :: kpc_ph2 = 0.05/SCALE    ! PleCphos + DivK -> PleCph2
  double precision, parameter :: kph2_pc = 5.0       !                 <-

  double precision, parameter :: kph2_p22 = 0.016/SCALE	!PleCph2 + DivK -> PleCkin22
  double precision, parameter :: kp22_ph2 = 1.6e-08   !              <-

  double precision, parameter :: kp22_pk4 = 5.0	     ! PleCkin22 -> PleCkin4
  double precision, parameter :: kpk4_p22 = 5.0      !           <-

  double precision, parameter :: kpt4_pk3h = 0.16  ! PleCpt4 -> PleCkin02p + DivKp
  double precision, parameter :: kpk3h_pt4 = 5.0/SCALE     !         <-

  double precision, parameter :: kpk3_pk3h = 5.0     ! PleCkin3 -> PleCkin02p
  double precision, parameter :: kpk3h_pk3 = 5.0     !          <-

  double precision, parameter :: kpk1p_p3h = 5.0/SCALE      !PleC1p + DivK -> PleCkin02p
  double precision, parameter :: kp3h_pk1p = 0.0016 !              <-

  double precision, parameter :: kph1_p12 = 1.6e-02/SCALE ! PleCph1 + DivK -> PleCkin12
  double precision, parameter :: kp12_ph1 = 1.6e-04 !                <-

  double precision, parameter :: kp12_pk2 = 5.0	      ! PleCkin12 -> PleCkin2
  double precision, parameter :: kpk2_p12 = 5.0      !           <-

  double precision, parameter :: kph2_p12 = 1.6/SCALE	!PleCph2 + DivKp -> Plekin12
  double precision, parameter :: kp12_ph2 = 1.6e-04   !                <-

  double precision, parameter :: kh1_h2 = 1.6e-02/SCALE !               <-
  double precision, parameter :: kh2_h1 = 1.6/SCALE      !PleCph2 + DivKp -> PleCph1 + DivK

  double precision, parameter :: kp11_pt4 = 0.0755	!PleCkin11 -> PleCpt4
  double precision, parameter :: kpt4_p11 = 5.0         !          <-

  double precision, parameter :: kph1_ph2 = 10.0        !PleCph1 -> PleCph2
  double precision, parameter :: kph2_ph1 = 5e-03     !        <-

  ! CtrA phosphorylation
  double precision, parameter :: kctr_kin = 600.0/SCALE
  double precision, parameter :: kctr_phos = 600.0/SCALE

  !two new add in reactions
  double precision :: Kphos_DivK = 0.0
  double precision :: Kdeph_DivKp = 0.0

  !**********************
  !*                    *
  !* changing variables *
  !*                    *
  !**********************
  double precision :: h
  double precision :: a0
  double precision, dimension(RNUM) :: a
  double precision, dimension(SNUM_SSA) :: p
  double precision, dimension(MM, SNUM_SSA) :: pp
  double precision:: mu
  integer(KIND=8), dimension(10) :: fire

  !**********************
  !*                    *
  !* Param for SSA rate *
  !*                    *
  !**********************
  double precision, dimension(RNUM) :: RATE
  
  integer, dimension(DPNUM_ODE), parameter :: ODEDEPEND = (/ 1,	2,	3,	4,	5,	6 /)
    
  integer, dimension(RSPACE, RNUM), parameter :: NETWORK = reshape( (/ &
10, 1,  0,  0,  0,  0, &
10, 2,  0,  0,  0,  0, &
10, 3,  0,  0,  0,  0, &
10, 4,  0,  0,  0,  0, &
10, 5,  0,  0,  0,  0, &
10, 6,  0,  0,  0,  0, &
1,  0,  0,  1,  0,  7, &
1,  0,  0,  1,  0,  8, &
2,  1,  0,  0,  0,  0, &
1,  0,  0,  2,  0,  7, &
1,  0,  0,  2,  0,  8, &
2,  2,  0,  0,  0,  0, &
1,  0,  0,  3,  0,  7, &
1,  0,  0,  3,  0,  8, &
2,  3,  0,  0,  0,  0, &
1,  0,  0,  4,  0,  7, &
1,  0,  0,  4,  0,  8, &
2,  4,  0,  0,  0,  0, &
1,  0,  0,  5,  0,  7, &
1,  0,  0,  5,  0,  8, &
2,  5,  0,  0,  0,  0, &
1,  0,  0,  6,  0,  7, &
1,  0,  0,  6,  0,  8, &
2,  6,  0,  0,  0,  0 /), (/ RSPACE, RNUM /) )

  integer, dimension(DPNUM_SSA,RNUM) :: SSADEPEND


contains

  subroutine init_para()

    RATE = (/ 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.625,  0.625,  0.25, 0.625,  0.625,  0.25, &
      0.625,  0.625,  0.25, 0.625,  0.625,  0.25, 0.625,  0.625,  0.25, 0.625,  0.625,  0.25 /)

SSADEPEND(:,1) = (/ 1, 9, 0, 0 /)
SSADEPEND(:,2) = (/ 1, 12, 0, 0 /)
SSADEPEND(:,3) = (/ 1, 15, 0, 0 /)
SSADEPEND(:,4) = (/ 1, 18, 0, 0 /)
SSADEPEND(:,5) = (/ 1, 21, 0, 0 /)
SSADEPEND(:,6) = (/ 1, 24, 0, 0 /)
SSADEPEND(:,7) = (/ 3, 1, 7, 9 /)
SSADEPEND(:,8) = (/ 3, 1, 8, 9 /)
SSADEPEND(:,9) = (/ 2, 1, 9, 0 /)
SSADEPEND(:,10) = (/ 3, 2, 10, 12 /)
SSADEPEND(:,11) = (/ 3, 2, 11, 12 /)
SSADEPEND(:,12) = (/ 2, 2, 12, 0 /)
SSADEPEND(:,13) = (/ 3, 3, 13, 15 /)
SSADEPEND(:,14) = (/ 3, 3, 14, 15 /)
SSADEPEND(:,15) = (/ 2, 3, 15, 0 /)
SSADEPEND(:,16) = (/ 3, 4, 16, 18 /)
SSADEPEND(:,17) = (/ 3, 4, 17, 18 /)
SSADEPEND(:,18) = (/ 2, 4, 18, 0 /)
SSADEPEND(:,19) = (/ 3, 5, 19, 21 /)
SSADEPEND(:,20) = (/ 3, 5, 20, 21 /)
SSADEPEND(:,21) = (/ 2, 5, 21, 0 /)
SSADEPEND(:,22) = (/ 3, 6, 22, 24 /)
SSADEPEND(:,23) = (/ 3, 6, 23, 24 /)
SSADEPEND(:,24) = (/ 2, 6, 24, 0 /)

  end subroutine init_para

end module parameters
