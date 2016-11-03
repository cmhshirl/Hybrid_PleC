module parameters
  implicit none

  !************************
  !*                      *
  !* Parameters for system*
  !*                      *
  !************************
  integer, parameter:: MM = 50
  integer, parameter:: SNUM = 45
  integer, parameter:: SNUM_ODE = 8
  integer, parameter:: RNUM = 141
  integer, parameter:: DPNUM_ODE = 35
  integer, parameter:: DPNUM_SSA = 33

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

  integer, parameter :: mRNA_DivJ	= 9
  integer, parameter :: mRNA_DivK	= 10
  integer, parameter :: mRNA_PleC	= 11
  integer, parameter :: mRNA_DivL	= 12
  integer, parameter :: mRNA_CckA	= 13
  integer, parameter :: mRNA_CtrA	= 14
  
  integer, parameter :: DivL	    = 15
  integer, parameter :: DivLKp	    = 16
  integer, parameter :: CckA_phos	= 17
  integer, parameter :: CckA_kin	= 18
  
  integer, parameter :: DivJ_st		= 40
  integer, parameter :: DivL_st		= 41
  integer, parameter :: CckA_st		= 42
  integer, parameter :: PleC_st		= 43
  
  integer, parameter :: gene_ori	= 44
  integer, parameter :: gene_dup	= 45

  integer, parameter :: syn_OrimrDivJ = 30
  integer, parameter :: syn_DupmrDivJ = 31
  integer, parameter :: phos_DivK = 59
  integer, parameter :: deph_DivKp = 60

  !************************
  !*                      *
  !*Para for ODE RATE     *
  !*                      *
  !************************
  double precision, parameter :: SCALE = 1.0e+3
  double precision, parameter :: D_SP1           = 1.0e+1
  double precision, parameter :: D_SP2           = 1.0e+0  
!  double precision, parameter :: D_mRNA         = 5.0e-2
  double precision, parameter :: K_bndKp_DivL   = 2.5/SCALE
  double precision, parameter :: K_ubdKp_DivLKp = 5.0e-1
  double precision, parameter :: K_phos_CtrA    = 6.0e+2/SCALE
  double precision, parameter :: K_deph_CtrAp   = 6.0e+2/SCALE
  double precision, parameter :: K_ph2k_CckA    = 1.0e+1
  double precision, parameter :: K_deg_DivLKp   = 5.0e-2
  double precision, parameter :: K_syn_DivLf    = 1.25e-2*SCALE
  double precision, parameter :: K_deg_DivLf    = 5.0e-2
  double precision, parameter :: K_bnd_DivL     = 1.0e+0
  double precision, parameter :: K_ubd_DivL     = 1.0e-1
  double precision, parameter :: K_deg_DivL     = 5.0e-2

  double precision, parameter :: K_syn_DivJf = 1.25e+1
  double precision, parameter :: K_syn_CckA	 = 1.25e+1
  double precision, parameter :: K_syn_CtrA	 = 2.5e+1
  double precision, parameter :: K_deg_CtrA	 = 5.0e-2
  double precision, parameter :: K_deg_CtrAp = 5.0e-2
  double precision, parameter :: K_syn_DivK	 = 2.5e+1
  double precision, parameter :: K_syn_PleCf = 5.0e+1
  
  double precision, parameter :: K_deg_PleCf = 5.0e-2
  double precision, parameter :: K_k2ph_CckA = 1.0e+0
  
  double precision, parameter :: K_sticky	 = 1.0e+0
  double precision, parameter :: K_syn_mRNA	 = 6.25e-1  
  double precision, parameter :: K_deg_mRNA	 = 2.5e-1

  !**********************
  !*                    *
  !* changing variables *
  !*                    *
  !**********************
  double precision :: h
  double precision :: a0
  double precision, dimension(RNUM) :: a
  double precision, dimension(RNUM) :: totPro
  double precision, dimension(MM, RNUM) :: prePro
  double precision, dimension(MM) :: auxt
  double precision, dimension(SNUM) :: p
  double precision, dimension(MM, SNUM) :: pp
  double precision:: mu
  integer(KIND=8), dimension(10) :: fire

  !**********************
  !*                    *
  !* Param for SSA rate *
  !*                    *
  !**********************
  double precision, dimension(RNUM) :: RATE
  
  integer, dimension(DPNUM_ODE), parameter :: ODEDEPEND = (/ 1,	2,	3,	4,	5,	6,  7,	8,	9,	10,	11,	12,	13,	14,	15,	16,	17,	18,	19,	20,	21,	22,	23,	24,	25,	26,	27,	28,	29,	50,	63,	69,	70,	73,	83 /)
    
  integer, dimension(RSPACE, RNUM), parameter :: NETWORK = reshape( (/ &
10,	9,	0,	0,	0,	0, &
10,	10,	0,	0,	0,	0, &
10,	11,	0,	0,	0,	0, &
10,	12,	0,	0,	0,	0, &
10,	13,	0,	0,	0,	0, &
10,	14,	0,	0,	0,	0, &
8,	19,	2,	20,	0,	0, &
8,	19,	1,	21,	0,	0, &
8,	15,	1,	16,	0,	0, &
8,	22,	1,	23,	0,	0, &
8,	22,	2,	24,	0,	0, &
8,	23,	1,	25,	0,	0, &
8,	23,	2,	26,	0,	0, &
8,	24,	1,	26,	0,	0, &
8,	24,	2,	27,	0,	0, &
8,	29,	1,	28,	0,	0, &
8,	29,	2,	30,	0,	0, &
8,	31,	1,	30,	0,	0, &
8,	31,	2,	32,	0,	0, &
8,	34,	1,	29,	0,	0, &
8,	34,	2,	31,	0,	0, &
8,	33,	1,	36,	0,	0, &
8,	33,	1,	35,	0,	0, &
8,	33,	2,	37,	0,	0, &
8,	35,	1,	38,	0,	0, &
8,	37,	1,	39,	0,	0, &
7,	23,	2,	24,	1,	0, &
7,	24,	1,	23,	2,	0, &
9,	17,	0,	18,	0,	15, &
1,	0,	0,	9,	0,	44, &
1,	0,	0,	9,	0,	45, &
2,	9,	0,	0,	0,	0, &
1,	0,	0,	10,	0,	44, &
1,	0,	0,	10,	0,	45, &
2,	10,	0,	0,	0,	0, &
1,	0,	0,	11,	0,	44, &
1,	0,	0,	11,	0,	45, &
2,	11,	0,	0,	0,	0, &
1,	0,	0,	12,	0,	44, &
1,	0,	0,	12,	0,	45, &
2,	12,	0,	0,	0,	0, &
1,	0,	0,	13,	0,	44, &
1,	0,	0,	13,	0,	45, &
2,	13,	0,	0,	0,	0, &
1,	0,	0,	14,	0,	44, &
1,	0,	0,	14,	0,	45, &
2,	14,	0,	0,	0,	0, &
1,	0,	0,	3,	0,	9, &
2,	3,	0,	0,	0,	0, &
5,	3,	0,	19,	0,	40, &
3,	19,	0,	3,	0,	0, &
2,	19,	0,	0,	0,	0, &
4,	20,	0,	19,	2,	0, &
3,	20,	0,	21,	0,	0, &
3,	21,	0,	20,	0,	0, &
4,	21,	0,	19,	1,	0, &
2,	20,	0,	0,	0,	0, &
2,	21,	0,	0,	0,	0, &
3,	2,	0,	1,	0,	0, &
3,	1,	0,	2,	0,	0, &
1,	0,	0,	4,	0,	13, &
2,	4,	0,	0,	0,	0, &
5,	4,	0,	17,	0,	42, &
3,	17,	0,	4,	0,	0, &
3,	18,	0,	17,	0,	0, &
2,	18,	0,	0,	0,	0, &
2,	17,	0,	0,	0,	0, &
1,	0,	0,	5,	0,	14, &
2,	5,	0,	0,	0,	0, &
2,	6,	0,	0,	0,	0, &
1,	0,	0,	7,	0,	12, &
2,	7,	0,	0,	0,	0, &
5,	7,	0,	15,	0,	41, &
3,	15,	0,	7,	0,	0, &
2,	15,	0,	0,	0,	0, &
4,	16,	0,	15,	1,	0, &
2,	16,	0,	0,	0,	0, &
1,	0,	0,	2,	0,	10, &
2,	2,	0,	0,	0,	0, &
2,	1,	0,	0,	0,	0, &
1,	0,	0,	8,	0,	11, &
2,	8,	0,	0,	0,	0, &
5,	8,	0,	22,	0,	43, &
3,	22,	0,	8,	0,	0, &
4,	23,	0,	22,	1,	0, &
4,	24,	0,	22,	2,	0, &
3,	24,	0,	23,	0,	0, &
3,	23,	0,	24,	0,	0, &
4,	25,	0,	23,	1,	0, &
4,	26,	0,	24,	1,	0, &
4,	26,	0,	23,	2,	0, &
4,	27,	0,	24,	2,	0, &
3,	25,	0,	28,	0,	0, &
3,	26,	0,	30,	0,	0, &
3,	27,	0,	32,	0,	0, &
3,	28,	0,	25,	0,	0, &
3,	30,	0,	26,	0,	0, &
3,	32,	0,	27,	0,	0, &
4,	28,	0,	29,	1,	0, &
4,	30,	0,	31,	1,	0, &
4,	30,	0,	29,	2,	0, &
4,	32,	0,	31,	2,	0, &
4,	29,	0,	34,	1,	0, &
4,	31,	0,	34,	2,	0, &
3,	34,	0,	22,	0,	0, &
4,	36,	0,	33,	1,	0, &
4,	35,	0,	33,	1,	0, &
4,	37,	0,	33,	2,	0, &
3,	33,	0,	22,	0,	0, &
3,	31,	0,	36,	0,	0, &
3,	30,	0,	38,	0,	0, &
3,	32,	0,	39,	0,	0, &
3,	39,	0,	25,	0,	0, &
3,	36,	0,	31,	0,	0, &
3,	38,	0,	30,	0,	0, &
3,	39,	0,	32,	0,	0, &
3,	25,	0,	39,	0,	0, &
3,	29,	0,	35,	0,	0, &
3,	31,	0,	37,	0,	0, &
3,	35,	0,	29,	0,	0, &
3,	37,	0,	31,	0,	0, &
4,	38,	0,	35,	1,	0, &
4,	39,	0,	37,	1,	0, &
2,	22,	0,	0,	0,	0, &
3,	23,	0,	1,	0,	0, &
2,	24,	0,	0,	0,	0, &
4,	25,	0,	1,	1,	0, &
3,	26,	0,	1,	0,	0, &
2,	27,	0,	0,	0,	0, &
4,	28,	0,	1,	1,	0, &
3,	30,	0,	1,	0,	0, &
2,	32,	0,	0,	0,	0, &
3,	29,	0,	1,	0,	0, &
2,	31,	0,	0,	0,	0, &
2,	34,	0,	0,	0,	0, &
2,	33,	0,	0,	0,	0, &
2,	37,	0,	0,	0,	0, &
3,	36,	0,	1,	0,	0, &
3,	35,	0,	1,	0,	0, &
3,	39,	0,	1,	0,	0, &
4,	38,	0,	1,	1,	0 /), (/ RSPACE, RNUM /) )

  integer, dimension(DPNUM_SSA,RNUM) :: SSADEPEND


contains

  subroutine init_para()

RATE = (/ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  0.005,  0.005,  0.0025, &
0.005,  5e-05,  0.005,  1.6e-05,  0.0016, 1.6e-05,  0.005,  0.005,  0.005,  0.005, &
0.005,  0.005,  0.005,  0.005,  0.005,  0.005,  0.005,  1.6e-05,  0.0016, 10.0, 0.625, &
0.625,  0.25, 0.625,  0.625,  0.25, 0.625,  0.625,  0.25, 0.625,  0.625,  0.25, 0.625, &
0.625,  0.25, 6.25,  6.25,  0.25, 12.5, 0.05, 1.0,  0.0,  0.05, 0.0016, 5.0,  0.16, &
1.0,  0.05, 0.05, 0.0, 0.0, 12.5, 0.05, 1.0,  0.1,  1.0,  0.05, 0.05, 25.0, 0.05, 0.05, &
12.5, 0.05, 1.0,  0.1,  0.05, 0.5,  0.05, 25.0, 0.005,  0.005,  50.0, 0.05, 1.0,  0.5, &
5.0,  5.0,  0.005,  10.0, 2.5,  0.00016,  0.00016,  1.6e-08,  2.5,  5.0,  5.0,  5.0, &
5.0,  5.0,  0.16, 0.16, 0.0016, 0.0016, 0.16, 0.0016, 5.0,  0.16, 0.16, 0.0016, 5.0, &
5.0,  5.0,  5.0,  5.0,  0.16, 0.16, 0.16, 0.0755, 5.0,  5.0,  5.0,  5.0,  0.16, 0.16, &
0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 /)

SSADEPEND(:,1) = (/ 3, 1, 32, 48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,2) = (/ 3, 2, 35, 78, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,3) = (/ 3, 3, 38, 81, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,4) = (/ 3, 4, 41, 71, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,5) = (/ 3, 5, 44, 61, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,6) = (/ 3, 6, 47, 68, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,7) = (/ 17, 7, 8, 11, 13, 15, 17, 19, 21, 24, 27, 51, 52, 53, 54, 57, 59, 79, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,8) = (/ 21, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 25, 26, 28, 51, 52, 55, 56, 58, 60, 80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,9) = (/ 20, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 25, 26, 28, 29, 60, 74, 75, 76, 77, 80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,10) = (/ 23, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 22, 23, 25, 26, 27, 28, 60, 80, 84, 85, 88, 124, 125, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,11) = (/ 19, 7, 10, 11, 13, 14, 15, 17, 19, 21, 24, 27, 28, 59, 79, 84, 86, 87, 124, 126, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,12) = (/ 24, 8, 9, 10, 12, 13, 14, 16, 18, 20, 22, 23, 25, 26, 27, 28, 60, 80, 85, 88, 89, 93, 117, 125, 127, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,13) = (/ 19, 7, 11, 12, 13, 15, 17, 19, 21, 24, 27, 59, 79, 85, 88, 90, 91, 94, 125, 128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,14) = (/ 23, 8, 9, 10, 12, 14, 15, 16, 18, 20, 22, 23, 25, 26, 28, 60, 80, 86, 87, 90, 91, 94, 126, 128, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,15) = (/ 19, 7, 11, 13, 14, 15, 17, 19, 21, 24, 27, 28, 59, 79, 86, 87, 92, 95, 126, 129, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,16) = (/ 22, 8, 9, 10, 12, 14, 16, 17, 18, 20, 22, 23, 25, 26, 28, 60, 80, 96, 99, 103, 118, 130, 133, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,17) = (/ 20, 7, 11, 13, 15, 16, 17, 19, 21, 24, 27, 59, 79, 97, 100, 101, 103, 111, 118, 131, 133, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,18) = (/ 25, 8, 9, 10, 12, 14, 16, 18, 19, 20, 22, 23, 25, 26, 28, 60, 80, 97, 100, 101, 104, 110, 111, 119, 131, 134, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,19) = (/ 20, 7, 11, 13, 15, 17, 18, 19, 21, 24, 27, 59, 79, 98, 102, 104, 110, 112, 119, 132, 134, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,20) = (/ 22, 8, 9, 10, 12, 14, 16, 17, 18, 20, 21, 22, 23, 25, 26, 28, 60, 80, 103, 105, 118, 133, 135, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,21) = (/ 19, 7, 11, 13, 15, 17, 18, 19, 20, 21, 24, 27, 59, 79, 104, 105, 110, 119, 134, 135, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,22) = (/ 21, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 24, 25, 26, 28, 60, 80, 106, 109, 114, 136, 138, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,23) = (/ 21, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 24, 25, 26, 28, 60, 80, 107, 109, 120, 136, 139, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,24) = (/ 19, 7, 11, 13, 15, 17, 19, 21, 22, 23, 24, 26, 27, 59, 79, 108, 109, 121, 136, 137, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,25) = (/ 21, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 25, 26, 28, 60, 80, 107, 115, 120, 122, 139, 141, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,26) = (/ 22, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 25, 26, 28, 60, 80, 108, 113, 116, 121, 123, 137, 140, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,27) = (/ 32, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 59, 60, 79, 80, 85, 86, 87, 88, 125, 126 /)
SSADEPEND(:,28) = (/ 32, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 59, 60, 79, 80, 85, 86, 87, 88, 125, 126 /)
SSADEPEND(:,29) = (/ 5, 29, 64, 65, 66, 67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,30) = (/ 4, 1, 30, 32, 48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,31) = (/ 4, 1, 31, 32, 48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,32) = (/ 3, 1, 32, 48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,33) = (/ 4, 2, 33, 35, 78, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,34) = (/ 4, 2, 34, 35, 78, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,35) = (/ 3, 2, 35, 78, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,36) = (/ 4, 3, 36, 38, 81, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,37) = (/ 4, 3, 37, 38, 81, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,38) = (/ 3, 3, 38, 81, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,39) = (/ 4, 4, 39, 41, 71, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,40) = (/ 4, 4, 40, 41, 71, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,41) = (/ 3, 4, 41, 71, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,42) = (/ 4, 5, 42, 44, 61, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,43) = (/ 4, 5, 43, 44, 61, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,44) = (/ 3, 5, 44, 61, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,45) = (/ 4, 6, 45, 47, 68, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,46) = (/ 4, 6, 46, 47, 68, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,47) = (/ 3, 6, 47, 68, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,48) = (/ 3, 48, 49, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,49) = (/ 2, 49, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,50) = (/ 6, 7, 8, 49, 50, 51, 52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,51) = (/ 6, 7, 8, 49, 50, 51, 52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,52) = (/ 4, 7, 8, 51, 52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,53) = (/ 17, 7, 8, 11, 13, 15, 17, 19, 21, 24, 27, 51, 52, 53, 54, 57, 59, 79, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,54) = (/ 6, 53, 54, 55, 56, 57, 58, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,55) = (/ 6, 53, 54, 55, 56, 57, 58, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,56) = (/ 21, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 25, 26, 28, 51, 52, 55, 56, 58, 60, 80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,57) = (/ 3, 53, 54, 57, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,58) = (/ 3, 55, 56, 58, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,59) = (/ 26, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 59, 60, 79, 80, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,60) = (/ 26, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 59, 60, 79, 80, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,61) = (/ 3, 61, 62, 63, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,62) = (/ 2, 62, 63, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,63) = (/ 5, 29, 62, 63, 64, 67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,64) = (/ 5, 29, 62, 63, 64, 67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,65) = (/ 5, 29, 64, 65, 66, 67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,66) = (/ 2, 65, 66, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,67) = (/ 3, 29, 64, 67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,68) = (/ 2, 68, 69, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,69) = (/ 1, 69, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,70) = (/ 1, 70, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,71) = (/ 3, 71, 72, 73, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,72) = (/ 2, 72, 73, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,73) = (/ 6, 9, 29, 72, 73, 74, 75, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,74) = (/ 6, 9, 29, 72, 73, 74, 75, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,75) = (/ 4, 9, 29, 74, 75, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,76) = (/ 20, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 25, 26, 28, 29, 60, 74, 75, 76, 77, 80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,77) = (/ 2, 76, 77, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,78) = (/ 12, 7, 11, 13, 15, 17, 19, 21, 24, 27, 59, 78, 79, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,79) = (/ 11, 7, 11, 13, 15, 17, 19, 21, 24, 27, 59, 79, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,80) = (/ 15, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 25, 26, 28, 60, 80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,81) = (/ 3, 81, 82, 83, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,82) = (/ 2, 82, 83, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,83) = (/ 6, 10, 11, 82, 83, 84, 124, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,84) = (/ 6, 10, 11, 82, 83, 84, 124, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,85) = (/ 23, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 22, 23, 25, 26, 27, 28, 60, 80, 84, 85, 88, 124, 125, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,86) = (/ 19, 7, 10, 11, 13, 14, 15, 17, 19, 21, 24, 27, 28, 59, 79, 84, 86, 87, 124, 126, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,87) = (/ 12, 12, 13, 14, 15, 27, 28, 85, 86, 87, 88, 125, 126, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,88) = (/ 12, 12, 13, 14, 15, 27, 28, 85, 86, 87, 88, 125, 126, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,89) = (/ 24, 8, 9, 10, 12, 13, 14, 16, 18, 20, 22, 23, 25, 26, 27, 28, 60, 80, 85, 88, 89, 93, 117, 125, 127, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,90) = (/ 23, 8, 9, 10, 12, 14, 15, 16, 18, 20, 22, 23, 25, 26, 28, 60, 80, 86, 87, 90, 91, 94, 126, 128, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,91) = (/ 19, 7, 11, 12, 13, 15, 17, 19, 21, 24, 27, 59, 79, 85, 88, 90, 91, 94, 125, 128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,92) = (/ 19, 7, 11, 13, 14, 15, 17, 19, 21, 24, 27, 28, 59, 79, 86, 87, 92, 95, 126, 129, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,93) = (/ 7, 89, 93, 96, 99, 117, 127, 130, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,94) = (/ 9, 90, 91, 94, 97, 100, 101, 111, 128, 131, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,95) = (/ 7, 92, 95, 98, 102, 112, 129, 132, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,96) = (/ 7, 89, 93, 96, 99, 117, 127, 130, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,97) = (/ 9, 90, 91, 94, 97, 100, 101, 111, 128, 131, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,98) = (/ 7, 92, 95, 98, 102, 112, 129, 132, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,99) = (/ 22, 8, 9, 10, 12, 14, 16, 17, 18, 20, 22, 23, 25, 26, 28, 60, 80, 96, 99, 103, 118, 130, 133, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,100) = (/ 25, 8, 9, 10, 12, 14, 16, 18, 19, 20, 22, 23, 25, 26, 28, 60, 80, 97, 100, 101, 104, 110, 111, 119, 131, 134, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,101) = (/ 20, 7, 11, 13, 15, 16, 17, 19, 21, 24, 27, 59, 79, 97, 100, 101, 103, 111, 118, 131, 133, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,102) = (/ 20, 7, 11, 13, 15, 17, 18, 19, 21, 24, 27, 59, 79, 98, 102, 104, 110, 112, 119, 132, 134, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,103) = (/ 22, 8, 9, 10, 12, 14, 16, 17, 18, 20, 21, 22, 23, 25, 26, 28, 60, 80, 103, 105, 118, 133, 135, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,104) = (/ 19, 7, 11, 13, 15, 17, 18, 19, 20, 21, 24, 27, 59, 79, 104, 105, 110, 119, 134, 135, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,105) = (/ 8, 10, 11, 20, 21, 84, 105, 124, 135, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,106) = (/ 21, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 24, 25, 26, 28, 60, 80, 106, 109, 114, 136, 138, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,107) = (/ 21, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 24, 25, 26, 28, 60, 80, 107, 109, 120, 136, 139, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,108) = (/ 19, 7, 11, 13, 15, 17, 19, 21, 22, 23, 24, 26, 27, 59, 79, 108, 109, 121, 136, 137, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,109) = (/ 9, 10, 11, 22, 23, 24, 84, 109, 124, 136, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,110) = (/ 9, 18, 19, 104, 106, 110, 114, 119, 134, 138, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,111) = (/ 8, 97, 100, 101, 111, 115, 122, 131, 141, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,112) = (/ 8, 98, 102, 112, 113, 116, 123, 132, 140, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,113) = (/ 8, 89, 93, 113, 116, 117, 123, 127, 140, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,114) = (/ 9, 18, 19, 104, 106, 110, 114, 119, 134, 138, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,115) = (/ 8, 97, 100, 101, 111, 115, 122, 131, 141, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,116) = (/ 8, 98, 102, 112, 113, 116, 123, 132, 140, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,117) = (/ 8, 89, 93, 113, 116, 117, 123, 127, 140, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,118) = (/ 9, 16, 17, 25, 103, 107, 118, 120, 133, 139, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,119) = (/ 10, 18, 19, 26, 104, 108, 110, 119, 121, 134, 137, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,120) = (/ 9, 16, 17, 25, 103, 107, 118, 120, 133, 139, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,121) = (/ 10, 18, 19, 26, 104, 108, 110, 119, 121, 134, 137, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,122) = (/ 21, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 25, 26, 28, 60, 80, 107, 115, 120, 122, 139, 141, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,123) = (/ 22, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 25, 26, 28, 60, 80, 108, 113, 116, 121, 123, 137, 140, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,124) = (/ 4, 10, 11, 84, 124, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,125) = (/ 20, 8, 9, 10, 12, 13, 14, 16, 18, 20, 22, 23, 25, 26, 27, 28, 60, 80, 85, 88, 125, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,126) = (/ 6, 14, 15, 28, 86, 87, 126, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,127) = (/ 19, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 25, 26, 28, 60, 80, 89, 93, 117, 127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,128) = (/ 19, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 25, 26, 28, 60, 80, 90, 91, 94, 128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,129) = (/ 3, 92, 95, 129, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,130) = (/ 18, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 25, 26, 28, 60, 80, 96, 99, 130, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,131) = (/ 20, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 25, 26, 28, 60, 80, 97, 100, 101, 111, 131, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,132) = (/ 4, 98, 102, 112, 132, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,133) = (/ 19, 8, 9, 10, 12, 14, 16, 17, 18, 20, 22, 23, 25, 26, 28, 60, 80, 103, 118, 133, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,134) = (/ 6, 18, 19, 104, 110, 119, 134, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,135) = (/ 4, 20, 21, 105, 135, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,136) = (/ 5, 22, 23, 24, 109, 136, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,137) = (/ 4, 26, 108, 121, 137, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,138) = (/ 18, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 25, 26, 28, 60, 80, 106, 114, 138, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,139) = (/ 18, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 25, 26, 28, 60, 80, 107, 120, 139, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,140) = (/ 19, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 25, 26, 28, 60, 80, 113, 116, 123, 140, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
SSADEPEND(:,141) = (/ 18, 8, 9, 10, 12, 14, 16, 18, 20, 22, 23, 25, 26, 28, 60, 80, 115, 122, 141, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
  end subroutine init_para

end module parameters
