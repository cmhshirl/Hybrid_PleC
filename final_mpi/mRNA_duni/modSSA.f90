# define IJth(y,i,j) y((i-1)*50+j)

module stochastic
  use deterministic
  integer::bin_diff

contains

  !****************************
  !*                          *
  !* ssa System of Cell Cycle *
  !*                          *
  !****************************
  subroutine ssa_alg()
    implicit none

    integer :: i, ita, bin_ita, inf_ita, j,k, react1_ita, react2_ita, aux_ita
    double precision :: sum, sum_a
    double precision:: random, propensity, pFalse
    
    !***************************************
    !choose reaction channel
    !***************************************
    random = DUNI()
    sum = a0 * random
    !print *, 'sum:',sum
    ita=1
    sum_a=a(ita)
    do while(sum_a<sum)
      ita = ita + 1

      !print error------------------------------
      if (ita>RNUM) then
        print *, 'a0:', a0
        do k=1,RNUM
            call cal_propensity_ssa(k, pFalse)

            if (ABS(pFalse-a(k)) .GT. 1.0e-6) then
                print *, k, 'a(i):', a(k), pFalse
                react1_ita = NETWORK(REACT1, k)
                react2_ita = NETWORK(REACT2, k)
                aux_ita = NETWORK(AUX, k)

                if (react1_ita>0) then
                    print *, 'p:', p(react1_ita)
                    do j= 1, MM
                        print *, j, 'pp(j):', pp(j,react1_ita)
                    end do
                end if

                if (react2_ita>0) then
                    print *, 'react2, p:', p(react1_ita),p(react2_ita)
                    do j= 1, MM
                        print *, j, 'pp(j):', pp(j,react1_ita),pp(j,react2_ita)
                    end do
                end if
    
                if (aux_ita>0) then
                    print *, 'aux, p:', p(react1_ita),p(aux_ita)
                    do j= 1, MM
                        print *, j, 'pp(j):',pp(j,react1_ita),pp(j,aux_ita)  
                    end do
                end if

            end if
        end do
      end if

      sum_a = sum_a + a(ita)
    end do

    !***************************************
    !update population
    !***************************************
    call population_update(bin_ita, ita)

!print *, 'CtrA, p ', p(5), 'CtrAp, p ', p(6)
!print *, 'phos:',a(7),p(5),p(18),'dephos:',a(8),p(6),p(17)
!print log
!do i=1,SNUM
!   if (p(i) .LT. 0) then
!        print *, 'ita', ita
!        print *, 'a', a(ita)
!        print *, 'sum', a0,random, sum
!        print *, 'i ', i, 'p ', p(i)
!        do j=1, MM
!            print *, 'pp', pp(j,i)
!        end do
!
!do k=1,RNUM
!print *, k, 'a(i)', a(k)
!end do
!        stop -1
!    end if
!end do

    do i=2, SSADEPEND(1,ita)+1
      !cout<<", i:"<<i<<", ssadepend:"<<SSADEPEND(ita,i)<<", type:"<<NETWORK(ita,TYPE)<<endl
      inf_ita = SSADEPEND(i,ita)
      a0 = a0 - a(inf_ita)
      call propensity_update_ssa(bin_ita, inf_ita, propensity)
      a(inf_ita) = propensity
      a0 = a0 + a(inf_ita)
    end do
  end subroutine ssa_alg


  !****************************
  !*                          *
  !* Initialize ssa           *
  !*                          *
  !****************************
  subroutine init_ssa()
    implicit none
    integer :: i, j

    h = 1.3/MM
    do i=1, SNUM_SSA
        call initZero(i, 0)
    end do

    do i=MM*9/10+1, MM
        IJth(y, PleC_s, i) = 1.0
    end do

    call uniFill(DivL_s)
    call uniFill(CckA_s)
  end subroutine init_ssa

  subroutine DNAorigination()
      call initZero(gene_ori, 0)
      call initZero(gene_dup, 0)
      p(gene_ori) = 1
      pp(MM*9/10, gene_ori) = 1
  end subroutine DNAorigination
  
  subroutine DNAreplication()
      p(gene_dup) = 1
      pp(MM/10, gene_dup) = 1
  end subroutine DNAreplication
 
 subroutine cellFixed()
    mu = 0.0
    RATE(syn_OrimrDivJ) = 0.0
    RATE(syn_DupmrDivJ) = 0.0
  end subroutine cellFixed

  subroutine cellGrow()
    mu = 0.005577
    RATE(syn_OrimrDivJ) = 0.625
    RATE(syn_DupmrDivJ) = 0.625
  end subroutine cellGrow

  subroutine introDivJ()
    implicit none
    integer :: i

    call initZero(DivJ_s, 1)
    do i=MM*9/10+1, MM
      IJth(y, DivJ_s, i) = 1.0
    end do
  
    Kphos_DivK = 0.05
    Kdeph_DivKp = 0.01
  end subroutine introDivJ
  
  
  subroutine uniFill(ita)
    implicit none
    integer :: i
    integer, intent(in) :: ita

    do i=1, MM
      IJth(y, ita, i) = 1.0
    end do

  end subroutine uniFill
  
  subroutine initZero(ita, flag)
    implicit none
    integer :: i
    integer, intent(in) :: ita, flag

    if (flag .EQ. 0) then! initialize SSA species population
      do i=1, MM
        pp(i,ita) = 0.0
      end do 
      p(ita) = 0.0
    else ! initialize ODE species population
      do i=1, MM
        IJth(y, ita, i) = 0.0
      end do 
    end if
  end subroutine initZero
  
  subroutine clearFire()
    implicit none
    integer :: i

    do i=1, 10
      fire(i) = 0
    end do 
  end subroutine clearFire
  
  subroutine clearPleC()
    implicit none
    integer :: i

    call initZero(PleC_s, 1)
    call initZero(CckA_s, 1)
    do i=MM*9/10+1, MM
      IJth(y, CckA_s, i) = 1.0
    end do 
  end subroutine clearPleC
  
  subroutine checkPoint3()
    implicit none
    integer :: i

    call initZero(PleC_s, 1) 
    call initZero(DivL_s, 1) 
    call initZero(CckA_s, 1) 
    
    do i=1, MM/10
      IJth(y, PleC_s, i) = 1.0
      IJth(y, DivL_s, i) = 1.0
      IJth(y, CckA_s, i) = 1.0
    end do

      
    do i=MM*9/10+1, MM
      IJth(y, CckA_s, i) = 1.0
    end do 
    
  end subroutine checkPoint3


  !/* update propensity influneced by SSA */
  subroutine propensity_update_ssa(bin_ita, ita, propensity)
    implicit none
    integer :: i, j
    integer :: react_type, react1_ita, react2_ita, aux_ita
    double precision :: rt, temp, Kmdl, avg, tempH, propH
    integer, intent(in) :: ita, bin_ita
    double precision, intent(out) :: propensity

    react_type = NETWORK(RTYPE, ita)
    react1_ita = NETWORK(REACT1, ita)
    react2_ita = NETWORK(REACT2, ita)
    aux_ita = NETWORK(AUX, ita)
    rt = RATE(ita)
    temp = 0.0

    select case (react_type)

        case (diff)
            propensity= 2*rt/(h*h)*p(react1_ita)

        case (syn)
            propensity= rt * p(aux_ita)
            
        case (deg)
            propensity= rt * p(react1_ita)

        case default
            print *, "Reaction_type in propensity_update_ssa: ", react_type, " not found "
            stop -1

    end select

  end subroutine propensity_update_ssa


  subroutine population_update(bin_ita, ita)
    implicit none
    integer :: react_type, react1_ita, react2_ita, prod1_ita, prod2_ita, aux_ita
    double precision :: ADD = 1.0, REMOVE = -1.0, random
    integer, intent(in) :: ita
    integer, intent(out) :: bin_ita 
    integer:: i,j

    react_type = NETWORK(RTYPE, ita)
    react1_ita = NETWORK(REACT1, ita)
    react2_ita = NETWORK(REACT2, ita)
    prod1_ita = NETWORK(PROD1, ita)
    prod2_ita = NETWORK(PROD2, ita)
    aux_ita = NETWORK(AUX, ita)
    bin_diff = 0
    fire(react_type) = fire(react_type) + 1
!print *, 'react_type',react_type,'react1_ita',react1_ita,'prod1_ita',prod1_ita,'aux_ita',aux_ita

    select case (react_type)

      case (diff) !diffusion
        call selectBin1st(bin_ita, react1_ita)
        random = DUNI()
        if (random>0.5) then
            if (bin_ita>1) then
                bin_diff = bin_ita-1
            end if
        else
            if (bin_ita<MM) then
                bin_diff = bin_ita+1
            end if
        end if
        if (bin_diff .NE. 0) then
            pp(bin_ita, react1_ita) = pp(bin_ita, react1_ita) + REMOVE
            pp(bin_diff,react1_ita) = pp(bin_diff,react1_ita) + ADD
        end if

      case (syn) ! null -> A aux
        call selectBin1st(bin_ita, aux_ita)
        p(prod1_ita) = p(prod1_ita) + ADD
        pp(bin_ita, prod1_ita) = pp(bin_ita, prod1_ita) + ADD

      case (deg) ! A -> null
        call selectBin1st(bin_ita, react1_ita)
        p(react1_ita) = p(react1_ita) + REMOVE
        pp(bin_ita, react1_ita) = pp(bin_ita, react1_ita) + REMOVE

      case default
        print *, "Reaction type in population_update: ", react_type, " not found !"
        stop -1
    
    end select

  end subroutine population_update

  subroutine cal_propensity_ssa(ita, propensity)
  
    implicit none
    integer, intent(in) :: ita
    double precision, intent(out) :: propensity
    integer :: i, j
    integer :: react_type, react1_ita, react2_ita, aux_ita
    double precision :: rt, temp, Kmdl, avg
    
    react_type = NETWORK(RTYPE, ita)
    react1_ita = NETWORK(REACT1, ita)
    react2_ita = NETWORK(REACT2, ita)
    aux_ita = NETWORK(AUX, ita)
    rt = RATE(ita)
    temp = 0.0
    
    select case (react_type)

      case (diff)
        propensity= 2*rt/(h*h)*p(react1_ita)

      case (syn)
        propensity= rt * p(aux_ita)

      case (deg)
        propensity= rt * p(react1_ita)

      case default
        print *, "Reaction_type in ita: ", ita, " not found !"
        print *, "Reaction_type in cal-propensity: ", react_type, " not found !"
        stop -1
    
    end select
  
  
  end subroutine cal_propensity_ssa

!/* select which bin fires the 1st order reaction */

  subroutine selectBin1st(bin_ita, ita)
    implicit none
    integer, intent(in) :: ita
    integer, intent(out) :: bin_ita
    double precision :: bin_ra0, bin_sum, random

    random = DUNI()
    bin_ra0 = random*p(ita)

    bin_ita = 1
    bin_sum = pp(bin_ita, ita)
    do while(bin_sum<bin_ra0)
      bin_ita = bin_ita + 1

!print log
if (bin_ita>MM) then
print *, "selectBin1st exceed"
print *, 'ita ', ita

print *, 'random ', random, 'p', p(ita)
print *, pp(:,ita)
stop -1
end if

      bin_sum = bin_sum + pp(bin_ita, ita)
    end do

    if (bin_ita>MM) then
        print *, "selectBin1st exceed"
        stop -1
    end if

   end subroutine selectBin1st



end module stochastic

