# define IJth(y,i,j) y((i-1)*5+j)

module stochastic
  use deterministic


contains

  !****************************
  !*                          *
  !* ssa System of Cell Cycle *
  !*                          *
  !****************************
  subroutine ssa_alg()
    implicit none

    integer :: i, ita, bin_ita, inf_ita, j,k
    double precision :: sum, sum_a
    double precision:: random, propensity
    
    !***************************************
    !choose reaction channel
    !***************************************
    call random_number(random)
    sum = a0 * random
    
    ita=1
    sum_a=a(ita)
    do while(sum_a<sum)
      ita = ita + 1
      sum_a = sum_a + a(ita)
    end do
    
    !***************************************
    !update population
    !***************************************
    call population_update(bin_ita, ita)


do i=1,SNUM
    if (p(i) .LT. 0) then
        print *, 'ita', ita
        print *, 'a', a(ita)
        print *, 'sum', a0,random, sum
        print *, 'i ', i
        print *, 'p ', p(i)
        do j=1, MM
            print *, 'pp', pp(j,i)
        end do

do k=1,RNUM
print *, k, 'a(i)', a(k)
end do
        stop -1
    end if
end do
    
    do i=1, SSADEPEND(1,ita)
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

    do i=1, SNUM_ODE
        p(i) = 0.0
        do j=1, MM
            pp(j,i) = IJth(y,i,j)
            p(i) = p(i) + pp(j,i)
        end do
    end do

    h = 1.3/MM
    do i=SNUM_ODE+1, SNUM
        call initZero(i)
    end do

    do i=MM*4/5+1, MM
        pp(i,PleC_st) = 1.0
    end do
    p(PleC_st) = MM/5.0

    call uniFill(DivL_st)
    call uniFill(CckA_st)

  end subroutine init_ssa

  subroutine DNAorigination()
      call initZero(gene_ori)
      call initZero(gene_dup)
      p(gene_ori) = 1
      pp(MM*4/5, gene_ori) = 1
  end subroutine DNAorigination
  
  subroutine DNAreplication()
      p(gene_dup) = 1
      pp(MM/5, gene_dup) = 1
  end subroutine DNAreplication
  
  subroutine introDivJ()
    implicit none
    integer :: i

    call initZero(DivJ_st)
    do i=MM*4/5+1, MM
      pp(i,DivJ_st) = 1.0
    end do
    p(DivJ_st) = MM/5.0
  
  end subroutine introDivJ
  
  
  subroutine uniFill(ita)
    implicit none
    integer :: i
    integer, intent(in) :: ita

    do i=1, MM
      pp(i,ita) = 1
    end do
    p(ita) = MM
  
  end subroutine uniFill
  
  subroutine initZero(ita)
    implicit none
    integer :: i
    integer, intent(in) :: ita

    do i=1, MM
      pp(i,ita) = 0.0
    end do 
    p(ita) = 0.0
  
  end subroutine initZero
  
  subroutine clearPleC()
    implicit none
    integer :: i

    call initZero(PleC_st)
    call initZero(CckA_st)
    do i=MM*4/5+1, MM
      pp(i,CckA_st) = 1.0
    end do 
    p(CckA_st) = MM/5.0
  
  end subroutine clearPleC
  
  subroutine checkPoint3()
    implicit none
    integer :: i

    call initZero(PleC_st) 
    call initZero(DivL_st) 
    call initZero(CckA_st) 
    
    do i=1, MM/5
      pp(i,PleC_st) = 1.0
      pp(i,DivL_st) = 1.0
      pp(i,CckA_st) = 1.0
    end do
    p(PleC_st) = MM/5.0
    p(DivL_st) = MM/5.0
      
    do i=MM*4/5+1, MM
      pp(i,CckA_st) = 1.0
    end do 
    p(CckA_st) = 2*MM/5.0
  
  end subroutine checkPoint3



  !/* update propensity influneced by SSA */
  subroutine propensity_update_ssa(bin_ita, ita, propensity)
    implicit none
    integer :: react_type, react1_ita, react2_ita, aux_ita
    double precision :: rt, temp
    integer, intent(in) :: ita, bin_ita
    double precision, intent(out) :: propensity

    react_type = NETWORK(RTYPE, ita)
    react1_ita = NETWORK(REACT1, ita)
    react2_ita = NETWORK(REACT2, ita)
    aux_ita = NETWORK(AUX, ita)
    rt = RATE(ita)
    temp = 0.0

    select case (react_type)
            
        case (syn)
            propensity= rt * p(aux_ita)
            
        case (deg : decomposition)
            propensity= rt * p(react1_ita)

        case (sticky)
            totPro(ita) = totPro(ita) - prePro(bin_ita, ita)
            prePro(bin_ita, ita) = pp(bin_ita, react1_ita)*pp(bin_ita, aux_ita)
            totPro(ita) = totPro(ita) + prePro(bin_ita, ita)
            propensity= rt*totPro(ita)

        case (displacement : combination)
            totPro(ita) = totPro(ita) - prePro(bin_ita, ita)
            prePro(bin_ita, ita) = pp(bin_ita, react1_ita)*pp(bin_ita, react2_ita)
            totPro(ita) = totPro(ita) + prePro(bin_ita, ita)
            propensity= rt*totPro(ita)/h
            
        case default
            print *, "Reaction_type in propensity_update_ssa: ", react_type, " not found "
            stop -1
    end select

  end subroutine propensity_update_ssa


  subroutine population_update(bin_ita, ita)
    implicit none
    integer :: react_type, react1_ita, react2_ita, prod1_ita, prod2_ita, aux_ita
    double precision :: ADD = 1.0, REMOVE = -1.0
    integer, intent(in) :: ita
    integer, intent(out) :: bin_ita 

    react_type = NETWORK(RTYPE, ita)
    react1_ita = NETWORK(REACT1, ita)
    react2_ita = NETWORK(REACT2, ita)
    prod1_ita = NETWORK(PROD1, ita)
    prod2_ita = NETWORK(PROD2, ita)
    aux_ita = NETWORK(AUX, ita)

!print *, 'react_type: ', react_type
    select case (react_type)
            
      case (syn) ! null -> A aux
        call selectBin1st(aux_ita, bin_ita)
        p(prod1_ita) = p(prod1_ita) + ADD
        pp(bin_ita, prod1_ita) = pp(bin_ita, prod1_ita) + ADD
        if (prod1_ita<SNUM_ODE) then 
          IJth(y, prod1_ita, bin_ita) = IJth(y, prod1_ita, bin_ita) + ADD
        end if


      case (deg) ! A -> null
        call selectBin1st(bin_ita, react1_ita)
        p(react1_ita) = p(react1_ita) + REMOVE
        pp(bin_ita, react1_ita) = pp(bin_ita, react1_ita) + REMOVE
        if (react1_ita<SNUM_ODE) then
          IJth(y, react1_ita, bin_ita) = IJth(y, react1_ita, bin_ita) + REMOVE
        end if
            
      case (first_order) ! A -> B
        call selectBin1st(bin_ita, react1_ita)
        p(react1_ita) = p(react1_ita) + REMOVE
        pp(bin_ita, react1_ita) = pp(bin_ita, react1_ita) + REMOVE
        p(prod1_ita) = p(prod1_ita) + ADD
        pp(bin_ita, prod1_ita) = pp(bin_ita, prod1_ita) + ADD
        if (react1_ita<SNUM_ODE) then
          IJth(y, react1_ita, bin_ita) = IJth(y, react1_ita, bin_ita) + REMOVE
        end if
        if (prod1_ita<SNUM_ODE) then
          IJth(y, prod1_ita, bin_ita) = IJth(y, prod1_ita, bin_ita) + ADD
        end if
    
        
      case (decomposition) ! A -> C + D
        call selectBin1st(bin_ita, react1_ita)
        p(react1_ita) = p(react1_ita) + REMOVE
        p(prod1_ita) = p(prod1_ita) + ADD
        p(prod2_ita) = p(prod2_ita) + ADD
        pp(bin_ita, react1_ita) = pp(bin_ita, react1_ita) + REMOVE
        pp(bin_ita, prod1_ita) = pp(bin_ita, prod1_ita) + ADD
        pp(bin_ita, prod2_ita) = pp(bin_ita, prod2_ita) + ADD
        if (react1_ita<SNUM_ODE) then
          IJth(y, react1_ita, bin_ita) = IJth(y, react1_ita, bin_ita) + REMOVE
        end if
        if (prod1_ita<SNUM_ODE) then
          IJth(y, prod1_ita, bin_ita) = IJth(y, prod1_ita, bin_ita) + ADD
        end if
        if (prod2_ita<SNUM_ODE) then
          IJth(y, prod2_ita, bin_ita) = IJth(y, prod2_ita, bin_ita) + ADD
        end if
  
      case (sticky) ! A --> C aux
        call selectBin2nd(bin_ita, ita)
        p(react1_ita) = p(react1_ita) + REMOVE
        p(prod1_ita) = p(prod1_ita) + ADD
        pp(bin_ita, react1_ita) = pp(bin_ita, react1_ita) + REMOVE
        pp(bin_ita, prod1_ita) = pp(bin_ita, prod1_ita) + ADD
        if (react1_ita<SNUM_ODE) then
          IJth(y, react1_ita, bin_ita) = IJth(y, react1_ita, bin_ita) + REMOVE
        end if
        if (prod1_ita<SNUM_ODE) then
          IJth(y, prod1_ita, bin_ita) = IJth(y, prod1_ita, bin_ita) + ADD
        end if
          
      case (combination) ! A + B -> C
        call selectBin2nd(bin_ita, ita)
        p(react1_ita) = p(react1_ita) + REMOVE
        p(react2_ita) = p(react2_ita) + REMOVE
        p(prod1_ita) = p(prod1_ita) + ADD
        pp(bin_ita, react1_ita) = pp(bin_ita, react1_ita) + REMOVE
        pp(bin_ita, react2_ita) = pp(bin_ita, react2_ita) + REMOVE
        pp(bin_ita, prod1_ita) = pp(bin_ita, prod1_ita) + ADD
        if (react1_ita<SNUM_ODE) then
          IJth(y, react1_ita, bin_ita) = IJth(y, react1_ita, bin_ita) + REMOVE
        end if
        if (react2_ita<SNUM_ODE) then
          IJth(y, react2_ita, bin_ita) = IJth(y, react2_ita, bin_ita) + REMOVE
        end if
        if (prod1_ita<SNUM_ODE) then
          IJth(y, prod1_ita, bin_ita) = IJth(y, prod1_ita, bin_ita) + ADD
        end if
       
      case (displacement) ! A + B -> C + D
        call selectBin2nd(bin_ita, ita)
        p(react1_ita) = p(react1_ita) + REMOVE
        p(react2_ita) = p(react2_ita) + REMOVE
        p(prod1_ita) = p(prod1_ita) + ADD
        p(prod2_ita) = p(prod2_ita) + ADD
        pp(bin_ita, react1_ita) = pp(bin_ita, react1_ita) + REMOVE
        pp(bin_ita, react2_ita) = pp(bin_ita, react2_ita) + REMOVE
        pp(bin_ita, prod1_ita) = pp(bin_ita, prod1_ita) + ADD
        pp(bin_ita, prod2_ita) = pp(bin_ita, prod2_ita) + ADD
        if (react1_ita<SNUM_ODE) then
          IJth(y, react1_ita, bin_ita) = IJth(y, react1_ita, bin_ita) + REMOVE
        end if
        if (react2_ita<SNUM_ODE) then
          IJth(y, react2_ita, bin_ita) = IJth(y, react2_ita, bin_ita) + REMOVE
        end if
        if (prod1_ita<SNUM_ODE) then
          IJth(y, prod1_ita, bin_ita) = IJth(y, prod1_ita, bin_ita) + ADD
        end if
        if (prod2_ita<SNUM_ODE) then
          IJth(y, prod2_ita, bin_ita) = IJth(y, prod2_ita, bin_ita) + ADD
        end if
    
      case default
        print *, "Reaction type in population_update: ", react_type, " not found !"
        stop -1
    
    end select

  end subroutine population_update


!/* select which bin fires the 1st order reaction */

  subroutine selectBin1st(bin_ita, ita)
    implicit none
    integer, intent(in) :: ita
    integer, intent(out) :: bin_ita
    double precision :: bin_ra0, bin_sum, random
    
    call random_number(random)
    bin_ra0 = random*p(ita)

    bin_ita = 1
    bin_sum = pp(bin_ita, ita)
    do while(bin_sum<bin_ra0) 
      bin_ita = bin_ita + 1

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


  !/* select which bin fires the 2nd order reaction */
  subroutine selectBin2nd(bin_ita, ita)
    implicit none
    integer, intent(in) :: ita
    integer, intent(out) :: bin_ita
    double precision :: bin_ra0, bin_sum, random
    
    call random_number(random)
    bin_ra0 = random*totPro(ita)

    bin_ita = 1
    bin_sum = prePro(bin_ita, ita)
    do while(bin_sum<bin_ra0) 
      bin_ita = bin_ita + 1
      bin_sum = bin_sum + prePro(bin_ita, ita)
    end do
    
    if (bin_ita>MM) then
        print *, "selectBin2nd exceed"
        stop -1
    end if 
    
 end subroutine selectBin2nd



end module stochastic

