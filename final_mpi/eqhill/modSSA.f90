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
    double precision:: random, propensity, pFalse, pTemp
    
    !***************************************
    !choose reaction channel
    !***************************************
    call random_number(random)
    sum = a0 * random
    !print *, 'sum:',sum
    ita=1
    sum_a=a(ita)
    do while(sum_a<sum)
      ita = ita + 1

      !print error------------------------------
      if (ita>RNUM) then
        pTemp = 0.0
        do k=1,RNUM
          print *, 'check: ', k, 'a(i):', a(k)
          pTemp = pTemp + a(k)
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

        if (ABS(pTemp-a0) .GT. 1.0e-6) then
            print *, 'a0: ', a0, 'aSum: ', pTemp
        end if

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

    do i=MM*9/10+1, MM
        pp(i,PleC_st) = 1.0
    end do
    p(PleC_st) = MM/10.0

    call uniFill(DivL_st)
    call uniFill(CckA_st)
  end subroutine init_ssa

  subroutine DNAorigination()
      call initZero(gene_ori)
      call initZero(gene_dup)
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

    call initZero(DivJ_st)
    do i=MM*9/10+1, MM
      pp(i,DivJ_st) = 1.0
    end do
    p(DivJ_st) = MM/10.0
  
    RATE(phos_DivK) = 0.05
    RATE(deph_DivKp) = 0.01
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

    call initZero(PleC_st)
    call initZero(CckA_st)
    do i=MM*9/10+1, MM
      pp(i,CckA_st) = 1.0
    end do 
    p(CckA_st) = MM/10.0
  end subroutine clearPleC
  
  subroutine checkPoint3()
    implicit none
    integer :: i

    call initZero(PleC_st) 
    call initZero(DivL_st) 
    call initZero(CckA_st) 
    
    do i=1, MM/10
      pp(i,PleC_st) = 1.0
      pp(i,DivL_st) = 1.0
      pp(i,CckA_st) = 1.0
    end do
    p(PleC_st) = MM/10.0
    p(DivL_st) = MM/10.0
      
    do i=MM*9/10+1, MM
      pp(i,CckA_st) = 1.0
    end do 
    p(CckA_st) = 2*MM/10.0
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
            
        case (deg : decomposition)
            propensity= rt * p(react1_ita)

        case (sticky)
            totPro(ita) = totPro(ita) - prePro(bin_ita, ita)
            prePro(bin_ita, ita) = pp(bin_ita, react1_ita)*pp(bin_ita, aux_ita)
            totPro(ita) = totPro(ita) + prePro(bin_ita, ita)
            if (bin_diff.NE.0) then
                totPro(ita) = totPro(ita) - prePro(bin_diff, ita)
                prePro(bin_diff, ita) = pp(bin_diff, react1_ita)*pp(bin_diff, aux_ita)
                totPro(ita) = totPro(ita) + prePro(bin_diff, ita)
            end if
            propensity= rt*totPro(ita)

!        case (catalytic)
!            totPro(ita) = totPro(ita) - prePro(bin_ita, ita)
!            prePro(bin_ita, ita) = pp(bin_ita, react1_ita)*pp(bin_ita, aux_ita)
!            totPro(ita) = totPro(ita) + prePro(bin_ita, ita)
!            if (bin_diff.NE.0) then
!                totPro(ita) = totPro(ita) - prePro(bin_diff, ita)
!                prePro(bin_diff, ita) = pp(bin_diff, react1_ita)*pp(bin_diff, aux_ita)
!                totPro(ita) = totPro(ita) + prePro(bin_diff, ita)
!            end if
!            propensity= rt*totPro(ita)/h

        case (displacement : combination)
            totPro(ita) = totPro(ita) - prePro(bin_ita, ita)
            prePro(bin_ita, ita) = pp(bin_ita, react1_ita)*pp(bin_ita, react2_ita)
            totPro(ita) = totPro(ita) + prePro(bin_ita, ita)
            if (bin_diff.NE.0) then
                totPro(ita) = totPro(ita) - prePro(bin_diff, ita)
                prePro(bin_diff, ita) = pp(bin_diff, react1_ita)*pp(bin_diff, react2_ita)
                totPro(ita) = totPro(ita) + prePro(bin_diff, ita)
            end if
            propensity= rt*totPro(ita)/h

        case (hill)
            do i=0, HMEAN-1
                avg = 0.0
                do j=1, MM/HMEAN
                    avg = avg + pp(i*MM/HMEAN+j,aux_ita)
                end do
                avg = avg/(MM/HMEAN)
                avg = avg*avg*avg*avg
                do j=1, MM/HMEAN
                    auxt(i*MM/HMEAN+j) = avg
                end do
            end do

            Kmdl = 0.5*h*SCALE
            Kmdl = Kmdl*Kmdl*Kmdl*Kmdl
            do j=1,MM
                prePro(j,ita) = pp(j,react1_ita)*(auxt(j)/(auxt(j)+Kmdl))
                temp = temp + prePro(j,ita)
            end do
            totPro(ita) = temp
            propensity = rt*temp

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
        call random_number(random)
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
        if (prod1_ita .LE. SNUM_ODE) then
          IJth(y, prod1_ita, bin_ita) = IJth(y, prod1_ita, bin_ita) + ADD
        end if

      case (deg) ! A -> null
        call selectBin1st(bin_ita, react1_ita)
        p(react1_ita) = p(react1_ita) + REMOVE
        pp(bin_ita, react1_ita) = pp(bin_ita, react1_ita) + REMOVE
        if (react1_ita .LE. SNUM_ODE) then
          IJth(y, react1_ita, bin_ita) = IJth(y, react1_ita, bin_ita) + REMOVE
        end if
            
      case (first_order) ! A -> B
        call selectBin1st(bin_ita, react1_ita)
        p(react1_ita) = p(react1_ita) + REMOVE
        pp(bin_ita, react1_ita) = pp(bin_ita, react1_ita) + REMOVE
        p(prod1_ita) = p(prod1_ita) + ADD
        pp(bin_ita, prod1_ita) = pp(bin_ita, prod1_ita) + ADD
        if (react1_ita .LE. SNUM_ODE) then
          IJth(y, react1_ita, bin_ita) = IJth(y, react1_ita, bin_ita) + REMOVE
        end if
        if (prod1_ita .LE. SNUM_ODE) then
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
        if (react1_ita .LE. SNUM_ODE) then
          IJth(y, react1_ita, bin_ita) = IJth(y, react1_ita, bin_ita) + REMOVE
        end if
        if (prod1_ita .LE. SNUM_ODE) then
          IJth(y, prod1_ita, bin_ita) = IJth(y, prod1_ita, bin_ita) + ADD
        end if
        if (prod2_ita .LE. SNUM_ODE) then
          IJth(y, prod2_ita, bin_ita) = IJth(y, prod2_ita, bin_ita) + ADD
        end if
  
      case (sticky) ! A --> C aux
        call selectBin2nd(bin_ita, ita)
        p(react1_ita) = p(react1_ita) + REMOVE
        p(prod1_ita) = p(prod1_ita) + ADD
        pp(bin_ita, react1_ita) = pp(bin_ita, react1_ita) + REMOVE
        pp(bin_ita, prod1_ita) = pp(bin_ita, prod1_ita) + ADD
        if (react1_ita .LE. SNUM_ODE) then
          IJth(y, react1_ita, bin_ita) = IJth(y, react1_ita, bin_ita) + REMOVE
        end if
        if (prod1_ita .LE. SNUM_ODE) then
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
        if (react1_ita .LE. SNUM_ODE) then
          IJth(y, react1_ita, bin_ita) = IJth(y, react1_ita, bin_ita) + REMOVE
        end if
        if (react2_ita .LE. SNUM_ODE) then
          IJth(y, react2_ita, bin_ita) = IJth(y, react2_ita, bin_ita) + REMOVE
        end if
        if (prod1_ita .LE. SNUM_ODE) then
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
        if (react1_ita .LE. SNUM_ODE) then
          IJth(y, react1_ita, bin_ita) = IJth(y, react1_ita, bin_ita) + REMOVE
        end if
        if (react2_ita .LE. SNUM_ODE) then
          IJth(y, react2_ita, bin_ita) = IJth(y, react2_ita, bin_ita) + REMOVE
        end if
        if (prod1_ita .LE. SNUM_ODE) then
          IJth(y, prod1_ita, bin_ita) = IJth(y, prod1_ita, bin_ita) + ADD
        end if
        if (prod2_ita .LE. SNUM_ODE) then
          IJth(y, prod2_ita, bin_ita) = IJth(y, prod2_ita, bin_ita) + ADD
        end if

    case (hill)
        call selectBin2nd(bin_ita, ita)
        p(react1_ita) = p(react1_ita) + REMOVE
        p(prod1_ita) = p(prod1_ita) + ADD
        pp(bin_ita, react1_ita) = pp(bin_ita, react1_ita) + REMOVE
        pp(bin_ita, prod1_ita) = pp(bin_ita, prod1_ita) + ADD

      case default
        print *, "Reaction type in population_update: ", react_type, " not found !"
        stop -1
    
    end select
!if ( (react1_ita .ne. 0) ) then
!if (p(react1_ita) .LT. 0.0) then
!print *, 'react1_ita',react1_ita, 'p' ,p(react1_ita)
!print *, 'pp', pp
!end if
!end if
!if ( (react2_ita .ne. 0) ) then
!if (p(react2_ita) .LT. 0.0) then
!print *, 'react2_ita',react2_ita, 'p' ,p(react2_ita)
!end if
!end if



!check PleC_f
!random = 0.0
!i=5
!do j=1,MM
!random = random + pp(j,i)
!end do
!if ( ABS(random-p(i)) .GT. 0.01 ) then
!print *,'p(i)',p(i)
!do j=1,MM
!print *, pp(j,i)
!end do
!end if


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
    
      case (sticky)
        do j=1, MM
          prePro(j,ita) = pp(j,react1_ita)*pp(j,aux_ita)
          temp = temp + prePro(j,ita)
        end do
        totPro(ita) = temp
        propensity = rt*temp
    
!      case (catalytic)
!        do j=1, MM
!          prePro(j,ita) = pp(j,react1_ita)*pp(j,aux_ita)
!          temp = temp + prePro(j,ita)
!        end do
!        totPro(ita) = temp
!        propensity = rt*temp/h
      
      case (displacement : combination)
        do j=1, MM
          prePro(j,ita) = pp(j,react1_ita)*pp(j,react2_ita)
          temp = temp + prePro(j,ita)
        end do
        totPro(ita) = temp
        propensity = rt*temp/h
      
      case (hill)
        do i=0, HMEAN-1
            avg = 0.0
            do j=1, MM/HMEAN
                avg = avg + pp(i*MM/HMEAN+j,aux_ita)
            end do
            avg = avg/(MM/HMEAN)
            avg = avg*avg*avg*avg
            do j=1, MM/HMEAN
                auxt(i*MM/HMEAN+j) = avg
            end do
        end do

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

      case (syn)
        propensity= rt * p(aux_ita)

      case (deg : decomposition)
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

    call random_number(random)
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

