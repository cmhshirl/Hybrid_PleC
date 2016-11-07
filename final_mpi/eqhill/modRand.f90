module rand
contains
  subroutine init_random_seed(rank)
      implicit none
      
      integer, intent(in) :: rank
      integer :: i, n, clock
      integer, dimension(:), allocatable :: seed
      call random_seed(size = n)
      allocate(seed(n))

      call system_clock(count=clock)

      seed = abs( mod( ((clock + 37 * (/ (i - 1, i = 1, n) /))*181)*((rank-83)*359), 104729) )
      call random_seed(put = seed)
      print *, 'rank:',rank, seed
      deallocate(seed)
  end subroutine init_random_seed
end module rand
