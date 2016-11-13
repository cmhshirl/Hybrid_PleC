program main
  use simulation
  use mpi
  use rand

  implicit none

  integer ierr, rank, num_procs
  double precision :: start, finish

  REAL(KIND=R8)::ISEED,USEED
  
  !Initialize MPI.
  call MPI_Init ( ierr )

  call MPI_COMM_RANK (MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
  
  call cpu_time(start)

  !====Read Input====
  call set_evn(rank)

  !====Init====
  call init_random_seed(rank)

  CALL RANDOM_NUMBER(ISEED)
  !WRITE(*,*)'SEED',ISEED
  USEED = DUSTAR(INT(2.0_R8**24*ISEED))

  call init_solver_parameters()
  call init_para()
  call init_ode()
  call init_ssa()

CALL RANDOM_NUMBER(ISEED)
  WRITE(*,*)'SEED',ISEED

  !====Main Loop====
  call DNAorigination()
  call cellFixed()
  call clearFire()
  call simulate(150.0d0)

  call clearFire()
  call cellGrow()
  call simulate(30.0d0)

  call introDivJ()
  call simulate(20.0d0)

  call DNAreplication()
  call clearPleC()
  call simulate(40.0d0)

  call checkPoint3();
  call simulate(30.0d0)

  call cpu_time(finish)
  print '("Time = ", f12.6, "hours.")', (finish-start)/3600.0

  !====Deallocate Memory====
  call deallocate_parameters()

  call MPI_Finalize ( ierr )

end program main

