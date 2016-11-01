program main
  use simulation


  implicit none
  
  double precision :: start, finish
  call cpu_time(start)

  !====Read Input====
  call set_evn()

  !====Init====
  call init_random_seed()
  call init_solver_parameters()
  call init_para()
  call init_ode()
  call init_ssa()


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

end program main

