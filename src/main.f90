program main

  use parameters
  use arrays
  use starts
  use dynamics
  implicit none

  !Read input parameters
  call read_input_parameters()

  !Allocate variables
  allocate(U(L,L,L,L))

  !Initialie variables
  call cold_start(U)

  ! Visualize
  print*, U(L,L,L,L)%link(4)

  print*, "Create update"
  !call create_update(U(L,L,L,L)%link(4))
  call sweeps(U,L)
  print*, U(L,L,L,L)%link(4)


end program main
