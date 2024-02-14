program main

  use periodic_boundary_contidions_mod
  use parameters
  use arrays
  use starts
  use dynamics
  implicit none

  !Read input parameters
  call read_input_parameters()

  !Allocate variables
  allocate(U(L,L))
  call set_periodic_bounds(L)

  !Initialie variables
  call cold_start(U)

  ! Visualize
  !print*, U(L,L,L,L)%link(4)

  print*, "Create update"
  !call create_update(U(L,L,L,L)%link(4))
  !call sweeps(U,L)
  !print*, U(L,L)%link(2)

  !print*, "Trace of 1_2x2 = ", tr(U(L,L)%link(2)%matrix)

  print*, "Action of a SU(2) 2d Configuration before gauge transformation:"
  call sweeps(U,L)
  call take_measurements(U)

  print*, "Action of a SU(2) 2d Configuration after gauge transformation: "
  call gauge_transformation(U)
  call take_measurements(U) 

end program main
