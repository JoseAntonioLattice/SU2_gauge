program main

  use iso_fortran_env, only : dp => real64
  use periodic_boundary_contidions_mod
  use parameters
  use arrays
  use starts
  use dynamics
  use data_types_observables
  implicit none

  type(complex_2x2_matrix) :: Up, Uold 
  integer :: i
  real(dp) :: S
  
  !Read input parameters
  call read_input_parameters()

  !Allocate variables
  allocate(U(L,L))
  call set_periodic_bounds(L)


  open(unit = 100, file = 'action_glauber.dat')
  
  !Initialie variables
  call cold_start(U)
  call take_measurements(U,S)
  write(100,*) S

  !Thermalization
  do i = 1, 1000
     call sweeps(U,L)
     call take_measurements(U,S)
     write(100,*) S
  end do

  
  ! Visualize
  !print*, U(L,L,L,L)%link(4)

  !print*, "Create update"
  !call create_update(Up)

  !Uold = U(L,L)%link(2)
  !U(L,L)%link(2) = Up
  
  !print*, "Action of a SU(2) 2d Configuration before gauge transformation:"
  !call sweeps(U,L)
  !call take_measurements(U)

  !U(L,L)%link(2) = Uold

  !print*, DS(U,2,Up,L,L)
  
  !print*, "Action of a SU(2) 2d Configuration after gauge transformation: "
  !call gauge_transformation(U)
  !call take_measurements(U) 

end program main
