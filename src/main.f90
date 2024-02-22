program main

  use iso_fortran_env, only : dp => real64
  use periodic_boundary_conditions_mod
  use parameters
  use arrays
  use starts
  use dynamics
  use data_types_observables
  use statistics
  implicit none

  integer :: i, itemp
 
  !Read input parameters
  call read_input_parameters()

  !Allocate variables
  allocate(U(L**d))
  allocate(E_p%array(N_measurements))
  call set_periodic_bounds(L)

  beta = [(i*1.00_dp, i = 1, 100)]
  !beta = 1/temperature

  open(unit = 100, file = 'data/action_'//trim(algorithm)//'.dat')

  do itemp = 1, size(beta)
     call hot_start(U)
     call thermalization(U,L,beta(itemp),N,d,algorithm,N_thermalization)
     call measurements_sweeps(U,L,beta(itemp),N,d,algorithm,N_measurements,N_skip, E_p%array)
     call std_err(E_p%array,E_p%avr,E_p%err)
     write(100,*) beta(itemp), E_p%avr, E_p%err
  end do


end program main
