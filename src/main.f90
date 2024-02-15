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
  integer :: i, itemp
  real(dp) :: S, suma_S
  real(dp), allocatable, dimension(:) :: beta
  
  !Read input parameters
  call read_input_parameters()

  !Allocate variables
  allocate(U(L,L))
  call set_periodic_bounds(L)

  beta = [(i*0.1, i = 0, 100)]


  open(unit = 100, file = 'action.dat')

  do itemp = 1, size(beta)
     suma_S = 0.0_dp
     !beta = 1/temperature(itemp)
     !Initialie variables
     call cold_start(U)
     
     !Thermalization
     do i = 1, 1000
        call sweeps(U,L,beta(itemp),N)
     end do

     do i = 1, 1000
        call sweeps(U,L,beta(itemp),2)
        if( mod(i,10) == 0)then
           call take_measurements(U,L,beta(itemp),N,2,S)
           suma_s = suma_s + S
        end if
     end do
     write(100,*) beta(itemp), suma_S/100
  end do
  
  
end program main
