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

  !beta = [(i*0.1, i = 1, 100)]
  beta = (/1,2,3,4,5/)

  open(unit = 100, file = 'action_heatbath.dat')

  do itemp = 1, size(beta)
     suma_S = 0.0_dp
     !beta = 1/temperature(itemp)
     !Initialie variables
     call cold_start(U)
     
     !Thermalization
     print*,"Initializing thermalization at beta = ", beta(itemp)
     do i = 1, 100
        print*, i
        call sweeps(U,L,beta(itemp),N)
     end do
     print*, "Thermalization done at beta = ", beta(itemp)
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
