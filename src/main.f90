program main

  use iso_fortran_env, only : dp => real64
  use periodic_boundary_contidions_mod
  use parameters
  use arrays
  use starts
  use dynamics
  use data_types_observables
  use statistics
  implicit none

  integer :: i, itemp
  real(dp), allocatable, dimension(:) :: beta
  type observable
    real(dp), allocatable :: array(:)
    real(dp) :: avr
    real(dp) :: err
  end type observable
  type(observable):: S

  !Read input parameters
  call read_input_parameters()

  !Allocate variables
  allocate(U(L,L))
  allocate(S%array(N_measurements))
  call set_periodic_bounds(L)

  beta = [(i*1.0_dp, i = 1, 100)]
  !beta = 1/temperature

  open(unit = 100, file = 'data/action_'//trim(algorithm)//'.dat')

  !Initialie variables
  !call take_measurements(U,L,N,d,S%array(1))
  !print*, S%array(1)
  do itemp = 1, size(beta)
  call hot_start(U)
     !Thermalization
      !print*, "Thermalization at beta = ", beta(itemp)
     do i = 1, N_thermalization
        !print*, "sweep ",i
        call sweeps(U,L,beta(itemp),N,d,algorithm)
     end do

     do i = 1, N_measurements*N_skip
      ! print*, i
        call sweeps(U,L,beta(itemp),N,d,algorithm)
        if( mod(i,N_skip) == 0)then
           !call take_measurements(U,L,N,d,S%array(i/N_skip))
          S%array(i/N_skip) = action(U,-1.0_dp/N,d)/(L**d)
        end if
     end do

      call std_err(S%array,S%avr,S%err)
     write(100,*) beta(itemp), S%avr, S%err
  end do


end program main
