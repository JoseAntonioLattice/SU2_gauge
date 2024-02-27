module data_types_observables


  use iso_fortran_env, only : dp => real64, i4 => int32

  implicit none

  type complex_2x2_matrix
     complex(dp), dimension(2,2) :: matrix
  end type complex_2x2_matrix

  type link_variable
     type(complex_2x2_matrix), dimension(6) :: link
  end type link_variable

  type observable
     real(dp), allocatable :: array(:)
     real(dp) :: avr
     real(dp) :: err
  end type observable


  type(observable) :: E_p
  real(dp), allocatable, dimension(:) :: beta


  private
  public :: link_variable, complex_2x2_matrix, beta, E_p


end module data_types_observables
