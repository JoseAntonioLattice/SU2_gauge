module data_types_observables


  use iso_fortran_env, only : dp => real64, i4 => int32

  implicit none

  type complex_2x2_matrix
     complex(dp), dimension(2,2) :: matrix
  end type complex_2x2_matrix

  type link_variable
     type(complex_2x2_matrix), dimension(4) :: link
  end type link_variable


  private
  public :: link_variable, complex_2x2_matrix


end module data_types_observables
