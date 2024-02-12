module arrays

  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none

  private
  
  type complex_2x2_matrix
     complex(dp), dimension(2,2) :: matrix
  end type complex_2x2_matrix

  type link_variable
     type(complex_2x2_matrix), dimension(4) :: link
  end type link_variable

  type(link_variable), allocatable, dimension(:,:,:,:) :: U

  public :: U, link_variable
  
end module arrays
