module arrays

  use data_types_observables
  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none

  type(link_variable), allocatable, dimension(:) :: U
  real(dp), allocatable, dimension(:) :: temperature


  private
  public :: U, temperature

end module arrays
