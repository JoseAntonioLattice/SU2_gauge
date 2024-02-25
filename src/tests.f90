program main
  use get_index_mod
  implicit none

  integer :: x = 4, d = 3, L = 2

  print*, "Array of ",x, get_index_array(x,d,L)

end program main
