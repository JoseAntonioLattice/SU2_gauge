program main

  use parameters
  use arrays
  use starts
  implicit none

  !Allocate variables
  allocate(U(L,L,L,L))


  call cold_start(U)
  
  ! Visualize
  print*, U(L,L,L,L)%link(4)
  
  
  
end program main
