program main
  use data_types_observables
  use matrix_operations
  implicit none

  type(link_variable) :: u1(2,2), u2, u3, u4
  type(complex_2x2_matrix) :: a, b, c,d

  a%matrix = reshape([(1.0_dp,1.0_dp),(2.0_dp,1.0_dp),(3.0_dp,1.0_dp),(4.0_dp,1.0_dp)],[2,2])
  b%matrix = reshape([(5.0_dp,1.0_dp),(6.0_dp,1.0_dp),(7.0_dp,1.0_dp),(8.0_dp,1.0_dp)],[2,2])
  c%matrix = reshape([(9.0_dp,1.0_dp),(10.0_dp,1.0_dp),(11.0_dp,1.0_dp),(12.0_dp,1.0_dp)],[2,2])

  d = a*b*c

  u1(1,1)%link(1) = a
  u1(1,2)%link(2) = b


  !print*, matmul(a%matrix,b%matrix)
  !print*, a*b*c
  print*,"complex_2x2_matrix = ", a
  print*, "link_variable     = ", u1%link(1)
  print*,tr(u1(1,1)%link(1))
  print*, dagger(u1(1,1)%link(1))
  print*, u1(1,1)%link(1) * u1(1,2)%link(2)


end program main
