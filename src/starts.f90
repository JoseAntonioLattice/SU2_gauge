module starts

  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters, only : d
  implicit none

  private
  public :: cold_start, hot_start

contains


  subroutine cold_start(U)
    use data_types_observables, only : link_variable
    type(link_variable), intent(out), dimension(:) :: U
    integer(i4) :: L,x,mu!,z,w

    L = size(U)

    do x = 1, L
       do mu = 1, 3
          U(x)%link(mu)%matrix = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp], [2,2])
       end do
    end do

  end subroutine cold_start

  subroutine hot_start(U)
    use data_types_observables, only : link_variable
    type(link_variable), intent(out), dimension(:) :: U
    integer(i4) :: V,x,mu!,z,w
    real(dp) :: r(4)
    complex(dp) :: a, b

    V = size(U)

    do x = 1, V
       do mu = 1, 3
          call  random_number(r)
          r = r - 0.5_dp
          r = r/norm2(r)
          a = cmplx(r(1),r(2),dp)
          b = cmplx(r(3),r(4),dp)
          U(x)%link(mu)%matrix = reshape([a,-conjg(b),b,conjg(a)], [2,2])
       end do
    end do

  end subroutine hot_start

end module starts
