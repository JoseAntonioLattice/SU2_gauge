module starts

  use iso_fortran_env, only : dp => real64, i4 => int32

  implicit none

  private
  public :: cold_start, hot_start

contains


  subroutine cold_start(U)
    use data_types_observables, only : link_variable
    type(link_variable), intent(out), dimension(:,:) :: U
    integer(i4) :: L,x,y,mu!,z,w

    L = size(U(:,1))

    do x = 1, L
      do y = 1, L
        !do z = 1, L
        !  do w = 1, L
            do mu = 1, 2
               U(x,y)%link(mu)%matrix = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp], [2,2])
            end do
        !  end do
        !end do
      end do
    end do

  end subroutine cold_start

  subroutine hot_start(U)
    use data_types_observables, only : link_variable
    type(link_variable), intent(out), dimension(:,:) :: U
    integer(i4) :: L,x,y,mu!,z,w
    real(dp) :: r(4)
    complex(dp) :: a, b

    L = size(U(:,1))

    do x = 1, L
      do y = 1, L
        !do z = 1, L
        !  do w = 1, L
            do mu = 1, 2
              call  random_number(r)
              r = r - 0.5_dp
              r = r/norm2(r)
              a = cmplx(r(1),r(2))
              b = cmplx(r(3),r(4))
              U(x,y)%link(mu)%matrix = reshape([a,-conjg(b),b,conjg(a)], [2,2])
            end do
        !  end do
        !end do
      end do
    end do

  end subroutine hot_start

end module starts
