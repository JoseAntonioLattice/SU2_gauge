module dynamics

  use iso_fortran_env, only: dp => real64, i4 => int32
  use data_types_observables, only : link_variable, complex_2x2_matrix

  implicit none

  private !:: dp, i4, link_variable
  public :: sweeps, create_update, sgn, drand

contains

  subroutine sweeps(U,L)
    type(link_variable), intent(inout), dimension(:,:,:,:) :: U
    integer(i4), intent(in)  :: L
    type(complex_2x2_matrix) :: Up
    integer(i4) :: x, y, z, w, mu

    do x = 1, L
       do y = 1, L
          do z = 1, L
             do w = 1, L
                do mu = 1, 4
                   call create_update(Up)
                   U(x,y,z,w)%link(mu) = Up
                   !call metropolis(DS,U(x,y,z,w)%link(mu),Up)
                end do
             end do
          end do
       end do
    end do
  end subroutine sweeps

  subroutine create_update(Up)
    type(complex_2x2_matrix), intent(out) :: Up
    real(dp) :: x0, x1, x2, x3
    complex(dp) :: a, b
    real(dp), parameter :: eps = 0.5
    real(dp), dimension(4) :: r
    real(dp) :: norm_r

    r(1) = drand(-0.5_dp,0.5_dp)
    r(2) = drand(-0.5_dp,0.5_dp)
    r(3) = drand(-0.5_dp,0.5_dp)
    r(4) = drand(-0.5_dp,0.5_dp)

    norm_r = norm2(r)

    x1 = eps*r(1)/norm_r
    x2 = eps*r(2)/norm_r
    x3 = eps*r(3)/norm_r
    x0 = sgn(r(4))*sqrt(1.0_dp - eps**2)

    a = cmplx(x0,x1,dp)
    b = cmplx(x2,x3,dp)

    Up%matrix(1,1) = a
    Up%matrix(1,2) = b
    Up%matrix(2,1) = -conjg(b)
    Up%matrix(2,2) =  conjg(a)

  end subroutine create_update


  function drand(a,b) result(y)
    real(dp), intent(in) :: a, b
    real(dp) :: y
    real(dp) :: r

    call random_number(r)
    y = a + r*(b-a)

  end function drand

  pure function sgn(x) result(y)
    real(dp), intent(in) :: x
    integer(i4) :: y

    if( x > 0.0_dp)then
       y = 1
    elseif(x < 0.0_dp)then
       y = -1
    else
       y = 0
    end if
  end function sgn

end module dynamics
