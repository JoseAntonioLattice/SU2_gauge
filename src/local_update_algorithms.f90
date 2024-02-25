module local_update_algorithms

  use data_types_observables
  use matrix_operations
  use iso_fortran_env, only : dp => real64, i4 => int32
  use periodic_boundary_conditions_mod, only : ip_func, im_func
  use get_index_mod

  implicit none

contains

  subroutine metropolis(Delta_S,U,Up)
    real(dp), intent(in) :: Delta_S
    complex(dp), dimension(:,:), intent(inout) :: U
    complex(dp), dimension(:,:), intent(in) :: Up
    real(dp) :: r, prob

    prob = min(1.0_dp,exp(-Delta_S))

    call random_number(r)
    if( prob >= r )then
       U = Up
    end if

  end subroutine metropolis


 subroutine glauber(Delta_S,U,Up)
    real(dp), intent(in) :: Delta_S
    complex(dp), dimension(:,:), intent(inout) :: U
    complex(dp), dimension(:,:), intent(in) :: Up
    real(dp) :: r

    call random_number(r)
    if ( 1/(exp(Delta_S) + 1.0_dp) > r)then
       U = Up
    end if

  end subroutine glauber

  subroutine heatbath_gattringer(U,x,mu,beta)

    type(link_variable), dimension(:), intent(inout) :: U
    integer(i4), intent(in) :: x, mu
    real(dp), intent(in) :: beta
    type(complex_2x2_matrix) :: A, V, XX
    complex(dp) :: c_a, b
    real(dp) :: det_A, sqrt_det_A
    real(dp) :: x0, x_vec(3), norm_x
    real(dp) :: r(3), s
    real(dp) :: lambdasq
    real(dp), parameter :: pi = acos(-1.0_dp)
    complex(dp) :: c1,c2, det_staple_complex
    real(dp) :: d1, d2, d3, d4

    A = staples(U,x,mu)

    det_A = det(A)
    if (det_A <= 0.0_dp)then
       call create_unbiased_update(U(x)%link(mu))
       return
    end if

    sqrt_det_A = sqrt(det_A)
    V%matrix = A%matrix/sqrt_det_A

    call generate_lambdasq(sqrt_det_A,beta,lambdasq,s)
    do while ( s**2 > 1.0_dp - lambdasq )
       call generate_lambdasq(det_A,beta,lambdasq,s)
    end do

    x0 = 1.0_dp - 2 * lambdasq
    norm_x = sqrt(1.0_dp - x0*x0)

    x_vec = norm_x * random_vector()

    c1 = cmplx(x0,x_vec(1),dp)
    c2 = cmplx(x_vec(2),x_vec(3),dp)

    XX = SU2_matrix(c1,c2)

    U(x)%link(mu) = XX * V
    !print*, det(XX), det(V), det(XX * V)

  end subroutine heatbath_gattringer


  subroutine heatbath(U,x,mu,beta)
    type(link_variable), dimension(:), intent(inout) :: U
    integer(i4), intent(in) :: x, mu
    real(dp), intent(in) :: beta

    real(dp) :: a, b, a0, a_vec(3), norm_a
    type(complex_2x2_matrix) :: staple, XX, V, prod
    logical :: boolean

    complex(dp) :: c1,c2
    real(dp) :: det_staple, sqrt_det_staple, d1, d2, d3, d4,r

    staple = staples(U,x,mu)
    det_staple = det(staple)
    sqrt_det_staple = sqrt(det_staple)
    V%matrix = staple%matrix/sqrt_det_staple

    a = exp(-2*beta*sqrt_det_staple)
    b = 1.0_dp

    boolean = .false.
    do while( boolean .eqv. .false.)
      r = random_uniform(a,b)
      a0 = 1.0_dp + log(r)/(beta * sqrt_det_staple)
      call random_number(r)
      if( r > 1.0_dp - sqrt(1.0_dp - a0**2)) boolean = .true.
    end do

    norm_a = sqrt(1.0_dp - a0**2)

    a_vec = norm_a * random_vector()

    c1 = cmplx(a0,a_vec(1),dp)
    c2 = cmplx(a_vec(2),a_vec(3),dp)

    XX = SU2_matrix(c1,c2)

    U(x)%link(mu) = XX * V

  end subroutine heatbath

  function random_vector() result(y)
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: theta, phi
    real(dp), dimension(3) :: y

    theta = random_uniform(0.0_dp, pi)
    phi   = random_uniform(0.0_dp, 2*pi)

    y(1) = cos(phi) * sin(theta)
    y(2) = sin(phi) * sin(theta)
    y(3) = cos(theta)

  end function random_vector

  function SU2_matrix(a,b) result(matrix)
    complex(dp), intent(in)  :: a, b
    type(complex_2x2_matrix) :: matrix
      matrix%matrix(1,1) = a
      matrix%matrix(1,2) = b
      matrix%matrix(2,1) = -conjg(b)
      matrix%matrix(2,2) =  conjg(a)
  end function SU2_matrix

  subroutine generate_lambdasq(det_A,beta,lambdasq,s)
    real(dp), intent(in) :: det_A,beta
    real(dp), intent(out):: lambdasq, s
    real(dp) :: r(3)
    real(dp), parameter :: pi = acos(-1.0_dp)

    call random_number(r)
    r = 1.0_dp - r
    lambdasq = -(1/(2*det_A*beta)) * ( log(r(1)) + (cos(2*pi*r(2)))**2 * log(r(3)) )
    call random_number(s)
  end subroutine generate_lambdasq


  subroutine create_update(Up)
    type(complex_2x2_matrix), intent(out) :: Up
    complex(dp) :: a, b
    real(dp), dimension(0:3) :: r, x
    real(dp), parameter :: eps = 0.1_dp
    real(dp) :: norm_r

    call random_number(r)
    r = r - 0.5_dp
    norm_r = sqrt(r(1)**2 + r(2)**2 + r(3)**2)

    x(1:3) = eps*r(1:3)/norm_r
    x(0) = sgn(r(0)) * sqrt(1.0_dp - eps**2)

    a = cmplx(x(0),x(1),dp)
    b = cmplx(x(2),x(3),dp)

    Up = SU2_matrix(a,b)

  end subroutine create_update

  subroutine create_unbiased_update(Up)
    type(complex_2x2_matrix), intent(out) :: Up
    complex(dp) :: a, b
    real(dp), dimension(0:3) :: r

    call random_number(r)
    r = r - 0.5_dp
    r = r/norm2(r)
    a = cmplx(r(0),r(1),dp)
    b = cmplx(r(2),r(3),dp)

    Up = SU2_matrix(a,b)

  end subroutine create_unbiased_update

  pure function sgn(x)
    real(dp), intent(in) :: x
    integer(i4) :: sgn

    if( x > 0.0_dp )then
      sgn = 1
    elseif( x < 0.0_dp)then
      sgn = -1
    else
      sgn = 0
    end if

  end function sgn

  function staples(U,x,mu) result(A)

    use parameters, only : L, d

    type(link_variable), dimension(:), intent(in) :: U
    integer(i4), intent(in) :: x, mu
    integer(i4) :: nu
    type(complex_2x2_matrix) :: A

    integer(i4) :: ipx_mu, ipx_nu, imx_nu, ipx_mu_imx_nu

    A%matrix = 0.0_dp
    do nu = 1, d
       if(nu .ne. mu)then
          !print*, "inside staples", "mu=",mu,"nu=",nu,x, get_index_array(x,d,L), &
          !     ip_func(get_index_array(x,d,L),mu), ip_func(get_index_array(x,d,L),nu), im_func(get_index_array(x,d,L),nu)
          ipx_mu = get_index(ip_func(get_index_array(x,d,L),mu),d,L)
          ipx_nu = get_index(ip_func(get_index_array(x,d,L),nu),d,L)
          imx_nu = get_index(im_func(get_index_array(x,d,L),nu),d,L)
          !print*, 'imx_nu',imx_nu, get_index_array(imx_nu,d,L)
   ipx_mu_imx_nu = get_index(ip_func(get_index_array(imx_nu,d,L),mu),d,L)
          !print*, "before computing staples"

          A = A +    U(   x  )%link(nu)  * U(ipx_nu)%link(mu) * dagger(U(ipx_mu       )%link(nu)) +  &
              dagger(U(imx_nu)%link(nu)) * U(imx_nu)%link(mu) *        U(ipx_mu_imx_nu)%link(nu)
       end if
    end do
    end function staples

  function DS(U,mu,Up,x)
    type(link_variable), dimension(:), intent(in) :: U
    integer(i4), intent(in) :: x, mu
    type(complex_2x2_matrix), intent(out) :: Up
    real(dp) :: DS

    call create_unbiased_update(Up)
    !call create_update(Up); Up = Up * U(x,y)%link(mu)

    DS = - real( tr( (Up - U(x)%link(mu)) * dagger(staples(U,x,mu)) ) )

  end function DS


  function random_uniform(a,b) result(y)
    real(dp), intent(in) :: a, b
    real(dp) :: y, r

    call random_number(r)

    y = a + r * ( b - a )

  end function random_uniform

end module local_update_algorithms
