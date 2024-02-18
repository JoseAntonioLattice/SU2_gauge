module local_update_algorithms

  use data_types_observables
  use matrix_operations
  use iso_fortran_env, only : dp => real64, i4 => int32

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

  subroutine heatbath_gattringer(U,x,y,mu,beta)

    type(link_variable), dimension(:,:), intent(inout) :: U
    integer(i4), intent(in) :: x, y, mu
    real(dp), intent(in) :: beta
    type(complex_2x2_matrix) :: A, V, XX
    complex(dp) :: c_a, b
    real(dp) :: det_A, sqrt_det_A
    real(dp) :: x0, x_vec(3), norm_x
    real(dp) :: r(3), s
    real(dp) :: lambdasq

    A = staples(U,x,y,mu)

    det_A = det(A)
    if (det_A <= 0.0_dp)then
       call create_unbiased_update(U(x,y)%link(mu))
       return
    end if

    sqrt_det_A = sqrt(det_A)
    V%matrix = A%matrix/det_A

    call generate_lambdasq(sqrt_det_A,beta,lambdasq,s)
    do while ( s**2 > 1.0_dp - lambdasq)
       !call generate_lambdasq(det_A,beta,lambdasq,s)
       CALL random_number(s)
    end do

    x0 = 1.0_dp - 2*lambdasq
    norm_x = sqrt(1.0_dp - x0**2)

    call random_number(r)
    r = 2*r - 1.0_dp
    !if( (norm2(r))**2 > 1.0_dp) return
    do while ( (norm2(r))**2 > 1.0_dp )
        call random_number(r)
        r = 2*r - 1.0_dp
    end do
    x_vec = norm_x * r/norm2(r)

    c_a = cmplx(x0,x_vec(1),dp)
    b = cmplx(x_vec(2),x_vec(3),dp)

    XX%matrix(1,1) = c_a
    XX%matrix(1,2) = b
    XX%matrix(2,1) = -conjg(b)
    XX%matrix(2,2) =  conjg(c_a)

    U(x,y)%link(mu) = XX * dagger(V)

  end subroutine heatbath_gattringer


  subroutine heatbath(U,x,y,mu,beta)
    type(link_variable), dimension(:,:), intent(inout) :: U
    integer(i4), intent(in) :: x, y, mu
    real(dp), intent(in) :: beta

    real(dp) :: det_staple,sqrt_det_staple, a, b, a0, a1, a2, a3, norm_a, theta, phi
    type(complex_2x2_matrix) :: staple, XX, V
    logical :: boolean
    real(dp), parameter :: pi = acos(-1.0_dp)
    complex(dp) :: c1,c2


    staple = staples(U,x,y,mu)
    det_staple = det(staple)
    if (det_staple <= 0.0_dp)then
       call create_unbiased_update(U(x,y)%link(mu))
       return
    end if
    sqrt_det_staple = sqrt(det_staple)

    V%matrix = staple%matrix/sqrt_det_staple

    a = exp(-2*beta*sqrt_det_staple)
    b = 1.0_dp
    print*, det_staple, sqrt_det_staple, a
    boolean = .false.
    do while( boolean .eqv. .false.)
      r = random_uniform(a,b)
      a0 = 1.0_dp + log(r)/(beta * sqrt_det_staple)
      call random_number(r)
      if( r > 1.0_dp - sqrt(1.0_dp - a0**2)) boolean = .true.
    end do

    norm_a = sqrt(1.0_dp - a0**2)
    theta = random_uniform(0.0_dp, pi)
    phi = random_uniform(0.0_dp, 2*pi)
    a1 = norm_a * cos(phi) * sin(theta)
    a2 = norm_a * sin(phi) * sin(theta)
    a3 = norm_a * cos(theta)

    c1 = cmplx(a0,a1,dp)
    c2 = cmplx(a2,a3,dp)

    XX%matrix(1,1) = c1
    XX%matrix(1,2) = c2
    XX%matrix(2,1) = -conjg(c2)
    XX%matrix(2,2) =  conjg(c1)

    U(x,y)%link(mu) = XX * dagger(V)

  end subroutine heatbath


!  subroutine heatbath2(U,x,y,mu,beta)
!    type(link_variable), dimension(:,:), intent(inout) :: U
!    integer(i4), intent(in) :: x, y, mu
!    real(dp), intent(in) :: beta
!    type(complex_2x2_matrix) :: V, sigma, hat_sigma
!    real(dp) :: sqrt_det_sigma, rho, y0, Y, W3
!
!    sigma = staples(U,x,y,mu)
!    sqrt_det_sigma = sqrt(det(sigma))
!    hat_sigma = sigma/sqrt_det_sigma
!
!    rho = beta * sqrt_det_sigma
!
!    !y0 = rho * ( 1.0_dp - a(0) )
!
!    call random_number(W3)
!    !do while( 2* W3**2 > 2.0_dp - Y)
!
!    !end do
!
!  end subroutine heatbath2

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

    Up%matrix(1,1) = a
    Up%matrix(1,2) = b
    Up%matrix(2,1) = -conjg(b)
    Up%matrix(2,2) =  conjg(a)

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

    Up%matrix(1,1) = a
    Up%matrix(1,2) = b
    Up%matrix(2,1) = -conjg(b)
    Up%matrix(2,2) =  conjg(a)

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

  pure function staples(U,x,y,mu) result(A)
    use periodic_boundary_contidions_mod, only : ip, im
    use parameters, only : d

    type(link_variable), dimension(:,:), intent(in) :: U
    integer(i4), intent(in) :: x, y, mu
    integer(i4) :: nu
    type(complex_2x2_matrix) :: A

    !A%matrix = 0.0_dp
    if (mu == 1)then
      A =        U(x,   y )%link(2)  * U(x,ip(y))%link(mu) * dagger(U(ip(x),   y )%link(2))  +  &
          dagger(U(x,im(y))%link(2)) * U(x,im(y))%link(mu) *        U(ip(x),im(y))%link(2)! + A
    else if(mu == 2)then
      A =        U(x,    y)%link(1)  * U(ip(x),y)%link(mu) * dagger(U(   x, ip(y))%link(1))  +  &
          dagger(U(im(x),y)%link(1)) * U(im(x),y)%link(mu) *        U(im(x),ip(y))%link(1)! + A1
    end if
  end function staples

  function DS(U,mu,Up,x,y)
    type(link_variable), dimension(:,:), intent(in) :: U
    integer(i4), intent(in) :: x, y, mu
    type(complex_2x2_matrix), intent(out) :: Up
    real(dp) :: DS

    call create_unbiased_update(Up)
    !call create_update(Up); Up = Up * U(x,y)%link(mu)

    DS = - real( tr( (Up - U(x,y)%link(mu)) * dagger(staples(U,x,y,mu)) ),dp )

  end function DS


  function random_uniform(a,b) result(y)
    real(dp), intent(in) :: a, b
    real(dp) :: y, r

    call random_number(r)

    y = a + r * ( b - a )

  end function random_uniform


end module local_update_algorithms
