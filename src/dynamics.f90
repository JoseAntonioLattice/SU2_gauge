module dynamics

  use iso_fortran_env, only: dp => real64, i4 => int32
  use data_types_observables, only : link_variable, complex_2x2_matrix

  implicit none

  private !:: dp, i4, link_variable
  public :: sweeps, create_update, sgn, drand, take_measurements, dagger, tr, gauge_transformation, DS

contains

  subroutine sweeps(U,L,beta,N)
    type(link_variable), intent(inout), dimension(:,:) :: U
    integer(i4), intent(in)  :: L
    integer(i4), intent(in) :: N
    real(dp), intent(in) :: beta
    type(complex_2x2_matrix) :: Up
    integer(i4) :: x, y, mu
    real(dp) :: Delta_S

    do x = 1, L
       do y = 1, L
          do mu = 1, 2
             call create_update(Up)
             Delta_S = (beta/N) * DS(U,mu,Up,x,y)
             call metropolis(Delta_S,U(x,y)%link(mu)%matrix,Up%matrix)
             !call heatbath(U,x,y,mu,beta)
          end do
       end do
    end do
  end subroutine sweeps

  subroutine metropolis(Delta_S,U,Up)
    real(dp), intent(in) :: Delta_S
    complex(dp), dimension(:,:), intent(inout) :: U
    complex(dp), dimension(:,:), intent(in) :: Up
    real(dp) :: r

    call random_number(r)
    if (Delta_S <= 0.0_dp)then
       U = Up
    elseif( (Delta_S > 0.0_dp) .and. (exp(-Delta_S) > r) )then
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

  subroutine heatbath(U,x,y,mu,beta)

    type(link_variable), dimension(:,:), intent(inout) :: U
    integer(i4), intent(in) :: x, y, mu
    real(dp), intent(in) :: beta
    complex(dp), dimension(2,2) :: A, V, XX
    complex(dp) :: c_a, b
    real(dp) :: det_A
    real(dp) :: x0, x_vec(3), norm_x
    real(dp) :: r(3), s
    real(dp) :: lambdasq
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: i, j

    A = staples(U,x,y,mu)

    det_A = sqrt(det(A))
    !print*, det_A
    if (det_A <= 0.0_dp)then
       call create_update(U(x,y)%link(mu))
       return
    end if

    V = A/det_A

    call generate_lambdasq(det_A,beta,lambdasq,s)

    do while ( s**2 > 1.0_dp - lambdasq)
       call generate_lambdasq(det_A,beta,lambdasq,s)
    end do

    x0 = 1 - 2*lambdasq
    norm_x = sqrt(1.0_dp - x0**2)

    call random_number(r)
    r = 2*r - 1.0_dp
    do while ( sum(r**2) > 1.0_dp )
        call random_number(r)
        r = 2*r - 1.0_dp
    end do
    x_vec = norm_x * r/norm2(r)

    c_a = cmplx(x0,x_vec(1),dp)
    b = cmplx(x_vec(2),x_vec(3),dp)

    XX(1,1) = c_a
    XX(1,2) = b
    XX(2,1) = -conjg(b)
    XX(2,2) =  conjg(c_a)

    U(x,y)%link(mu)%matrix = matmul(XX,dagger(V))

  end subroutine heatbath

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
    real(dp), dimension(4) :: r
    real(dp) :: norm_r

    call random_number(r)
    r = r - 0.5_dp

    norm_r = norm2(r)

    r = r/norm_r

    a = cmplx(r(1),r(2),dp)
    b = cmplx(r(3),r(4),dp)

    Up%matrix(1,1) = a
    Up%matrix(1,2) = b
    Up%matrix(2,1) = -conjg(b)
    Up%matrix(2,2) =  conjg(a)

  end subroutine create_update

  pure function det(A)
    complex(dp), dimension(2,2), intent(in) :: A
    real(dp) :: det

    det = real(A(1,1)*A(2,2) - A(1,2)*A(2,1),dp)

  end function det


  subroutine take_measurements(U,L,beta,N,d,action)

    type(link_variable), dimension(:,:), intent(in) :: U
    integer(i4), intent(in) :: N,d,L
    real(dp), intent(in) :: beta
    integer(i4) :: x,y

    real(dp), intent(out) :: action
    action = 0.0_dp
    do x = 1, L
       do y = 1, L
          action = action + real(tr(plaquette(U,x,y,1,2)))
       end do
    end do

    action =  action/(N*L**d)

  end subroutine take_measurements

  function plaquette(U,x,y,mu,nu)
    use periodic_boundary_contidions_mod, only : ip

    type(link_variable), dimension(:,:), intent(in) :: U
    integer(i4), intent(in) :: x,y , mu, nu
    complex(dp), dimension(2,2) :: plaquette, prod1, prod2


    prod1 = matmul(U(x,y)%link(mu)%matrix,U(ip(x),y)%link(nu)%matrix)
    prod2 = matmul(dagger(U(x,ip(y))%link(mu)%matrix),dagger(U(x,y)%link(nu)%matrix))
    plaquette = matmul(prod1,prod2)
  end function plaquette

  pure function dagger(U) result(U_res)
    complex(dp), dimension(:,:), intent(in) :: U
    complex(dp), dimension(size(U(:,1)),size(U(1,:))) :: U_res

    U_res = conjg(U)
    U_res = transpose(U_res)

  end function dagger

  pure function tr(U) result(trace)
    complex(dp), dimension(:,:), intent(in) :: U
    complex(dp) :: trace
    integer :: i

    trace = sum([(U(i,i),i = 1, size(U(1,:)))])

  end function tr

  subroutine gauge_transformation(U)


    use periodic_boundary_contidions_mod, only : ip
    use parameters, only : L

    type(link_variable), dimension(L,L), intent(inout) :: U
    type(complex_2x2_matrix), dimension(L,L) :: V

    integer(i4) :: x, y, mu

    do x = 1, L
       do y = 1, L
          call create_update(V(x,y))
       end do
    end do


    do x = 1, L
       do y = 1, L
          U(x,y)%link(1)%matrix = matmul(V(x,y)%matrix,U(x,y)%link(1)%matrix)
          U(x,y)%link(1)%matrix = matmul(U(x,y)%link(1)%matrix,dagger(V(ip(x),y)%matrix))
          U(x,y)%link(2)%matrix = matmul(V(x,y)%matrix, U(x,y)%link(2)%matrix)
          U(x,y)%link(2)%matrix = matmul(U(x,y)%link(2)%matrix,dagger(V(x,ip(y))%matrix))
       end do
    end do

  end subroutine gauge_transformation

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


  pure function staples(U,x,y,mu)

    use periodic_boundary_contidions_mod, only : ip, im

    type(link_variable), dimension(:,:), intent(in) :: U
    integer(i4), intent(in) :: x, y, mu
    integer(i4) :: nu
    complex(dp), dimension(2,2) :: staples, prod1, prod2, prod3, prod4

    staples = 0.0_dp
    !do nu = 1, 2
    !   if( mu .ne. nu)then
    if(mu == 1) nu = 2
    if(mu == 2) nu = 1
          prod1 = matmul(U(ip(x),y)%link(nu)%matrix,dagger(U(x,ip(y))%link(mu)%matrix))
          prod2 = matmul(prod1,dagger(U(x,y)%link(nu)%matrix))

          prod3 = matmul(dagger(U(ip(x),im(y))%link(nu)%matrix),&
                         dagger(U(x,im(y))%link(mu)%matrix))
          prod4 = matmul(prod3,U(x,im(y))%link(nu)%matrix)

          staples = prod2 + prod4! + staples
    !   end if
    !end do

  end function staples

  pure function DS(U,mu,Up,x,y)
    type(link_variable), dimension(:,:), intent(in) :: U
    integer(i4), intent(in) :: x, y, mu
    integer(i4) :: nu
    type(complex_2x2_matrix), intent(in) :: Up
    real(dp) :: DS

    !call create_update(Up)

    DS = -real(tr(matmul(Up%matrix - U(x,y)%link(mu)%matrix,staples(U,x,y,mu))))

  end function DS


end module dynamics
