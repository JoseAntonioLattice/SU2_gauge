module dynamics

  use iso_fortran_env, only: dp => real64, i4 => int32
  use data_types_observables, only : link_variable, complex_2x2_matrix

  implicit none

  private !:: dp, i4, link_variable
  public :: sweeps, create_update, sgn, drand, take_measurements, dagger, tr, gauge_transformation, DS

contains

  subroutine sweeps(U,L)
    type(link_variable), intent(inout), dimension(:,:) :: U
    integer(i4), intent(in)  :: L
    type(complex_2x2_matrix) :: Up
    integer(i4) :: x, y, z, w, mu
    real(dp) :: Delta_S

    do x = 1, L
       do y = 1, L
          do mu = 1, 2
             call create_update(Up)
             Delta_S = DS(U,mu,Up,x,y)
             call glauber(Delta_S,U(x,y)%link(mu)%matrix,Up%matrix)
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
  
  subroutine create_update(Up)
    type(complex_2x2_matrix), intent(out) :: Up
    real(dp) :: x0, x1, x2, x3
    complex(dp) :: a, b
    real(dp), parameter :: eps = 0.5
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


  subroutine take_measurements(U,action)
    use parameters, only : L
    

    type(link_variable), dimension(:,:), intent(in) :: U
    integer(i4) :: x,y!,z,w,mu
    complex(dp), dimension(2,2), parameter :: one = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp], [2,2])
    real(dp), intent(out) :: action
    action = 0.0_dp
    do x = 1, L
       do y = 1, L
          action = action - real(tr(plaquette(U,x,y,1,2)))
       end do
    end do

    action = 2.0_dp*L**2 + action

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


  pure function staples(U,x,y,mu,nu)

    use periodic_boundary_contidions_mod, only : ip, im
    
    type(link_variable), dimension(:,:), intent(in) :: U
    integer(i4), intent(in) :: x, y, mu, nu
    complex(dp), dimension(2,2) :: staples, prod1, prod2, prod3, prod4

    prod1 = matmul(U(ip(x),y)%link(nu)%matrix,dagger(U(x,ip(y))%link(mu)%matrix))
    prod2 = matmul(prod1,dagger(U(x,y)%link(nu)%matrix))

    prod3 = matmul(dagger(U(ip(x),im(y))%link(nu)%matrix),&
                   dagger(U(x,im(y))%link(mu)%matrix))
    prod4 = matmul(prod3,U(x,im(y))%link(nu)%matrix)
    
    staples = prod2 + prod4
    
  end function staples

  pure function DS(U,mu,Up,x,y)
    type(link_variable), dimension(:,:), intent(in) :: U
    integer(i4), intent(in) :: x, y, mu
    integer(i4) :: nu
    type(complex_2x2_matrix), intent(in) :: Up
    real(dp) :: DS
    
    !call create_update(Up)

    if(mu == 1) nu = 2
    if(mu == 2) nu = 1
    
    DS = -real(tr(matmul(Up%matrix - U(x,y)%link(mu)%matrix,staples(U,x,y,mu,nu))))
    
  end function DS
  

end module dynamics
