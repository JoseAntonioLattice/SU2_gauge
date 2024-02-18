module dynamics

  use iso_fortran_env, only: dp => real64, i4 => int32
  use data_types_observables, only : link_variable, complex_2x2_matrix
  use matrix_operations
  use local_update_algorithms

  implicit none

  !private !:: dp, i4, link_variable
  !public :: sweeps, create_update, sgn, drand, take_measurements, dagger, tr, gauge_transformation, DS

contains

  subroutine sweeps(U,L,beta,N,d,algorithm)

    type(link_variable), intent(inout), dimension(:,:) :: U
    integer(i4), intent(in)  :: L, N, d
    real(dp), intent(in) :: beta
    character(*) :: algorithm
    type(complex_2x2_matrix) :: Up
    integer(i4) :: x, y, mu!, z, w
    real(dp) :: Delta_S

    if( "metropolis" == trim(algorithm))then
      do x = 1, L
        do y = 1, L
        !  do z = 1, L
        !    do w = 1, L
              do mu = 1, d
                 Delta_S = (beta/N) * DS(U,mu,Up,x,y)
                 !Delta_S = DS2(U,x,y,mu,beta/N,d,Up)
                 call metropolis(Delta_S,U(x,y)%link(mu)%matrix,Up%matrix)
        !      end do
        !    end do
          enddo
        end do
      end do
    elseif( "glauber" == trim(algorithm))then
      do x = 1, L
        do y = 1, L
        !  do z = 1, L
        !    do w = 1, L
              do mu = 1, d
                Delta_S = (beta/N) * DS(U,mu,Up,x,y)
                !Delta_S = DS2(U,x,y,mu,beta/N,d,Up)
                call glauber(Delta_S,U(x,y)%link(mu)%matrix,Up%matrix)
              end do
        !    end do
        !  enddo
        end do
      end do
    elseif( "heatbath" == trim(algorithm))then
      do x = 1, L
        do y = 1, L
        !  do z = 1, L
        !    do w = 1, L
              do mu = 1, d
                 call heatbath(U,x,y,mu,beta)
              end do
        !    end do
        !  enddo
        end do
      end do
    else
      error stop "Not a valid algorithm"
    end if

  end subroutine sweeps

  subroutine take_measurements(U,L,N,d,action)

    type(link_variable), dimension(:,:), intent(in) :: U
    integer(i4), intent(in) :: L,N, d
    real(dp), intent(out) :: action
    integer(i4) :: x, y, mu, nu, number_of_planes!, z, w

    action = 0.0_dp
    do x = 1, L
       do y = 1, L
          !do z = 1, L
          !  do w = 1, L
              do mu = 1, d - 1
                do nu = mu + 1, d
                  action = action + real(tr(plaquette(U,x,y,mu,nu)),dp)
              !    action = action + real(tr(plaquette(U(:,y,:,w),x,y,mu,nu)),dp)
              !    action = action + real(tr(plaquette(U(:,y,z,:),x,y,mu,nu)),dp)
              !    action = action + real(tr(plaquette(U(x,:,:,w),x,y,mu,nu)),dp)
              !    action = action + real(tr(plaquette(U(x,:,z,:),x,y,mu,nu)),dp)
              !    action = action + real(tr(plaquette(U(x,y,:,:),x,y,mu,nu)),dp)
                end do
              end do
          !  end do
          !enddo
       end do
    end do

    number_of_planes = d*(d-1)/2

    action =  action!/(N*L**d*number_of_planes)

  end subroutine take_measurements


    function action(U,beta_N,d)

    type(link_variable), dimension(:,:), intent(in) :: U
    integer(i4), intent(in) ::  d
    real(dp) :: action, beta_N
    integer(i4) :: x, y, L!, z, w

    L = size(U(:,1))
    action = 0.0_dp
    do x = 1, L
      do y = 1, L
        action = action + real(tr(plaquette(U,x,y,1,2)))
      end do
    end do
    action =  - beta_N * action
  end function action

  subroutine gauge_transformation(U)


    use periodic_boundary_contidions_mod, only : ip
    use parameters, only : L

    type(link_variable), dimension(L,L), intent(inout) :: U
    type(complex_2x2_matrix), dimension(L,L) :: V

    integer(i4) :: x, y

    do x = 1, L
       do y = 1, L
          call create_update(V(x,y))
       end do
    end do

    do x = 1, L
       do y = 1, L
          U(x,y)%link(1) = V(x,y)*U(x,y)%link(1)*dagger(V(ip(x),y))
          U(x,y)%link(2) = V(x,y)*U(x,y)%link(2)*dagger(V(x,ip(y)))
       end do
    end do

  end subroutine gauge_transformation


 function plaquette(U,x,y,mu,nu)
    use periodic_boundary_contidions_mod, only : ip

    type(link_variable), dimension(:,:), intent(in) :: U
    integer(i4), intent(in) :: x, y, mu, nu
    type(complex_2x2_matrix) :: plaquette

    plaquette = U(x,y)%link(mu) * U(ip(x),y)%link(nu) * dagger(U(x,ip(y))%link(mu)) * dagger(U(x,y)%link(nu))
  end function plaquette

  function DS2(U,x,y,mu,beta_N,d,Up)
    type(link_variable), dimension(:,:), intent(inout) :: U
    integer(i4), intent(in) :: x, y, mu, d
    real(dp), intent(in) :: beta_N
    type(complex_2x2_matrix), intent(out) :: Up
    type(complex_2x2_matrix) :: Uold
    real(dp) :: DS2, Sold, Snew

    call create_unbiased_update(Up)

    Uold = U(x,y)%link(mu)

    Sold = action(U,beta_N,d)
    U(x,y)%link(mu) = Up
    Snew = action(U,beta_N,d)
    U(x,y)%link(mu) = Uold

    DS2 = Snew - Sold

  end function DS2

end module dynamics
