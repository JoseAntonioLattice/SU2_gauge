module dynamics

  use iso_fortran_env, only: dp => real64, i4 => int32
  use data_types_observables, only : link_variable, complex_2x2_matrix
  use matrix_operations
  use local_update_algorithms
  use periodic_boundary_conditions_mod
  use get_index_mod

  implicit none

  !private !:: dp, i4, link_variable
  !public :: sweeps, create_update, sgn, drand, take_measurements, dagger, tr, gauge_transformation, DS

contains

  subroutine thermalization(U,L,beta,N,d,algorithm,N_thermalization)
    type(link_variable), intent(inout), dimension(:) :: U
    integer(i4), intent(in)  :: L, N, d
    real(dp), intent(in) :: beta
    character(*), intent(in) :: algorithm
    integer(i4), intent(in) :: N_thermalization

    integer(i4) :: i

    do i = 1, N_thermalization
       call sweeps(U,L,beta,N,d,algorithm)
    end do
   end subroutine thermalization


   subroutine measurements_sweeps(U,L,beta,N,d,algorithm,N_measurements,N_skip,E_p)
     type(link_variable), intent(inout), dimension(:) :: U
     integer(i4), intent(in)  :: L, N, d
     real(dp), intent(in) :: beta
     character(*), intent(in) :: algorithm
     integer(i4), intent(in) :: N_measurements, N_skip
     real(dp), intent(out) :: E_p(:)
     integer(i4) :: i

     do i = 1, N_measurements*N_skip
        call sweeps(U,L,beta,N,d,algorithm)
        if( mod(i,N_skip) == 0)then
           E_p(i/N_skip) = action(U,-1.0_dp/N,d)/(L**d)
        end if
     end do

   end subroutine measurements_sweeps

  subroutine sweeps(U,L,beta,N,d,algorithm)
    type(link_variable), intent(inout), dimension(:) :: U
    integer(i4), intent(in)  :: L, N, d
    real(dp), intent(in) :: beta
    character(*), intent(in) :: algorithm
    type(complex_2x2_matrix) :: Up
    integer(i4) :: x, y, mu!, z, w
    real(dp) :: Delta_S

    do x = 1, L**d
       do mu = 1, d
          Delta_S = (beta/N) * DS(U,mu,Up,x)
          !Delta_S = DS2(U,x,y,mu,beta/N,d,Up)
          !print*, "Inside sweeps", x, mu, Delta_S
          call metropolis(Delta_S,U(x)%link(mu)%matrix,Up%matrix)
          !print*, "after metropolis"
       end do
    end do


  end subroutine sweeps

  subroutine take_measurements(U,L,N,d,action)

    type(link_variable), dimension(:), intent(in) :: U
    integer(i4), intent(in) :: L,N, d
    real(dp), intent(out) :: action
    integer(i4) :: x, y, mu, nu, number_of_planes!, z, w

    action = 0.0_dp
    do x = 1, L
       do mu = 1, d - 1
          do nu = mu + 1, d
             action = action + real(tr(plaquette(U,x,mu,nu)),dp)
          end do
       end do
    end do

    number_of_planes = d*(d-1)/2

    action =  action!/(N*L**d*number_of_planes)

  end subroutine take_measurements


    function action(U,beta_N,d)

    type(link_variable), dimension(:), intent(in) :: U
    integer(i4), intent(in) ::  d
    real(dp) :: action, beta_N
    integer(i4) :: x, V, mu,nu, number_of_planes

    V = size(U)
    action = 0.0_dp
    do x = 1, V
        do mu = 1, d - 1
           do nu = mu + 1, d
              !print*, "Inside action. Inside loop", x, mu, nu
             action = action + real(tr(plaquette(U,x,mu,nu)),dp)
          end do
       end do
    end do
    !print*, "Inside action. outside loop"
    number_of_planes = d*(d-1)/2
    action =  - beta_N * action/number_of_planes
  end function action

  !subroutine gauge_transformation(U)


    !use periodic_boundary_contidions_mod, only : ip
    !use parameters, only : L

   ! type(link_variable), dimension(L), intent(inout) :: U
  !  type(complex_2x2_matrix), dimension(L,L) :: V

 !   integer(i4) :: x, y

!    do x = 1, L
    !   do y = 1, L
    !      call create_update(V(x,y))
    !   end do
    !end do

   ! do x = 1, L
   !    do y = 1, L
    !      U(x,y)%link(1) = V(x,y)*U(x,y)%link(1)*dagger(V(ip(x),y))
   !       U(x,y)%link(2) = V(x,y)*U(x,y)%link(2)*dagger(V(x,ip(y)))
   !    end do
   ! end do

  !end subroutine gauge_transformation


 function plaquette(U,x,mu,nu)
    use parameters, only : d, L
    type(link_variable), dimension(:), intent(in) :: U
    integer(i4), intent(in) :: x, mu, nu
    type(complex_2x2_matrix) :: plaquette

    integer(i4) :: ipx_mu, ipx_nu
    !             x, mu           x + mu, nu                   x + nu, mu                    x, nu
    !print*, "Inside plaquette"
    ipx_mu = get_index(ip_func(get_index_array(x,d,L),mu),d,L)
    ipx_nu = get_index(ip_func(get_index_array(x,d,L),nu),d,L)
    plaquette = U(x)%link(mu) * U(ipx_mu)%link(nu) * dagger(U(ipx_nu)%link(mu)) * dagger(U(x)%link(nu))
  end function plaquette

!  function DS2(U,x,mu,beta_N,d,Up)
!    type(link_variable), dimension(:,:), intent(inout) :: U
!    integer(i4), intent(in) :: x, mu, d
!    real(dp), intent(in) :: beta_N
    !type(complex_2x2_matrix), intent(out) :: Up
   ! type(complex_2x2_matrix) :: Uold
  !  real(dp) :: DS2, Sold, Snew

 !   call create_unbiased_update(Up)

!    Uold = U(x)%link(mu)

    !Sold = action(U,beta_N,d)
    !U(x)%link(mu) = Up
    !Snew = action(U,beta_N,d)
    !U(x)%link(mu) = Uold

   ! DS2 = Snew - Sold

  !end function DS2

end module dynamics
