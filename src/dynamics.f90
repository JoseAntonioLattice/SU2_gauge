module dynamics

  use iso_fortran_env, only: dp => real64, i4 => int32
  use data_types_observables, only : link_variable, complex_2x2_matrix
  use matrix_operations
  use local_update_algorithms

  implicit none

  !private !:: dp, i4, link_variable
  !public :: sweeps, create_update, sgn, drand, take_measurements, dagger, tr, gauge_transformation, DS

contains

  subroutine sweeps(U,L,beta,N,algorithm)
    type(link_variable), intent(inout), dimension(:,:) :: U
    integer(i4), intent(in)  :: L
    integer(i4), intent(in) :: N
    real(dp), intent(in) :: beta
    character(*) :: algorithm
    type(complex_2x2_matrix) :: Up
    integer(i4) :: x, y, mu
    real(dp) :: Delta_S

    if( "metropolis" == trim(algorithm))then
      do x = 1, L
         do y = 1, L
            do mu = 1, 2
               Delta_S = (beta/N) * DS(U,mu,Up,x,y)
               call metropolis(Delta_S,U(x,y)%link(mu)%matrix,Up%matrix)
            end do
         end do
      end do
    elseif( "glauber" == trim(algorithm))then
      do x = 1, L
         do y = 1, L
            do mu = 1, 2
               Delta_S = (beta/N) * DS(U,mu,Up,x,y)
               call glauber(Delta_S,U(x,y)%link(mu)%matrix,Up%matrix)
            end do
         end do
      end do
    elseif( "heatbath" == trim(algorithm))then
      do x = 1, L
         do y = 1, L
            do mu = 1, 2
               call heatbath(U,x,y,mu,beta)
            end do
         end do
      end do

    else
      error stop "Not a valid algorithm"
    end if

  end subroutine sweeps

  subroutine take_measurements(U,L,N,d,action)

    type(link_variable), dimension(:,:), intent(in) :: U
    integer(i4), intent(in) :: L,N,d
    real(dp), intent(out) :: action
    integer(i4) :: x, y, mu, nu

    action = 0.0_dp
    do x = 1, L
       do y = 1, L
          do mu = 1, d - 1
            do nu = mu + 1, d
              action = action + real(tr(plaquette(U,x,y,mu,nu)),dp)
            end do
          enddo
       end do
    end do

    !action =  action/(N*L**d)

  end subroutine take_measurements

!  subroutine gauge_transformation(U)
!
!
!    use periodic_boundary_contidions_mod, only : ip
!    use parameters, only : L
!
!    type(link_variable), dimension(L,L), intent(inout) :: U
!    type(complex_2x2_matrix), dimension(L,L) :: V
!
!    integer(i4) :: x, y, mu
!
!    do x = 1, L
!       do y = 1, L
!          call create_update(V(x,y))
!       end do
!    end do
!
!
!    do x = 1, L
!       do y = 1, L
!          U(x,y)%link(1)%matrix = matmul(V(x,y)%matrix,U(x,y)%link(1)%matrix)
!          U(x,y)%link(1)%matrix = matmul(U(x,y)%link(1)%matrix,dagger(V(ip(x),y)%matrix))
!          U(x,y)%link(2)%matrix = matmul(V(x,y)%matrix, U(x,y)%link(2)%matrix)
!          U(x,y)%link(2)%matrix = matmul(U(x,y)%link(2)%matrix,dagger(V(x,ip(y))%matrix))
!       end do
!    end do
!
!  end subroutine gauge_transformation


 function plaquette(U,x,y,mu,nu)
    use periodic_boundary_contidions_mod, only : ip

    type(link_variable), dimension(:,:), intent(in) :: U
    integer(i4), intent(in) :: x, y, mu, nu
    type(complex_2x2_matrix) :: plaquette

    plaquette = U(x,y)%link(mu) * U(ip(x),y)%link(nu) * dagger(U(x,ip(y))%link(mu)) * dagger(U(x,y)%link(nu))
  end function plaquette

end module dynamics
