module starts

  use iso_fortran_env, only : dp => real64, i4 => int32

  implicit none

  private
  public :: cold_start
  
contains


  subroutine cold_start(U)
    use data_types_observables, only : link_variable
    type(link_variable), intent(out), dimension(:,:) :: U
    integer(i4) :: L,x,y,z,w,mu
    !Cold start

    L = size(U(:,1))
    
    do x = 1, L
       do y = 1, L
          !do z = 1, L
             !do w =1, L
                do mu = 1, 2
                   U(x,y)%link(mu)%matrix = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp], (/2,2/))
                end do
             !end do
          !end do
       end do
    end do

  end subroutine cold_start


end module starts
