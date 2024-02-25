module periodic_boundary_conditions_mod

  use iso_fortran_env, only : i4 => int32

  implicit none

  integer, allocatable, dimension(:) :: ip, im

contains

  subroutine set_periodic_bounds(L)

    integer(i4), intent(in) :: L

    allocate(ip(L),im(L))
    call initialize(L)

  end subroutine set_periodic_bounds

  subroutine initialize(L)

    integer(i4), intent(in) :: L
    integer(i4) :: i

    do i = 1, L
       ip(i) = i + 1
       im(i) = i - 1
    end do
    ip(L) = 1
    im(1) = L

  end subroutine initialize

  function ip_func(vector,d)
    integer, dimension(:), intent(in) :: vector
    integer, intent(in) :: d

    integer, dimension(size(vector)) :: ip_func


    ip_func = vector

    ip_func(d) = ip(vector(d))

  end function ip_func


  function im_func(vector,d)

    integer, dimension(:), intent(in) :: vector
    integer, intent(in) :: d

    integer, dimension(size(vector)) :: im_func


    im_func = vector

    im_func(d) = im(vector(d))

  end function im_func

end module periodic_boundary_conditions_mod
