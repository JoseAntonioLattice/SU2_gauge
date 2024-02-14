module periodic_boundary_contidions_mod

    use iso_fortran_env, only : i4 => int32

    implicit none

    integer(i4), allocatable, dimension(:) :: ip, im

    private

    public ip, im, set_periodic_bounds

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


end module periodic_boundary_contidions_mod
