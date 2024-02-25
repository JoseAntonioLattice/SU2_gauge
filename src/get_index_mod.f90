module get_index_mod

  implicit none

contains

  function get_index(vector,dimension_space,length_lattice)

    integer :: get_index

    integer, dimension(:), intent(in) :: vector
    integer, intent(in) :: dimension_space
    integer, intent(in) :: length_lattice

    integer :: suma
    integer :: i


    suma = vector(1)

    if( dimension_space > 1)then
       do i = 2, dimension_space
          suma = suma + (vector(i) - 1) * length_lattice**(i-1)
       end do
    end if

    get_index = suma


  end function get_index


  function get_index_array(idx,d,L) result(vector)

    integer, intent(in) :: idx
    integer, intent(in) :: d
    integer, intent(in) :: L

    integer, dimension(d) :: vector

    integer :: i, n, modx, suma

    modx = mod(idx,L)
    vector(1) = modx
    if(modx == 0) vector(1) = L

    suma = vector(1)
    do i = 2, d
      modx = mod(idx,L**i)
      if (i > 2) suma = suma + L**(i-2)*(vector(i-1)-1)
      vector(i) = (modx - suma)/L**(i-1) + 1
      if(modx == 0) vector(i) = L
    end do

  end function get_index_array


end module get_index_mod
