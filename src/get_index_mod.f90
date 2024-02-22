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


  function get_index_array(idx,d,L)

    integer, intent(in) :: idx
    integer, intent(in) :: d
    integer, intent(in) :: L
    
    integer, dimension(d) :: get_index_array
    
    integer :: i, n
    
    get_index_array(1) = mod(idx,L) 
    if(mod(idx,L) == 0) get_index_array(1) = L 
    
    do i = 2, d
       get_index_array(i) =  idx/L**(d-1) + 1
       if(mod(idx,L) == 0) get_index_array(i) =  idx/L**(d-1) 
    end do
    
  end function get_index_array
  
  
end module get_index_mod
