module parameters

  use iso_fortran_env, only: dp => real64, i4 => int32
  
  implicit none

  private :: dp, i4
  
  integer(i4) :: L ! Lattice length size



  namelist /input_parameters/ L
  
contains

  subroutine read_input_parameters()
    use iso_fortran_env, only : stdout => output_unit, stdin => input_unit, stderr => error_unit 
    character(100) :: parameters_file
    integer(i4) :: inunit
    
    
    write(stdout,'(a)') "Enter the parameters file, please."
    read(stdin,'(a)') parameters_file
    write(stdout,'(a)') "User typed: ", trim(parameters_file)

    open(newunit = inunit, file = trim(parameters_file), status = 'old')
    read(inunit, nml = input_parameters)
    if( L <= 0) error stop "Lattice length L must be > 0."
   
    write(stdout, nml = input_parameters)
    

  end subroutine read_input_parameters
  
end module parameters
