module parameters

  use iso_fortran_env, only: dp => real64, i4 => int32

  implicit none

  private :: dp, i4

  integer(i4) :: L ! Lattice length size
  integer(i4) :: N ! Dimension of gauge group
  integer(i4) :: d
  integer(i4) :: N_measurements
  integer(i4) :: N_thermalization
  integer(i4) :: N_skip
  character(20) :: algorithm


  namelist /input_parameters/ L,N, d,N_thermalization, N_measurements, N_skip, algorithm

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
    if( N <= 0) error stop "Dimension of gauge group N must be > 0."
    if( N_thermalization <= 0) error stop "Thermalization sweeps must be > 0."
    if( N_measurements <= 0) error stop "Number of measurements must be > 0."

    write(stdout, nml = input_parameters)


  end subroutine read_input_parameters

end module parameters
