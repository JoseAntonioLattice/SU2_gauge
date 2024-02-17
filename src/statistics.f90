module statistics

  use iso_fortran_env, dp => real64, i4 => int32


  implicit none

  private
  public :: std_err, std_err_scalar, jackknife, max_jackknife_error,max_jackknife_error_2

contains


  subroutine std_err(array, average, error)

    implicit none
    real(dp), intent(in), dimension(:) :: array
    real(dp), intent(out) :: average, error

    real(dp) :: variance2
    integer(i4) :: n

      n = size(array)

      average = sum(array)/n
      variance2 = (sum(array**2) - n*average**2)/(n-1)
      error = sqrt(variance2/n)

  end subroutine std_err


  subroutine jackknife(array,average,error,bins)

    real(dp), intent(in), dimension(:) :: array
    real(dp), intent(out) :: average, error
    integer(i4), intent(in) :: bins

    real(dp), dimension(bins) :: theta
    integer(i4) :: n,m
    integer(i4) :: i

      n = size(array)
      m = n/bins

      theta = [(sum(array) - sum(array(m*(i-1)+1:i*m)), i = 1, bins )]/(n-m)
      average = sum(array)/n
      error = sqrt((bins-1)/real(bins,dp)*sum((theta-average)**2))

  end subroutine jackknife


  subroutine max_jackknife_error(array,average,error,bins)

    real(dp), intent(in) :: array(:)
    real(dp), intent(out) :: average, error
    integer(i4), intent(out) :: bins
    integer(i4), allocatable, dimension(:) :: bin_array
    real(dp), allocatable, dimension(:) :: jackk_error_array
    integer(i4) :: i, n_bin, mxloc

    call multiples(size(array),bin_array)

    n_bin = size(bin_array)

    allocate(jackk_error_array(n_bin))

    do i = 1, n_bin
      call jackknife(array,average,error,bin_array(i))
      jackk_error_array(i) = error
      print*, i, average, error
    end do

    mxloc = maxloc(jackk_error_array,dim=1)
    bins  = bin_array(mxloc)
    error = maxval(jackk_error_array,dim=1)


    deallocate(bin_array)

  end subroutine max_jackknife_error


  subroutine max_jackknife_error_2(array,average,error,bins)

    real(dp), intent(in) :: array(:)
    real(dp), intent(out) :: average, error
    integer(i4), intent(out) :: bins

    real(dp),dimension(size(array)) :: jackk_error_array, jackk_average_array

    integer(i4) :: i, n, m

      n = size(array)

      do i = 2, n
        m = mod(n,i)
        call jackknife(array(m+1:n),average,error,i)
        jackk_error_array(i) = error
        jackk_average_array(i) = average
        print*, i, average, error, size(array(m+1:n))
      end do

      bins = maxloc(jackk_error_array,dim = 1)
      average = jackk_average_array(bins)
      error = maxval(jackk_error_array,dim = 1)

  end subroutine max_jackknife_error_2


  subroutine multiples(n,y)
    integer(i4), intent(in) :: n
    integer(i4), intent(out), allocatable, dimension(:) :: y
    integer(i4) :: i, counter

    counter = 0
    do i = 2, n
      if(mod(n,i) == 0) counter = counter + 1
    end do

   allocate(y(counter))
    counter = 0
    do i = 2, n
      if(mod(n,i) == 0)then
       counter = counter + 1
       y(counter) = i
      end if
    end do

  end subroutine multiples


    subroutine STD_ERR_scalar(summa,SUMMASQ,N,avg,std_ERROR)
    IMPLICIT NONE
    real(dp), INTENT(IN)  :: SUMMA, SUMMASQ
    REAL(dp), INTENT(OUT) :: AVG          ! SAMPLE MEAN
    REAL(dp), INTENT(OUT) :: STD_ERROR    ! standard error
    INTEGER(i4), INTENT(IN) :: N
    REAL(dp) :: VAR  ! VARIANCE
    real(dp) :: S    ! STANDARD DEVIATION OF THE SAMPLE MEAN

    AVG = SUMMA/DBLE(N)
    ! VariaNCE
    var = (SUMMASQ - N*AVG*AVG)/(N - 1)

    ! STANDARD ERROR
    S = sqrt(var)
    STD_ERROR = S/N**0.5D0

    end subroutine STD_ERR_scalar


end module statistics
