module matrix_operations

contains


  pure function det(A)
    complex(dp), dimension(2,2), intent(in) :: A
    real(dp) :: det

    det = real(A(1,1)*A(2,2) - A(1,2)*A(2,1),dp)

  end function det

  pure function dagger(U) result(U_res)
    complex(dp), dimension(:,:), intent(in) :: U
    complex(dp), dimension(size(U(:,1)),size(U(1,:))) :: U_res

    U_res = conjg(U)
    U_res = transpose(U_res)

  end function dagger

  pure function tr(U) result(trace)
    complex(dp), dimension(:,:), intent(in) :: U
    complex(dp) :: trace
    integer :: i

    trace = sum([(U(i,i),i = 1, size(U(1,:)))])

  end function tr


end module matrix_operations
