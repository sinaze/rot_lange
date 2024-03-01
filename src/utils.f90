module utils
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    implicit none
    
contains

    pure function matinv3(A) result(B)
        real(dp), intent(in) :: A(3,3)
        real(dp) :: B(3,3)
        real(dp) :: detinv

        ! Calculate the inverse determinant of the matrix
        detinv = 1.0_dp/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
                  - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
                  + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

        ! Calculate the inverse of the matrix
        B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
        B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
        B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
        B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
        B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
        B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
        B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
        B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
        B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    end function
    
    pure function star(a) result(M)
        ! Return the skew-symmetric matrix M of vector a
        real(dp), intent(in) :: a(3)
        real(dp) :: M(3,3)

        M(1,1) = 0
        M(1,2) = -a(3)
        M(1,3) = a(2)
        M(2,1) = a(3)
        M(2,2) = 0
        M(2,3) = -a(1)
        M(3,1) = -a(2)
        M(3,2) = a(1)
        M(3,3) = 0
    end function

end module utils
