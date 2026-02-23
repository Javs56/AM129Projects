! File: discreteFourierTransform.f90
! Purpose: Provides discrete Fourier transformation matrices

module discreteFourierTransform
    
    use utility, only : fp, pi
    
    implicit none
    
contains
    
    !! !! === 1. function: matvecprod
    function matvecprod(A,x) result(y)
        real (fp), intent(in) :: A(:,:) !arbitrary size, rank 2
        real (fp), intent(in) :: x(:) !arbitrary size, rank 1
        real (fp), dimension(size(A,1)) :: y !arbitrary size, rank 1
                                 !I think in this case y is a row vector?
        integer :: j             !or at least column vector spelt horizontally (array)  
                                   
        y = 0.0_fp !initialize y as a 0 vector
        do j = 1, size(A,2) ! for 1 to the number of columns in A
            y = y + A(:,j)*x(j) !! y = sum A(:,j)*x(j)
        end do

        
    end function matvecprod
    ! function: discreteFourierTransform_TransMat
    ! purpose: Fill transformation matrix for a discrete Fourier transform
    !          on a given domain, and return compatible set of wavenumbers
    subroutine discreteFourierTransform_TransMat(x,k,T)
        implicit none
        real (fp), intent(in)     :: x(:)
        real (fp), intent(out)    :: k(size(x))
        real (fp), intent(in out) :: T(size(x),size(x))
        
        ! Local variables
        integer :: N, i
        real (fp) :: om, dx
        
        ! Set sizes and base wavenumber
        N = size(x)
        print *, N
        dx = x(2) - x(1)
        om = 2*pi/(N*dx) !w = 2pi / Ndx
         
        
        ! Set wavenumbers (array)
        k(1) = 0.0_fp !first wavenumber = [0,....,]
        do i=2,N,2 !from 2 to len(k), steps of 2
            k(i) = i*om/2 !2om/2, 4om/2, 6om/2, ....
            if (i+1 <= N) then
                k(i+1) = k(i)
            end if
        end do ! k = (0, w, w, 2w, 2w, 3w, 3w, ....)
        
        !! !! === 2. Add your code to fill T here
        T(1,:) = 1.0/ N !row 1, all columns = 1/N

        do i=2,N,2 !from 2 to size(x), increment 2 == EVEN NUMBERS
            T(i,:) = (2.0/N) * cos(k(i) * x(:)) 
        end do

        do i = 3,N,2 !from 3 to size(x), increment 2 == ODD NUMBERS
            T(i,:) = (2.0/N) * sin(k(i) * x(:))
        end do

    end subroutine discreteFourierTransform_TransMat
    
    !! !! === 3. subroutine: discreteFourierTransform_InvTransMat
    subroutine discreteFourierTransform_InvTransMat(x,k,Tinv)
        implicit none
        real (fp), intent(in)     :: x(:)
        real (fp), intent(in)     :: k(size(x))
        real (fp), intent(in out) :: Tinv(size(x),size(x))

        integer :: N, j 
        N = size(x)
        Tinv(:,1) = 1.0_fp !all rows, first column = 1s
        do j=2,N,2
            Tinv(:,j) = cos(k(j) * x(:))
        end do

        do j=3,N,2
            Tinv(:,j) = sin(k(j) * x(:))
        end do

    end subroutine discreteFourierTransform_InvTransMat
end module discreteFourierTransform
