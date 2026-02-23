module leapfrog

    use utility, only: fp, maxFileLen, maxStrLen, pi

    public :: leap_frog !formula to move forward in time
    public :: leapstart !formula for x1 given time 0 info

contains                          
      
    subroutine leapstart(x0,dt,v0,x1)
        implicit none
        real (fp), intent(in) :: x0(:),v0(:), dt
        real (fp), intent(out) :: x1(:)
        integer :: i
        do i = 1,size(x0)
            x1(i) = x0(i) + dt * v0(i) !initial formula 1.8 our starting point
        end do
    end subroutine leapstart

    subroutine leap_frog(xn,xn1,xn_1,dt,K, alpha)               !xn   = x^n present step time 1
        implicit none                                           !xn1  = x^n+1 future step time +1
        real (fp), intent(in) :: xn(:),xn_1(:),K, dt, alpha     !xn_1 = x^n-1 past step time 0
        real (fp), intent(out) :: xn1(:)
        integer :: i

        xn1(1) = 0.0_fp !first wall = 0
        xn1(size(xn)) = 0.0_fp !last wall = 0
        do i = 2,size(xn) - 1 !I think size(xn) is right here, dongwook just put N
            xn1(i) = 2.0_fp*xn(i) - xn_1(i) + (K*(dt**2))*(xn(i+1) - 2.0_fp*xn(i) + xn(i-1))*(1.0_fp + alpha*(xn(i+1) - xn(i-1))) !formula 2.1
        end do                                                                        

    end subroutine leap_frog

    
end module leapfrog
