program fput 
    use utility, only: fp, pi
    use problemSetup, only: alpha, K, dt, setup_init, set_ics, nMasses, nSteps, outFile
    use leapfrog, only: leap_frog, leapstart
    
    implicit none
    integer :: i,j
    real(fp), allocatable :: xn(:),xn1(:),xn_1(:)
    real(fp), allocatable :: v0(:)

    call setup_init()
    allocate(xn_1(nMasses+2),xn1(nMasses+2),xn(nMasses+2),v0(nMasses+2))

    call set_ics(xn_1,v0)
    xn_1(1) = 0.0_fp
    xn_1(nMasses+2) = 0.0_fp

    open(20, file = outFile, status = "replace")
    write(20,*) 0.0_fp, (xn_1(i), i=1, nMasses+2) !timestep 0, all masses including walls
    

    call leapstart(xn_1,dt,v0,xn)
    xn(1) = 0.0_fp
    xn(nMasses+2) = 0.0_fp
    write(20,*) dt, (xn(i), i=1, nMasses+2) !timestep 1, all masses including walls
    do i = 2,nSteps
        call leap_frog(xn,xn1,xn_1,dt,K,alpha)
        write(20,*) i * dt, (xn1(j), j=1, nMasses+2) !timestep n, all masses including walls
        xn_1 = xn
        xn = xn1
        
    end do
    close(20)
    deallocate(xn_1,xn1,xn,v0)
    
    

    
    !open
    !open(20, file=outFile,status="replace")
    !write(20,*) (xn_1(i),i=1,nMasses)
    !number stuff

    
    !close
 


   
  
    
    
   


    





contains



end program fput