! File: Orbits.f90
! Author: Ian May

program orbits

  use utility, only: fp
  use timestep, only: take_step

  implicit none
  
  integer, parameter :: nSteps = 1000  ! Number of time steps to use
  integer :: nT                        ! Loop variable for time updates
  !!! Step 1

  ! Give all particles initial mass, position, and momentum

real(fp) :: tFinal !declaring tFinal, real
real(fp) :: dt !declaring dt, real
real(fp), dimension(2) :: mass !declaring mass rank 1 real array 2 elements
real(fp), dimension(2,2,nSteps) :: pos, mom !declaring pos,mom rank 3 real arrays 2x2

tFinal = 50
dt = tFinal/nSteps !50/1000 = 0.05 unit width
  call set_ics()

  ! Fill rest of the array by integrating in time
  do nT=1,nSteps-1
    call take_step(dt,mass,pos(:,:,nT),mom(:,:,nT),pos(:,:,nT + 1),mom(:,:,nT + 1))
  end do

  ! Write the arrays to a file for plotting
  call write_data()

contains

  subroutine set_ics()
    implicit none
    !!! Step 2
    ! Set particle masses
    mass = (/1.0_fp,0.01_fp/)
              !m1      m2
    
    ! First particle position and momentum
    pos(:,1,1) = (/0.0_fp,0.0_fp/) !at the origin
    mom(:,1,1) = mass(1)*(/0.0_fp,0.00_fp/) ! 0 momentum

    !(xy, no. of particle, starting timestep)

    ! Second particle position and momentum
    pos(:,2,1) = (/0.0_fp,-1.0_fp/) ! at x = (0,-1)
    mom(:,2,1) = mass(2)*(/1.0_fp,0.0_fp/) ! p2 = m2v2 where v2 = (1,0)
    
  end subroutine set_ics

  subroutine write_data()
    implicit none
    open(20,file = "sol.dat",status = "replace")
    do nT=1,nSteps
      write(20,*) pos(:,:,nT),mom(:,:,nT)
    end do
    close(20)
  end subroutine write_data

end program orbits
