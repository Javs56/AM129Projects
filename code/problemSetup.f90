module problemSetup

  use utility, only : fp, maxFileLen, maxStrLen, pi
  use read_initFile
  
  implicit none
  private

  integer, public :: nMasses, nSteps
  real (fp), public :: K, alpha, tFinal, dt, C
  character(len=maxStrLen), public :: runName, outFile
  character(len=maxStrLen) :: inFile

  public :: setup_init
  public :: set_ics

contains

  subroutine setup_init()
    implicit none

    ! Get name of input file
    call get_command_argument(1,inFile)
    print *, "Reading from ",inFile
    
    ! Fill in default values
    nMasses = 1
    K = 1.0_fp
    alpha = 0.0_fp
    C = 1.0
    
    ! Read problem settings from the input file
    nMasses = read_initFileInt(infile, 'num_masses')
    alpha = read_initFileReal(infile, 'alpha')
    C = read_initFileReal(infile, 'C')
    tFinal = pi * read_initFIleReal(infile, 'tFinal(pi)')
    !!! ====================== Calculate K, nSteps, and dt here ===============================
    K = 4.0_fp * (nMasses + 1)**2
    nSteps = ceiling((tFinal * SQRT(K)) / C)
    dt = tFinal / nSteps 


    ! Set the name of the run and echo it out
    runName = read_initFileChar(infile, 'run_name')
    print *, 'Running problem: ', runName
    
    ! Set the output file, note that // does string concatenation
    outFile = 'data/' // trim(runName) // '.dat'
  end subroutine setup_init

  !!! ============================= Add set_ics subroutine here ===============================

  subroutine set_ics(x,v)
    implicit none

    real (fp), dimension(:), intent(out) :: x,v !input variable x (an array), v array for velocity
    integer :: i !indexing variable
    
    x(1) = 0.0_fp
    x(size(x)) = 0.0_fp !first wall = x=0,v=0
    v(1) = 0.0_fp
    v(size(x)) = 0.0_fp !last wall x=0,v=0

    !ics from 2.4
    do i = 1,nMasses
      x(i+1) = 0.0_fp !maybe name x_0
    end do
    
    do i = 1,nMasses
      v(i+1) = sin(real((i*pi)) / real((nMasses + 1)))
    end do 
 
  end subroutine set_ics


end module problemSetup
