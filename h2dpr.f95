! 2D Heat Equation Using Peaceman-Rachford
!
! Dale Roberts <dale.o.roberts@gmail.com>
!
program h2dpr
  implicit none

  ! External Functions
  external sgtsv

  ! X Space Parameters
  integer, parameter :: nJ = 40
  real,    parameter :: mX = 0
  real,    parameter :: pX = 1

  ! Y Space Parameters
  integer, parameter :: nK = 40
  real,    parameter :: mY = 0
  real,    parameter :: pY = 1

  ! T Space Parameters
  integer, parameter :: nT = 20
  real,    parameter :: mT = 0
  real,    parameter :: pT = 0.4

  ! Variables
  integer :: t, s, r, info, fd
  real    :: U((nJ+1)*(nK+1)), b((nJ+1)*(nK+1))
  real    :: dl((nJ+1)*(nK+1)-1), dc((nJ+1)*(nK+1)), du((nJ+1)*(nK+1)-1)
  real    :: x, y, dx, dy, dt, mux, muy

  dx = (pX-mX) / real(nJ)
  dy = (pY-mY) / real(nK)
  dt = (pT-mT) / real(nT)

  mux = dt / (dx**2)
  muy = dt / (dy**2)

  ! Set Initial Data
  do s = 1, (nJ+1)
     do r = 1, (nK+1)
        x = mX + (s-1)*dx
        y = mY + (r-1)*dy
        U((s-1)*(nJ+1)+r) = max(0, 1.0/3.0 - sqrt((x-0.5)**2 + (y-0.5)**2))
     end do
  end do

  ! EXERCISE

  ! Output Data in Gnuplot format
  open(fd, FILE='data.0', STATUS='NEW')

  do s = 1, (nJ+1)
     x = mX + (s-1)*dx
     do r = 1, (nK+1)
        y = mY + (r-1)*dy
        write (fd, *)  x, y, U((s-1)*(nJ+1)+r)
     end do
     write (fd, *) ''
  end do



  close(fd)
end program h2dpr
