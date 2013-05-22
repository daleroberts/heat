! 2D Heat Equation Using LOD Method
!
! Dale Roberts <dale.o.roberts@gmail.com>
!
program h2dlod
  implicit none

  ! External Functions
  external sgtsv

  ! X Space Parameters
  integer, parameter :: nJ = 40
  real,    parameter :: mX = -5
  real,    parameter :: pX =  5

  ! Y Space Parameters
  integer, parameter :: nK = 40
  real,    parameter :: mY = -5
  real,    parameter :: pY =  5

  ! T Space Parameters
  integer, parameter :: nT = 20
  real,    parameter :: mT = 0
  real,    parameter :: pT = 0.4

  ! Variables
  integer :: t, s, r, info, fd
  real    :: U(nJ+1, nK+1), b(nJ+1)
  real    :: dl(nJ), dc(nJ+1), du(nJ)
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
        U(s, r) = exp( - (x**2 + y**2))
     end do
  end do

  ! Construct Tri-Diagonals
  dc(1) = (1+mux)
  du(1) = -0.5 * mux

  do s = 2, nJ
     dc(s)   = (1+mux)
     du(s)   = -0.5 * mux
     dl(s-1) = -0.5 * mux
  end do

  dc(nJ+1) = (1+mux)
  dl(nJ)   = -0.5 * mux



  do t = 1, nT
     ! Generate U^{n+*}
     do r = 1, (nK+1)
        ! Construct RHS of linear system
        b(1) = (1 - mux) * U(1,r) + 0.5 * mux * U(2,r)

        do s = 2, nJ
           b(s) = 0.5 * mux * U(s-1,r) + (1 - mux) * U(s,r) + 0.5 * mux * U(s+1,r)
        end do

        b(nJ+1) = 0.5 * mux * U(nJ,r) + (1 - mux) * U(nJ+1,r) 

        ! Solve Tri-Diagonal System
        call sgtsv(nJ+1, 1, dl, dc, du, b, nJ+1, info)

        do s = 1, (nJ+1)
           U(s, r) = b(s)
        end do
     end do

     ! Generate U^{n+1}
     do s = 1, (nJ+1)
        ! Construct RHS of linear system
        b(1) = (1 - muy) * U(s,1) + 0.5 * muy * U(s,2)

        do r = 2, nJ
           b(r) = 0.5 * muy * U(s,r-1) + (1 - muy) * U(s,r) + 0.5 * mux * U(s,r+1)
        end do

        b(nK+1) = 0.5 * muy * U(s,nK) + (1 - muy) * U(s,nK+1) 

        ! Solve Tri-Diagonal System
        call sgtsv(nK+1, 1, dl, dc, du, b, nK+1, info)

        do r = 1, (nJ+1)
           U(s, r) = b(r)
        end do
     end do

     ! Output Data in Gnuplot format
     open(fd, FILE='data.' // char(96 + t), STATUS='NEW')
     
     do s = 1, (nJ+1)
        x = mX + (s-1)*dx
        do r = 1, (nK+1)
           y = mY + (r-1)*dy
           write (fd, *)  x, y, U(s,r)
        end do
        write (fd, *) ''
     end do

     close(fd)

  end do
end program h2dlod
