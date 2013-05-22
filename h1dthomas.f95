! 1D Heat Equation Using Thomas Algorithm
! s.t.
!   u(1,t) = 0 for t=>0
!   u(0,t) = 0 for t=>0
!
! Dale Roberts <dale.o.roberts@gmail.com>
!
program h1dthomas
  implicit none

  ! External Functions
  external sgtsv

  ! Model Parameters
  real,    parameter :: th = 0.5

  ! X Space Parameters
  integer, parameter :: nJ = 20
  real,    parameter :: mX = 0
  real,    parameter :: pX = 1

  ! T Space Parameters
  integer, parameter :: nT = 10
  real,    parameter :: mT = 0
  real,    parameter :: pT = 0.1

  ! Variables
  integer :: t, s, info, fd
  real    :: U(nJ+1), b(nJ+1)
  real    :: dl(nJ), dc(nJ+1), du(nJ)
  real    :: x, dx, dt, mu

  dx = (pX-mX) / real(nJ)
  dt = (pT-mT) / real(nT)

  mu = dt / (dx**2)

  ! Set Initial Data
  do s = 1, (nJ+1)
     x = mX + (s-1)*dx
     if (x <= 0.5) then
        U(s) = 2 * x
     else
        U(s) = 2 - 2 * x   
     end if
  end do

  ! Output Data in Gnuplot format
  open(fd, FILE='data.0', STATUS='NEW')

  do s = 1, (nJ+1)
     x = mX + (s-1)*dx
     write (fd, *)  x, U(s)
  end do

  close(fd)

  ! Construct Tri-Diagonals
  dc(1) = 1
  du(1) = 0

  do s = 2, nJ
     du(s)   = -mu * th
     dc(s)   = (1 + 2 * mu * th)
     dl(s-1) = -mu * th
  end do

  dc(nJ+1) = 1
  dl(nJ)   = 0

  do t = 1, nT
     ! Construct RHS for system
     b(1) = 0

     do s = 2, (nJ+1)
        b(s) = (1-th)*mu * U(s-1) + (1-2*mu*(1-th))*U(2) + (1-th)*mu*U(s+1)
     end do

     b(nJ+1) = 0

     ! Solve Tri-Diagonal System
     call sgtsv(nJ+1, 1, dl, dc, du, b, nJ+1, info)

     U = b

     ! Output Data in Gnuplot format
     open(fd, FILE='data.' // char(96 + t), STATUS='NEW')

     do s = 1, (nJ+1)
        x = mX + (s-1)*dx
        write (fd, *)  x, U(s)
     end do

     close(fd)
  end do

end program h1dthomas
