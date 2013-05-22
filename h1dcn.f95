! 2D Heat Equation Using Crank-Nicholson
! s.t.
!   u(1,t) = 0 for t=>0
!   u(0,t) = 0 for t=>0
!
! Dale Roberts <dale.o.roberts@gmail.com
!
program h1dcn
  implicit none

  ! External Functions
  external sgtsv

  ! X Space Parameters
  integer, parameter :: nJ = 40
  real, parameter    :: mX = 0
  real, parameter    :: pX = 1

  ! T Space Parameters
  integer, parameter :: nT = 10
  real, parameter    :: mT = 0
  real, parameter    :: pT = 0.1

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
  do s = 1, nJ
     du(s)   = -mu * 0.5
     dc(s)   = 1 + mu
     dl(s-1) = -mu * 0.5
  end do

  dc(1) = 1
  du(1) = 0
  dl(1) = 0

  du(nJ)   = 0
  dc(nJ+1) = 1
  dl(nJ)   = 0

  do t = 1, nT
     ! Construct RHS for system
     b(1) = 0

     do s = 2, nJ
        b(s) = 0.5 * mu * U(s-1) + (1-mu)*U(s) + 0.5*mu*U(s+1)
     end do

     b(nJ+1) = 0

     ! Solve Tri-Diagonal System
     call sgtsv(nJ+1, 1, dl, dc, du, b, nJ+1, info)

     print '(10F8.2)/', b

     U = b

     ! Output Data in Gnuplot format
     open(fd, FILE='data.' // char(96 + t), STATUS='NEW')

     do s = 1, (nJ+1)
        x = mX + (s-1)*dx
        write (fd, *)  x, U(s)
     end do

     close(fd)
  end do

end program h1dcn
