program main
  use, intrinsic :: iso_fortran_env, only: dp=>real64
  use rot_lange, only: say_hello, run_sim
  use stdlib_io_npy, only: save_npy
  use stdlib_io, only: savetxt
  implicit none

  real(dp), allocatable :: t(:), y(:,:)
  integer :: nsteps

  ! call say_hello()
  
  nsteps = 10000

  allocate(t(nsteps+1))
  allocate(y(18, nsteps+1))

  call run_sim(t, y)
  ! call savetxt('t.dat', t)
  call savetxt('y.dat', y)
end program main
