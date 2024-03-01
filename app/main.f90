program main
  use, intrinsic :: iso_fortran_env, only: dp=>real64
  use rot_lange, only: say_hello, run_sim
  use stdlib_io_npy, only: save_npy
  use stdlib_io, only: savetxt
  implicit none

  real(dp), allocatable :: t(:), y(:,:)
  real(dp) :: dt
  integer :: nsteps, nsteps_out

  ! call say_hello()
  
  ! nsteps = 10000000
  nsteps = 100

  ! nsteps_out = 1000
  nsteps_out = 1

  dt = 0.0001_dp

  allocate(t(nsteps+1))
  allocate(y(18, nsteps+1))

  call run_sim(t, y, dt, nsteps_out)
  ! call savetxt('t.dat', t)
  ! call savetxt('y.dat', transpose(y))
end program main
