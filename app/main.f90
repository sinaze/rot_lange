program main
  use, intrinsic :: iso_fortran_env, only: dp=>real64
  use rot_lange, only: say_hello, run_sim
  use stdlib_io_npy, only: save_npy
  use stdlib_io, only: savetxt
  implicit none

  real(dp), allocatable :: t(:), y(:,:)
  real(dp) :: dt
  integer :: nsteps, nsteps_out
  character(len=100) :: filename

  ! call say_hello()
  
  ! 1 ns
  ! nsteps = 10000000
  ! 10 ns
  nsteps = 100000000
  ! nsteps = 1000

  ! 1 ps
  nsteps_out = 1000
  ! nsteps_out = 100

  dt = 0.0001_dp

  allocate(t(nsteps+1))
  allocate(y(18, nsteps+1))

  call get_command_argument(1, filename)

  call run_sim(t, y, dt, nsteps_out, filename)
  ! call savetxt('t.dat', t)
  ! call savetxt('y.dat', transpose(y))
end program main
