module rot_lange
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options
  use stdlib_linalg, only: eye
  use types, only: t_rigid_body
  use ode
  use rigid_body
  use integrate, only: rk_step, func
  use progressbar, only: progress_bar, remaining_time, progress_bar_time
  implicit none

contains

  subroutine say_hello
    logical :: print_version = .false.
    real(sp) :: a(9)
    real(sp) :: b(3, 3)
    
    print *, "sp: 1/3 = ", 1.0_sp/3.0_sp
    print *, "dp: 1/3 = ", 1.0_dp/3.0_dp

    if (print_version) then
      print '(4a)', 'This file was compiled by ', &
                   compiler_version(), ' using the options ', &
                   compiler_options()
    end if

    b = eye(3, 3)
    a = pack(b, .true.)
    print *, a
    print *, b(:, 1)
    print *, b(:, 2)
    print *, b(:, 3)

  end subroutine say_hello


  subroutine run_sim(t, y, dt, nsteps_out)
    real(dp), intent(inout) :: t(:), y(:,:)
    real(dp), intent(in) :: dt
    integer, intent(in) :: nsteps_out
    real(dp) :: m
    real(dp), dimension(3) :: x0, v0, w0
    real(dp), dimension(3, 3) :: I0, R0
    integer :: tshape(1)
    integer :: i, nsteps
    type(t_rigid_body) :: rb
    integer :: io, io_w, io_p, stat
    logical :: exists

    m = 1.0_dp
    I0(:, 1) = [0.012_dp, 0.0_dp, 0.0_dp]
    I0(:, 2) = [0.0_dp, 0.018_dp, 0.0_dp]
    I0(:, 3) = [0.0_dp, 0.0_dp, 0.006_dp]
    x0 = [0.0_dp, 0.0_dp, 0.0_dp]
    v0 = [0.0_dp, 0.0_dp, 0.0_dp]
    R0(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
    R0(:, 2) = [0.0_dp, 1.0_dp, 0.0_dp]
    R0(:, 3) = [0.0_dp, 0.0_dp, 1.0_dp]
    w0 = [0.0_dp, 0.0_dp, 0.0_dp]

    tshape = shape(t)
    nsteps = tshape(1) - 1

    call init_rigid_body(rb, m, I0, x0, v0, R0, w0, dt)

    inquire(file="out.dat", exist=exists)
    if (exists) then
      open(newunit=io, file="out.dat", iostat=stat)
      if (stat == 0) close(io, status='delete', iostat=stat)
    end if
    
    open(newunit=io, file="out.dat", status="new", action="write")

    t(1) = 0.0_dp
    y(:, 1) = state_to_y(rb)
    write(io, *) t(1), y(4:12, 1), y(16:18, 1)

    inquire(file="omega.dat", exist=exists)
    if (exists) then
      open(newunit=io_w, file="omega.dat", iostat=stat)
      if (stat == 0) close(io_w, status='delete', iostat=stat)
    end if
    
    open(newunit=io_w, file="omega.dat", status="new", action="write")
    write(io_w, *) t(1), rb%angular_velocity

    inquire(file="phi.dat", exist=exists)
    if (exists) then
      open(newunit=io_p, file="phi.dat", iostat=stat)
      if (stat == 0) close(io_p, status='delete', iostat=stat)
    end if
    
    open(newunit=io_p, file="phi.dat", status="new", action="write")
    write(io_p, *) t(1), rb%angular_displacement

    do i=2, nsteps+1
      call update_force(rb)
      t(i) = t(i-1) + dt
      y(:, i) = rk_step(func, y(:, i-1), rb)
      call y_to_state(y(:, i), rb)
      if (mod(i, nsteps_out) == 0) then
        write(io, *) t(i), y(4:12, i), y(16:18, i)
        write(io_w, *) t(i), rb%angular_velocity
        write(io_p, *) t(i), rb%angular_displacement
      end if
      call progress_bar_time(i, nsteps+1)
    end do

    close(io)
  end subroutine


end module rot_lange
