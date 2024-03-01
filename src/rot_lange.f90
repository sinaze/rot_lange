module rot_lange
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options
  use stdlib_linalg, only: eye
  use types, only: t_rigid_body !, t_params
  use ode
  use rigid_body
  use integrate, only: rk_step, func
  implicit none
  ! private

  public :: say_hello
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

  ! subroutine run_sim(params)
  subroutine run_sim(t, y)
    ! type(t_params), intent(in) :: params
    real(dp), intent(inout) :: t(:), y(:,:)
    real(dp) :: m, dt !, t_end
    real(dp), dimension(3) :: x0, v0, w0
    real(dp), dimension(3, 3) :: I0, R0
    integer :: i, nsteps !, nsteps_out
    type(t_rigid_body) :: rb

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

    ! dt = params%dt
    ! t_end = params%t_end
    ! nsteps = params%nsteps
    ! nsteps_out = params%nsteps_out

    dt = 0.001_dp
    ! nsteps = 10000

    ! allocate(t(nsteps+1))
    ! allocate(y(18, nsteps+1))

    call init_rigid_body(rb, m, I0, x0, v0, R0, w0, dt)

    t(1) = 0.0_dp
    y(:, 1) = state_to_y(rb)

    do i=2, nsteps+1
      call update_force(rb)
      t(i) = t(i-1) + dt
      y(:, i) = rk_step(func, y(:, i-1), dt)
      call y_to_state(y(:, i), rb)
    end do
  end subroutine


end module rot_lange
