module integrate
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    use types, only: t_rigid_body
    use ode
    use rigid_body
    implicit none
    
    abstract interface
        function int_func(y) result(f)
            import :: dp
            real(dp), dimension(18), intent(in) :: y
            real(dp), dimension(18) :: f
        end function
    end interface

contains

    function func(y) result(grad)
        real(dp), intent(in) :: y(18)
        type(t_rigid_body) :: rb
        real(dp) :: m, dt !, t_end
        real(dp), dimension(3) :: x0, v0, w0
        real(dp), dimension(3, 3) :: I0, R0
        real(dp), dimension(18) :: grad

        dt = 0.001_dp
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

        call init_rigid_body(rb, m, I0, x0, v0, R0, w0, dt)

        call y_to_state(y, rb)
        ! wrong!
        call update_force(rb)
        grad = grad_state_to_y(rb)
    end function func

    function rk_step(f, y, dt) result(step)
        real(dp), intent(in) :: y(18)
        real(dp), intent(in) :: dt
        real(dp) :: step(18)
        procedure(int_func) :: f
        real(dp), dimension(18) :: k1, k2, k3, k4

        k1 = f(y)
        k2 = f(y + 0.5_dp * dt * k1)
        k3 = f(y + 0.5_dp * dt * k2)
        k4 = f(y + dt * k3)
        step = y + dt / 6.0_dp * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4)
    end function

    ! subroutine rk4(f, y, t, omega, phi, rb, dt)
    !     real(dp), intent(inout) :: y(18), t
    !     real(dp), intent(in) :: dt
    !     real(dp), dimension(:,:), intent(inout) :: omega, phi
    !     type(t_rigid_body), intent(in) :: rb
    !     procedure(int_func) :: f
    !     integer :: i
        
    !     omega(:, 1) = rb%angular_velocity
    !     phi(:, 1) = rb%angular_displacement
    !     do i = 1, iend
    !         y(:, i+1) = rk_step(f, y(:, i), dt)
    !         omega(:, i+1) = rb%angular_velocity
    !         phi(:, i+1) = rb%angular_displacement
    !     end do
    ! end subroutine

    
end module integrate
