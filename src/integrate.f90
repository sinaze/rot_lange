module integrate
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    use types, only: t_rigid_body
    use ode
    use rigid_body
    implicit none
    
    abstract interface
        function int_func(y, rb) result(f)
            import :: dp, t_rigid_body
            real(dp), dimension(18), intent(in) :: y
            type(t_rigid_body), intent(in) :: rb
            real(dp), dimension(18) :: f
        end function
    end interface

contains

    function func(y, rb) result(grad)
        real(dp), intent(in) :: y(18)
        type(t_rigid_body), intent(in) :: rb
        type(t_rigid_body) :: rb_dummy
        real(dp) :: m, dt !, t_end
        real(dp), dimension(3) :: x0, v0, w0
        real(dp), dimension(3, 3) :: I0, R0
        real(dp), dimension(18) :: grad

        dt = rb%dt
        m = 1.0_dp
        I0(:, 1) = [0.006_dp, 0.0_dp, 0.0_dp]
        I0(:, 2) = [0.0_dp, 0.018_dp, 0.0_dp]
        I0(:, 3) = [0.0_dp, 0.0_dp, 0.012_dp]
        ! TEST
        ! I0(:, 1) = [60.0_dp, 0.0_dp, 0.0_dp]
        ! I0(:, 2) = [0.0_dp, 180.0_dp, 0.0_dp]
        ! I0(:, 3) = [0.0_dp, 0.0_dp, 120.0_dp]
        ! TEST
        x0 = [0.0_dp, 0.0_dp, 0.0_dp]
        v0 = [0.0_dp, 0.0_dp, 0.0_dp]
        R0(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
        R0(:, 2) = [0.0_dp, 1.0_dp, 0.0_dp]
        R0(:, 3) = [0.0_dp, 0.0_dp, 1.0_dp]
        w0 = [0.0_dp, 0.0_dp, 0.0_dp]

        call init_rigid_body(rb_dummy, m, I0, x0, v0, R0, w0, dt)

        call y_to_state(y, rb_dummy)
        rb_dummy%force = rb%force
        rb_dummy%torque = rb%torque
        grad = grad_state_to_y(rb_dummy)
    end function func

    function rk_step(f, y, rb) result(step)
        real(dp), intent(in) :: y(18)
        type(t_rigid_body), intent(in) :: rb
        real(dp) :: dt
        real(dp) :: step(18)
        procedure(int_func) :: f
        real(dp), dimension(18) :: k1, k2, k3, k4

        dt = rb%dt
        k1 = f(y, rb)
        k2 = f(y + 0.5_dp * dt * k1, rb)
        k3 = f(y + 0.5_dp * dt * k2, rb)
        k4 = f(y + dt * k3, rb)
        step = y + dt / 6.0_dp * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4)
    end function
    
end module integrate
