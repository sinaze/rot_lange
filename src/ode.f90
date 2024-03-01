module ode
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    use types, only: t_rigid_body
    use force, only: forcefield
    use utils, only: star
    implicit none
    
contains

    function state_to_y(rb) result(y)
        ! Write all state variables to single 1D array.
        type(t_rigid_body), intent(in) :: rb
        real(dp), dimension(18) :: y

        y(1:3) = rb%position
        y(4:12) = pack(rb%orientation, .true.)
        y(13:15) = rb%momentum
        y(16:18) = rb%angular_mom
    end function

    subroutine y_to_state(y, rb)
        ! Copy state variables from 1D array to object.
        ! Update velocity, inverse inertia tensor and angular velocity.
        real(dp), dimension(18), intent(in) :: y
        type(t_rigid_body), intent(inout) :: rb

        rb%position = y(1:3)
        rb%orientation = reshape(y(4:12), [3, 3])
        rb%momentum = y(13:15)
        rb%angular_mom = y(16:18)
        rb%velocity = rb%momentum / rb%mass
        rb%inv_inertia = matmul(rb%orientation, matmul(rb%inv_inertia_bf, transpose(rb%orientation)))
        rb%angular_velocity_old = rb%angular_velocity
        rb%angular_velocity = matmul(rb%inv_inertia, rb%angular_mom)
        rb%angular_displacement = rb%angular_displacement + (rb%angular_velocity_old + rb%angular_velocity) * rb%dt / 2.0_dp
    end subroutine
    
    subroutine update_force(rb)
        ! Update force and torque on rigid body.
        type(t_rigid_body), intent(inout) :: rb
        real(dp), dimension(6) :: ff

        ff = forcefield(rb)
        rb%force = ff(1:3)
        rb%torque = ff(4:6)
    end subroutine

    function grad_state_to_y(rb) result(grad)
        ! Return gradient of state variables.
        type(t_rigid_body), intent(in) :: rb
        real(dp), dimension(18) :: grad
        real(dp), dimension(3, 3) :: orientation_grad

        grad(1:3) = rb%velocity
        orientation_grad = matmul(star(rb%angular_velocity), rb%orientation)
        grad(4:12) = pack(orientation_grad, .true.)
        grad(13:15) = rb%force
        grad(16:18) = rb%torque
    end function

end module ode
