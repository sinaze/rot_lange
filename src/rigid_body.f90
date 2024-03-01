module rigid_body
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use types, only: t_rigid_body
    use utils, only: matinv3
    implicit none

contains

subroutine init_rigid_body(rb, m, I0, x0, v0, R0, w0, dt)
    type(t_rigid_body), intent(inout) :: rb
    real(dp), intent(in) :: m
    real(dp), dimension(3, 3), intent(in) :: I0
    real(dp), dimension(3), intent(in) :: x0
    real(dp), dimension(3), intent(in) :: v0
    real(dp), dimension(3, 3), intent(in) :: R0
    real(dp), dimension(3), intent(in) :: w0
    real(dp), intent(in) :: dt

    rb%mass = m
    rb%inertia_bf = I0
    rb%inv_inertia_bf = matinv3(I0)
    rb%position = x0
    rb%velocity = v0
    rb%momentum = m * v0
    rb%orientation = R0
    rb%angular_velocity = w0
    rb%angular_velocity_old = [0.0_dp, 0.0_dp, 0.0_dp]
    rb%angular_mom = matmul(I0, w0)
    rb%angular_displacement = [0.0_dp, 0.0_dp, 0.0_dp]
    rb%inv_inertia = matmul(R0, matmul(rb%inv_inertia_bf, transpose(R0)))
    rb%force = [0.0_dp, 0.0_dp, 0.0_dp]
    rb%torque = [0.0_dp, 0.0_dp, 0.0_dp]
    rb%dt = dt
end subroutine init_rigid_body

end module rigid_body
