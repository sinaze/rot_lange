module force
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    use types, only: t_rigid_body
    use stdlib_random, only: random_seed
    use stdlib_stats_distribution_normal, only: rvs_normal
    implicit none
    
contains

    function forcefield(rb) result(f)
        type(t_rigid_body), intent(in) :: rb
        real(dp), dimension(6) :: f
        real(dp), dimension(3) :: gamma, force, torque
        real(dp) :: k_B, temp

        gamma = [7.69_dp, 18.52_dp, 10.0_dp]
        k_B = 8.31446261815342e-3_dp
        temp = 300.0_dp
        force = [0.0_dp, 0.0_dp, 0.0_dp]
        torque = -gamma * rb%angular_velocity + sqrt(2.0_dp * k_B * temp * gamma / rb%dt) * rvs_normal(0.0_dp, 1.0_dp, 3)
        f(1:3) = force
        f(4:6) = torque
    end function
    
end module force
