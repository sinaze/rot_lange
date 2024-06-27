module force
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    use types, only: t_rigid_body
    use stdlib_random, only: random_seed
    use stdlib_stats_distribution_normal, only: rvs_normal
    use stdlib_linalg, only: cross_product
    implicit none
    
contains

    function forcefield(rb) result(f)
        type(t_rigid_body), intent(in) :: rb
        real(dp), dimension(6) :: f
        real(dp), dimension(3) :: gamma, force, torque
        real(dp), dimension(3) :: omega_i
        real(dp) :: k_B, temp
        
        omega_i = cross_product(rb%angular_velocity, rb%angular_mom)

        ! gamma = [10.15_dp, 18.84_dp, 7.64_dp]

        ! TEST
        ! gamma = [0.1015_dp, 0.1884_dp, 0.0764_dp]
        gamma = [0.0_dp, 0.0_dp, 0.0_dp]
        ! TEST

        k_B = 8.31446261815342e-3_dp
        temp = 300.0_dp
        force = [0.0_dp, 0.0_dp, 0.0_dp]
        ! torque = -gamma * rb%angular_velocity + sqrt(2.0_dp * k_B * temp * gamma / rb%dt) * rvs_normal(0.0_dp, 1.0_dp, 3)
        
        ! TEST
        torque = -gamma * rb%angular_velocity + sqrt(2.0_dp * k_B * temp * 0.1_dp / rb%dt) * rvs_normal(0.0_dp, 1.0_dp, 3)
        ! torque = -gamma * rb%angular_velocity - omega_i + sqrt(2.0_dp * k_B * temp * gamma / rb%dt) * rvs_normal(0.0_dp, 1.0_dp, 3)
        ! TEST

        f(1:3) = force
        f(4:6) = torque
    end function
    
end module force
