module types
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

type :: t_rigid_body
    real(dp) :: mass
    real(dp), dimension(3, 3) :: inertia_bf
    real(dp), dimension(3, 3) :: inv_inertia_bf
    real(dp), dimension(3) :: position
    real(dp), dimension(3) :: velocity
    real(dp), dimension(3) :: momentum
    real(dp), dimension(3, 3) :: orientation
    real(dp), dimension(3) :: angular_velocity
    real(dp), dimension(3) :: angular_velocity_old
    real(dp), dimension(3) :: angular_mom
    real(dp), dimension(3) :: angular_displacement
    real(dp), dimension(3, 3) :: inv_inertia
    real(dp), dimension(3) :: force
    real(dp), dimension(3) :: torque
    real(dp) :: dt
end type

type :: t_params
    real(dp) :: dt
    real(dp) :: t_end
    integer :: nsteps
    integer :: nsteps_out
end type

end module types
