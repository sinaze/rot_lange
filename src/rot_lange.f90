module rot_lange
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  implicit none
  private

  public :: say_hello
contains

  subroutine say_hello
    print *, "sp: 1/3 = ", 1.0_sp/3.0_sp
    print *, "dp: 1/3 = ", 1.0_dp/3.0_dp
  end subroutine say_hello

end module rot_lange
