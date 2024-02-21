module rot_lange
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, rot_lange!"
  end subroutine say_hello
end module rot_lange
