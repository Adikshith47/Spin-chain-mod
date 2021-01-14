program heiseprogram_test
  use heisenberg_model
  use lanczos_method
  implicit none
  !complex, allocatable :: r_state(:)
  !complex              :: a
  !integer , parameter            :: m = 2
  real, allocatable              :: tri_Hmlt(:,:)
  real                           :: eval(3)
  call create_triHmlt(2, tri_Hmlt)
  write(*,*) tri_Hmlt
  call calc_eigenval(tri_Hmlt, eval)
  write(*,*) eval
end program heiseprogram_test
