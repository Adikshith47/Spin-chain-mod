program mjg_test
  use mjd_ghosh
  use heisenberg_model
  implicit none

  real, allocatable :: mg_Hmlt(:,:)
  real              :: eval(2**num_atoms)
  call calc_mg_hmlt(mg_Hmlt)
  call calculate_eigenval(mg_Hmlt, eval)
  write(*,*) eval
end program mjg_test
