module mjd_ghosh
  use heisenberg_model
  implicit none

contains

  subroutine calc_mg_hmlt(mg_Hmlt)
    real, allocatable, intent(out) :: mg_Hmlt(:,:)
    integer                        :: i, j, k

    call create_spin_vector()
    allocate(mg_Hmlt, source = real(pauli_spin_final(1)%x_mat))
    mg_Hmlt = 0

    do i = 0, num_atoms-1
      j = mod(i+1, num_atoms)
      k = mod(i+2, num_atoms)
      mg_Hmlt = mg_Hmlt + dd_matmul(spin_vector(i+1), spin_vector(j+1)) &
                        + 0.5*dd_matmul(spin_vector(i+1), spin_vector(k+1))
    end do
  end subroutine calc_mg_hmlt

end module mjd_ghosh
