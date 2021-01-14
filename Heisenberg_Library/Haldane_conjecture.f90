program Haldane_conjecture
  use base_mod
  use arb_spin
  use spin_chain_models
  implicit none

  type(Hamiltonian_class) :: half_ham, int_ham
  call half_ham%init_params(0.5,5,1)
  call int_ham%init_params(1.0,5,1)

  call half_ham%calculate_Hamiltonian()
  call int_ham%calculate_Hamiltonian()

  call half_ham%calculate_eigenval()
  call int_ham%calculate_eigenval()

  write(*,*) half_ham%eval
  write(*,*) int_ham%eval
end program Haldane_conjecture
