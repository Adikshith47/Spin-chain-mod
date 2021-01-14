program test_site
  use base_mod
  use arb_spin
  use spin_chain_models
  implicit none

  type(Hamiltonian_class) :: test_Hamiltonian
  procedure(define_Hamiltonian), pointer :: Ham_ptr

  call test_Hamiltonian%init_params(0.5,7,1)
  !Ham_ptr => test_Ham1

  call test_Hamiltonian%calculate_Hamiltonian() !call test_Hamiltonian%calculate_Hamiltonian(Ham_ptr)
  call test_Hamiltonian%calculate_eigenval() !test_Hamiltonian%lanczos_eigenval(2)

  !write(*,*) test_Hamiltonian%Hmlt
  write(*,*) test_Hamiltonian%eval

contains
  pure function test_Ham1(H_class) result(custom_Hmlt)
    type(Hamiltonian_class), intent(in) :: H_class
    real, allocatable                    :: custom_Hmlt(:,:)
    allocate(custom_Hmlt(2**(H_class%num_sites),2**(H_class%num_sites)))
    custom_Hmlt = 1
  end function test_Ham1
end program test_site

!call test_Hamiltonian%calculate_Hamiltonian(Ham_ptr)
!call test_Hamiltonian%calculate_eigenval()

!write(*,*) test_Hamiltonian%Hmlt
!write(*,*) test_Hamiltonian%eval

!call test_Hamiltonian%lanczos_eigenval(2)
!write(*,*) test_Hamiltonian%eval
