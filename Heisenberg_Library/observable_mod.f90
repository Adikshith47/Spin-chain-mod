module observable_mod
  use spin_chain_models
  use base_mod
  use arb_spin
  implicit none

  type observables
    real, allocatable :: spin_corr(:)
  contains
    procedure :: calculate_spin_corr
    procedure :: calculate_ent_entropy
  end type observables

contains
  subroutine calculate_spin_corr(this, H_class, pos)
    class(observables)        :: this
    type(Hamiltonian_class) :: H_class
    complex, allocatable    :: spin_corr_op(:,:)
    integer                 :: i
    integer, intent(in)                 :: pos
    !spin_corr_op = H_class%spin_vector(1)%vector_matrix
    allocate(spin_corr_op, source=H_class%spin_vector(pos+1)%vector_matrix(1,:,:))
    allocate(this%spin_corr(H_class%num_sites))
    spin_corr_op = 0

    do i = 1, H_class%num_sites
      spin_corr_op = (H_class%spin_vector(pos+1).mm.H_class%spin_vector(i))
      this%spin_corr(i) = real(spin_corr_op(1,1))
    end do
  end subroutine calculate_spin_corr

  subroutine calculate_ent_entropy(this)
    class(observables) :: this
  end subroutine calculate_ent_entropy
end module observable_mod
