module spin_chain_models
  use base_mod
  use arb_spin
  use lanczos_mod
  implicit none

  !base derived data type for hamiltonian. contains the procedure to calculate
  !the matrix and eigenvalues.
  type :: Hamiltonian_class
    real    :: spin_value
    integer :: num_sites
    integer :: interaction_order
    real, allocatable, dimension(:,:) :: Hmlt, e_vec
    real, allocatable                 :: eval(:)
    real, allocatable                 :: rsoi(:)
    type(spin_base_matrix) :: spin_matrix
    type(final_spin), allocatable  :: site_spin_matrix(:)
    type(spin_vec), allocatable :: spin_vector(:)
    type(lanczos_base) :: lanczos
  contains
    procedure :: calculate_eigenval
    procedure :: init_params
    procedure :: calculate_Hamiltonian
    !procedure :: lanczos_eigenval
  end type Hamiltonian_class

  !function that allows to define a custom Hamiltonian, a pointer must be
  !created using this template and then point to the function required.
  abstract interface
    pure function define_Hamiltonian(H_class) result(custom_Hmlt)
      import Hamiltonian_class
      type(Hamiltonian_class), intent(in) :: H_class
      real, allocatable                    :: custom_Hmlt(:,:)
    end function define_Hamiltonian
  end interface

contains

  !Intialize the parameters for a particular system
  subroutine init_params(this, spin_val, sites, order)
    class(Hamiltonian_class) :: this
    integer                 :: i
    real                    :: spin_val
    integer                 :: sites, order

    this%spin_value = spin_val
    this%num_sites = sites
    this%interaction_order = order

    if(allocated(this%site_spin_matrix)) deallocate(this%site_spin_matrix)
    allocate(this%site_spin_matrix(this%num_sites))

    if(allocated(this%spin_vector)) deallocate(this%spin_vector)
    allocate(this%spin_vector(this%num_sites))

    do i = 1, this%num_sites
      this%site_spin_matrix(i)%num_atoms = this%num_sites
    end do

    if(allocated(this%Hmlt)) deallocate(this%Hmlt)
    allocate(this%Hmlt(2**(this%num_sites), 2**(this%num_sites)))
    this%Hmlt = 0

    if(allocated(this%rsoi)) deallocate(this%rsoi)
    allocate(this%rsoi(this%interaction_order))

    do i = 1, this%interaction_order
      this%rsoi(i) = i/2**(i-1)
    end do

    call this%spin_matrix%create_spin_matrices(this%spin_value)
  end subroutine init_params

  !Calculated the Hamiltonian. Possible to define a custom Hamiltonian, otherwise
  !calculates using default nearest neighbour interaction.
  subroutine calculate_Hamiltonian(this, def_ptr)
    class(Hamiltonian_class) :: this
    integer                  :: i=1, j, k
    !procedure(define_Hamiltonian), intent(in) :: def_Hmlt
    procedure(define_Hamiltonian), pointer, intent(in), optional  :: def_ptr

    do i = 1, this%num_sites
      call this%site_spin_matrix(i)%spin_at_position(i, this%spin_matrix)
      call this%spin_vector(i)%create_site_spin_vector(this%site_spin_matrix(i))
    end do

    if(present(def_ptr)) then
      this%Hmlt = def_ptr(this)
    else

      do i = 1, this%interaction_order
        do j = 0, this%num_sites-1
          k = mod(j+this%interaction_order, this%num_sites)
          this%Hmlt = this%Hmlt + this%rsoi(i)*(this%spin_vector(j+1).mm.this%spin_vector(k+1))
        end do
      end do
    end if
  end subroutine calculate_Hamiltonian

  !calculated the eigenvalue
  subroutine calculate_eigenval(this)
    class(Hamiltonian_class) :: this
    integer             :: n
    real, allocatable   :: d_Hmlt(:,:)
    integer             :: lda, ldvl, ldvr
    integer, parameter  :: lwmax = 1000
    real                :: work(lwmax)
    integer             :: lwork, info
    real, allocatable   :: vl(:,:), vr(:,:), wr(:), wi(:)
    if(allocated(this%Hmlt)) then

      d_Hmlt = this%Hmlt

      n = size(this%Hmlt,1)
      lda = n; ldvl = n; ldvr = n
      allocate(vr(ldvr, n), vl(ldvl, n))
      allocate(wr(n), wi(n))

      lwork = -1
      CALL SGEEV( 'V', 'V', N, d_Hmlt, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
      LWORK = MIN(LWMAX, INT(WORK(1)))
      CALL SGEEV( 'V', 'V', N, d_Hmlt, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
      this%eval = wr
      this%e_vec = vl
    else
      write(*,*) "Hamiltonian not created"
    end if
  end subroutine calculate_eigenval

  !subroutine lanczos_eigenval(this, num_eval)
    !class(Hamiltonian_class) :: this
    !integer                  :: num_eval
    !real, allocatable        :: tri_Hmlt(:,:)
    !call this%lanczos%create_triHmlt(this%spin_value, this%num_sites, this%Hmlt,&
    !                                num_eval, tri_Hmlt)
    !this%Hmlt = tri_Hmlt
    !call this%calculate_eigenval()
  !end subroutine lanczos_eigenval
end module spin_chain_models
