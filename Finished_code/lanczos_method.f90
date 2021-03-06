module lanczos_method
  use heisenberg_model
  implicit none

  complex :: states(2**num_atoms, 2**num_atoms)
contains

  function mat_element(bra, obs, ket) result(ex_val)
    real    :: obs(:,:)
    real    :: ex_val
    complex :: ket(:), bra(:)
    ex_val = real(dot_product(bra,matmul(obs,ket)))
  end function mat_element

  function norm(vector) result(n)
    complex :: vector(:)
    real    :: n

    n = real(dot_product(vector,vector))
  end function norm

  function create_states() result(states)
    complex :: states(2**num_atoms, 2**num_atoms)
    integer :: i

    do i = 1, 2**num_atoms
      states(i,i) = 1
    end do
  end function create_states

  function generate_random_state(states) result(r_state)
    complex, intent(in) :: states(:,:)
    complex             :: r_state(size(states,2))
    integer             :: i
    complex             :: c_val(size(states,1))
    real                :: r(2,size(states,1)), tc = 0
    call random_number(r)
    r_state = 0
    do i = 1, size(states,1)
      c_val(i) = cmplx(r(1,i),r(2,i))
      r_state = r_state+c_val(i)*states(i,:)
      tc = norm(r_state)
    end do
    r_state = r_state/tc**0.5
  end function generate_random_state

  subroutine create_krylov_space(m, krylov)
    integer, intent(in)               :: m
    complex, allocatable, intent(out) :: krylov(:,:)
    integer                           :: i
    complex                           :: r_state(2**num_atoms)

    states = create_states()
    r_state = generate_random_state(states)
    call calculate_hamiltonian(Hmlt)
    allocate(krylov(m+1,size(r_state)))
    krylov(1,:) = r_state
    do i = 1, m
      krylov(i+1,:) = matmul(Hmlt,krylov(i,:))
    end do
  end subroutine create_krylov_space

  subroutine create_triHmlt(m, tri_Hmlt)
    integer, intent(in)            :: m
    complex, allocatable           :: lanczos(:,:)
    real, allocatable, intent(out) :: tri_Hmlt(:,:)
    complex, allocatable           :: krylov(:,:)
    complex, allocatable           :: r_state(:)
    integer                        :: i
    real                           :: a(m+1)
    real                           :: N(m+1)

    call create_krylov_space(m, krylov)
    r_state = generate_random_state(krylov)

    allocate(lanczos(m+1,size(r_state)))

    lanczos(1,:) = r_state
    a(1) = mat_element(lanczos(1,:),Hmlt,lanczos(1,:))
    N(1) = norm(lanczos(1,:))

    lanczos(2,:) = matmul(Hmlt, lanczos(1,:)) - a(1)*lanczos(1,:)
    N(2) = norm(lanczos(2,:))
    lanczos(2,:) = lanczos(2,:)/N(2)**0.5
    a(2) = mat_element(lanczos(2,:),Hmlt,lanczos(2,:))

    do i = 2, m
      lanczos(i+1,:) = (matmul(Hmlt,lanczos(i,:)) - a(i)*lanczos(i,:) - N(i)*lanczos(i-1,:))
      N(i+1) = norm(lanczos(i+1,:))
      lanczos(i+1,:) = lanczos(i+1,:)/N(i+1)**0.5
      a(i+1) = mat_element(lanczos(i+1,:),Hmlt,lanczos(i+1,:))
    end do

    !write(*,*) N
    !write(*,*) a
    allocate(tri_Hmlt(m+1,m+1))
    tri_Hmlt = 0
    do i = 1, m+1
      tri_Hmlt(i,i) = a(i)
      if(i<m+1) then
        tri_Hmlt(i+1,i) = N(i+1)**0.5
        tri_Hmlt(i,i+1) = N(i+1)**0.5
      end if
    end do
  end subroutine create_triHmlt

  subroutine calc_eigenval(Hmlt, e_val)
    real, intent(inout)   :: Hmlt(:,:)
    integer, parameter    :: n = 3
    real, intent(out)     :: e_val(:)
    integer, parameter    :: lda = n, ldvl = n, ldvr = n
    integer, parameter    :: lwmax = 1000
    integer               :: lwork, info
    real                  :: vl(ldvl, n), vr(ldvr, n), wr(n), wi(n), work(lwmax)
    lwork = -1
    CALL SGEEV( 'V', 'V', N, Hmlt, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
    LWORK = MIN(LWMAX, INT(WORK(1)))
    CALL SGEEV( 'V', 'V', N, Hmlt, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
    e_val = wr
  end subroutine calc_eigenval
end module lanczos_method
