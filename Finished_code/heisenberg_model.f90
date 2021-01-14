module heisenberg_model
  implicit none

  integer            :: up_spin(2) = (/1, 0/), down_spin(2) = (/0, 1/)
  integer, parameter :: num_atoms = 7
  real, allocatable :: Hmlt(:,:)


  type :: spin_half_mat
    complex :: x(2,2) = reshape([(0.0,0.0), (1.0,0.0), (1.0,0.0), (0.0,0.0)], [2,2])
    complex :: y(2,2) = reshape((/(0.0,0.0), (0.0,-1.0), (0.0,1.0), (0.0,0.0)/), (/2,2/))
    complex :: z(2,2) = reshape((/(1.0,0.0), (0.0,0.0), (0.0,0.0), (-1.0,0.0)/), (/2,2/))
  end type spin_half_mat

  type :: spin_one_mat
    complex :: x(3,3) = (1.0/2**0.5)*reshape([(0.0,0.0),(1.0,0.0),(0.0,0.0),(1.0,0.0),&
                            (0.0,0.0),(1.0,0.0),(0.0,0.0),(1.0,0.0),(0.0,0.0)], [3,3])

    complex :: y(3,3) = (1.0/(2**0.5*(0.0,1.0)))*reshape([(0.0,0.0),(1.0,0.0),(0.0,0.0),(-1.0,0.0),&
                            (0.0,0.0),(1.0,0.0),(0.0,0.0),(-1.0,0.0),(0.0,0.0)], [3,3])

    complex :: Z(3,3) = reshape([(1.0,0.0),(0.0,0.0),(0.0,0.0),(0.0,0.0),(0.0,0.0),(0.0,0.0),&
                              (0.0,0.0),(0.0,0.0),(-1.0,0.0)], [3,3])
  end type spin_one_mat

  type :: total_spin_mat
    complex, allocatable :: x_mat(:,:)
    complex, allocatable :: y_mat(:,:)
    complex, allocatable :: z_mat(:,:)
  end type total_spin_mat

  type spin_vec
    complex, allocatable :: mat_vector(:,:,:)
  end type spin_vec

  real :: J_const(3) = [-1.0, -1.0, -1.0]

  type(spin_half_mat)  :: pauli_spin
  type(spin_one_mat)   :: one_spin
  type(total_spin_mat) :: pauli_spin_final(num_atoms)
  type(spin_vec)       :: spin_vector(num_atoms)
  complex              :: i_mat(2,2) = reshape([(1.0,0.0),(0.0,0.0), (0.0,0.0), (1.0,0.0)], [2,2])

  interface operator(.kp.)
    module procedure kronecker_product
  end interface

contains
  pure function kronecker_product(a,b) result(ab)
    complex, intent(in), dimension(:,:) :: a, b
    complex, allocatable                :: ab(:,:)
    integer :: r, ra, rb, c, ca, cb, i, j

    ra = ubound(a, dim=1)
    ca = ubound(a, dim=2)
    rb = ubound(b, dim=1)
    cb = ubound(b, dim=2)

    if(allocated(ab)) deallocate(ab)
    allocate(ab(ra*rb, ca*cb))
    r = 0
    do concurrent(i=1:ra)
      c = 0
      do concurrent(j=1:ca)
        ab(r+1:r + rb, c+1:c+cb) = A(i,j)*b
        c = c+cb
      end do
      r = r+rb
    end do
  end function kronecker_product

  function dd_matmul(ddmat1, ddmat2) result(r_mat)
    type(spin_vec), intent(in)  :: ddmat1
    type(spin_vec), intent(in)  :: ddmat2
    real, allocatable           :: r_mat(:,:)
    r_mat = real(matmul(ddmat1%mat_vector(1,:,:),ddmat2%mat_vector(1,:,:))) + &
            real(matmul(ddmat1%mat_vector(2,:,:),ddmat2%mat_vector(2,:,:))) + &
            real(matmul(ddmat1%mat_vector(3,:,:),ddmat2%mat_vector(3,:,:)))
  end function dd_matmul

  subroutine create_hamiltonian()
    integer                           :: i, j, counter
    complex, allocatable              :: temp(:,:)
    if(allocated(Hmlt)) deallocate(Hmlt)
    do i = 1, num_atoms
      if(allocated(pauli_spin_final(i)%x_mat)) deallocate(pauli_spin_final(i)%x_mat)
      pauli_spin_final(i)%x_mat = i_mat

      if(allocated(pauli_spin_final(i)%y_mat)) deallocate(pauli_spin_final(i)%y_mat)
      pauli_spin_final(i)%y_mat = i_mat

      if(allocated(pauli_spin_final(i)%z_mat)) deallocate(pauli_spin_final(i)%z_mat)
      pauli_spin_final(i)%z_mat = i_mat
    end do

    do concurrent(i=1:num_atoms)
      counter = 0
      do concurrent(j=1:num_atoms)
        if(i==1) then
          if(j==i) then

            pauli_spin_final(i)%x_mat = pauli_spin%x
            pauli_spin_final(i)%y_mat = pauli_spin%y
            pauli_spin_final(i)%z_mat = pauli_spin%z
          else

            if(allocated(temp)) deallocate(temp)
            temp = pauli_spin_final(i)%x_mat.kp.i_mat
            deallocate(pauli_spin_final(i)%x_mat)
            pauli_spin_final(i)%x_mat = temp

            if(allocated(temp)) deallocate(temp)
            temp = pauli_spin_final(i)%y_mat.kp.i_mat
            deallocate(pauli_spin_final(i)%y_mat)
            pauli_spin_final(i)%y_mat = temp

            if(allocated(temp)) deallocate(temp)
            temp = pauli_spin_final(i)%z_mat.kp.i_mat
            deallocate(pauli_spin_final(i)%z_mat)
            pauli_spin_final(i)%z_mat = temp
          end if

        else

          if(j+1==i) then

            if(allocated(temp)) deallocate(temp)
            temp = pauli_spin_final(i)%x_mat.kp.pauli_spin%x
            deallocate(pauli_spin_final(i)%x_mat)
            pauli_spin_final(i)%x_mat = temp

            if(allocated(temp)) deallocate(temp)
            temp = pauli_spin_final(i)%y_mat.kp.pauli_spin%y
            deallocate(pauli_spin_final(i)%y_mat)
            pauli_spin_final(i)%y_mat = temp

            if(allocated(temp)) deallocate(temp)
            temp = pauli_spin_final(i)%z_mat.kp.pauli_spin%z
            deallocate(pauli_spin_final(i)%z_mat)
            pauli_spin_final(i)%z_mat = temp

          else

            counter = counter+1
            if(counter<=num_atoms-2) then

              if(allocated(temp)) deallocate(temp)
              temp = pauli_spin_final(i)%x_mat.kp.i_mat
              deallocate(pauli_spin_final(i)%x_mat)
              pauli_spin_final(i)%x_mat = temp

              if(allocated(temp)) deallocate(temp)
              temp = pauli_spin_final(i)%y_mat.kp.i_mat
              deallocate(pauli_spin_final(i)%y_mat)
              pauli_spin_final(i)%y_mat = temp

              if(allocated(temp)) deallocate(temp)
              temp = pauli_spin_final(i)%z_mat.kp.i_mat
              deallocate(pauli_spin_final(i)%z_mat)
              pauli_spin_final(i)%z_mat = temp

            end if
          end if
        end if
      end do
    end do
  end subroutine create_hamiltonian

  subroutine create_spin_vector()
    integer :: i

    call create_hamiltonian()

    do i = 0, num_atoms-1
      if(allocated(spin_vector(i+1)%mat_vector)) deallocate(spin_vector(i+1)%mat_vector)
      allocate(spin_vector(i+1)%mat_vector(3, size(pauli_spin_final(i+1)%x_mat,1), &
                                          size(pauli_spin_final(i+1)%x_mat,2)))

     spin_vector(i+1)%mat_vector(1,:,:) = J_const(1)*pauli_spin_final(i+1)%x_mat
     spin_vector(i+1)%mat_vector(2,:,:) = J_const(2)*pauli_spin_final(i+1)%y_mat
     spin_vector(i+1)%mat_vector(3,:,:) = J_const(3)*pauli_spin_final(i+1)%z_mat
   end do
  end subroutine create_spin_vector

  subroutine calculate_hamiltonian(Hmlt)
    ! using periodic boundary conditions
    real, allocatable, intent(out)             :: Hmlt(:,:)
    integer                                       :: i, j

    call create_spin_vector()

    allocate(Hmlt, source=real(pauli_spin_final(1)%x_mat))
    Hmlt = 0

    do i = 0, num_atoms-1
      j = mod(i+1, num_atoms)
      Hmlt = Hmlt + dd_matmul(spin_vector(i+1), spin_vector(j+1))
    end do
  end subroutine calculate_hamiltonian

  subroutine calculate_eigenval(Hmlt, e_val)
    real, intent(inout)   :: Hmlt(:,:)
    integer, parameter    :: n = 2**num_atoms
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
  end subroutine calculate_eigenval

end module heisenberg_model
