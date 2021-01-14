module base_mod
  use arb_spin
  implicit none

  !contains the final spin matrices after calculating the KP
  type final_spin
    complex, allocatable :: x(:,:)
    complex, allocatable :: y(:,:)
    complex, allocatable :: z(:,:)
    integer              :: num_atoms
    real                 :: J_const(3) = [-1.0, -1.0, -1.0]
    integer              :: position
  contains
    procedure :: spin_at_position
  end type final_spin

  !stores the spin matrices in vector form as a 3D matrix
  type spin_vec
    complex, allocatable :: vector_matrix(:,:,:)
  contains
    procedure :: create_site_spin_vector
  end type spin_vec

  !operator for Kronecker product.
  interface operator(.kp.)
    module procedure kronecker_product
  end interface

  !operator for multiplying spin matrices
  interface operator(.mm.)
    module procedure dd_matmul
  end interface
contains

  !KP function
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

 !takes in derived data type spin_vec to multiply the spin matrices and give the
 !final result
  pure function dd_matmul(ddmat1, ddmat2) result(r_mat)
    type(spin_vec), intent(in)  :: ddmat1
    type(spin_vec), intent(in)  :: ddmat2
    real, allocatable           :: r_mat(:,:)
    r_mat = real(matmul(ddmat1%vector_matrix(1,:,:),ddmat2%vector_matrix(1,:,:))) + &
            real(matmul(ddmat1%vector_matrix(2,:,:),ddmat2%vector_matrix(2,:,:))) + &
            real(matmul(ddmat1%vector_matrix(3,:,:),ddmat2%vector_matrix(3,:,:)))
  end function dd_matmul

  !main subroutine that calculate the spin matrices at a site after taking all the kronecker products
  subroutine spin_at_position(this, pos, s_mat)
    class(final_spin)      :: this
    integer                :: pos
    type(spin_base_matrix) :: s_mat
    integer                ::   j
    complex, allocatable   :: temp(:,:), pre_prod(:,:), post_prod(:,:)

    !call s_mat%create_spin_matrices(spin_value)
    this%position = pos

    do j = 1, this%num_atoms

      if(j<this%position) then
        if(j==1) then
          if(allocated(pre_prod)) deallocate(pre_prod)
          pre_prod = s_mat%i_mat
        else
          if(allocated(temp)) deallocate(temp)
          temp = (pre_prod).kp.(s_mat%i_mat)
          deallocate(pre_prod)
          pre_prod = temp
        end if

      elseif(j>this%position) then
        if(j==this%position+1) then
          if(allocated(post_prod)) deallocate(post_prod)
          post_prod = s_mat%i_mat
        else
          if(allocated(temp)) deallocate(temp)
          temp = (post_prod).kp.(s_mat%i_mat)
          deallocate(post_prod)
          post_prod = temp
        end if
      end if
    end do

    if(.not.allocated(pre_prod)) then
      this%x = (s_mat%x).kp.(post_prod)
      this%y = (s_mat%y).kp.(post_prod)
      this%z = (s_mat%z).kp.(post_prod)
    elseif(.not.allocated(post_prod)) then
      this%x = (pre_prod).kp.(s_mat%x)
      this%y = (pre_prod).kp.(s_mat%y)
      this%z = (pre_prod).kp.(s_mat%z)
    else
      this%x = ((pre_prod).kp.(s_mat%x)).kp.(post_prod)
      this%y = ((pre_prod).kp.(s_mat%y)).kp.(post_prod)
      this%z = ((pre_prod).kp.(s_mat%z)).kp.(post_prod)
    end if

    if(allocated(pre_prod)) deallocate(pre_prod)
    if(allocated(post_prod)) deallocate(post_prod)
  end subroutine spin_at_position

  !puts tbe spin matrices in vector form
  subroutine create_site_spin_vector(this, tbt)
    class(spin_vec)       :: this
    class(final_spin)     :: tbt

    if(allocated(tbt%x).and.allocated(tbt%y).and.allocated(tbt%z)) then

      allocate(this%vector_matrix(3,size(tbt%x,1),size(tbt%x,2)))
      this%vector_matrix(1,:,:) = (tbt%J_const(1)**0.5)*(tbt%x)
      this%vector_matrix(2,:,:) = (tbt%J_const(2)**0.5)*(tbt%y)
      this%vector_matrix(3,:,:) = (tbt%J_const(3)**0.5)*(tbt%z)
    else

      write(*,*) "All spin matrices are not allocated"
    end if
  end subroutine create_site_spin_vector
end module base_mod
