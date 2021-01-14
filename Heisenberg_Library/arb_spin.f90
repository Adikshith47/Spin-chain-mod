module arb_spin
  implicit none

  !contains each of the spin matrices, the value of spin and a 2D matrix of all 2^n spin states
  !and the procedures to make them
  type spin_base_matrix
    !real :: spin_val
    complex, allocatable, dimension(:,:) :: x, y, z, raise, lower, i_mat
    real :: spin_value
    real, allocatable :: l_values(:)
    complex, allocatable :: spin_states(:,:)
  contains
    procedure, private :: create_lval
    procedure, private :: create_states
    procedure, private :: ladder_op
    procedure, public  :: create_spin_matrices
  end type spin_base_matrix

contains

  !acts as the raising and lowering opearator
  pure function ladder_op(this, state, op) result(new_state)
    class(spin_base_matrix), intent(in) :: this
    complex, intent(in) :: state(:)
    character(*), intent(in) :: op
    complex :: new_state(size(state))
    integer :: max_pos

    max_pos = maxloc(real(state),1)
    new_state = 0

    if(op=="raise") then

      if(max_pos==ubound(state,1)) then
        new_state = 0
      else
        new_state = ((this%spin_value - this%l_values(max_pos))*(this%spin_value + this%l_values(max_pos) + 1))**0.5&
                  *swap_pos(state,max_pos,max_pos+1)
      end if
    elseif(op=="lower") then

      if(max_pos==lbound(state,1)) then
        new_state = 0
      else
        new_state = ((this%spin_value + this%l_values(max_pos))*(this%spin_value - this%l_values(max_pos) + 1))**0.5&
                  *swap_pos(state,max_pos,max_pos-1)
      end if
    end if
  end function ladder_op

  !small function used for the ladder operator
  pure function swap_pos(state, p1, p2) result(swap_state)
    complex, intent(in) :: state(:)
    integer, intent(in) :: p1, p2
    complex :: swap_state(size(state))
    complex :: temp

    swap_state = state
    temp = swap_state(p1)
    swap_state(p1) = swap_state(p2)
    swap_state(p2) = temp
  end function swap_pos

  !creates the 2^n spin states
  subroutine create_states(this, d)
    class(spin_base_matrix)   :: this
    integer              :: d
    integer              :: i

    allocate(this%spin_states(d,d))
    this%spin_states = 0
    do i = 1, d
      this%spin_states(i,i) = 1
    end do
  end subroutine create_states

  !values from -l to l
  subroutine create_lval(this)
    class(spin_base_matrix) :: this
    real              :: d_spin
    integer           :: i

    d_spin = this%spin_value
    allocate(this%l_values(int(2*this%spin_value)+1))

    do i = size(this%l_values), 1, -1
      this%l_values(i) = d_spin
      d_spin = d_spin-1
    end do
  end subroutine create_lval

  !main subroutine. Creates all the spin matrices
  subroutine create_spin_matrices(this, spin_val)
    class(spin_base_matrix) :: this
    real, intent(inout) :: spin_val
    integer :: d, i, j
    complex, allocatable :: d1_state(:), d2_state(:)
    d = int(2*spin_val) + 1
    this%spin_value = (real(d)-1)/2
    spin_val = this%spin_value
    allocate(this%x(d,d),this%y(d,d),this%z(d,d),this%raise(d,d),&
              this%lower(d,d), this%i_mat(d,d))

    call this%create_states(d)
    call this%create_lval()

    this%z = 0
    this%raise = 0
    this%lower = 0
    this%i_mat = 0

    do concurrent(i=1:d)

      this%z(i,i) = this%l_values(i)
      this%i_mat(i,i) = 1.0
      d1_state = this%ladder_op(this%spin_states(i,:), "raise")
      d2_state = this%ladder_op(this%spin_states(i,:), "lower")
      do concurrent(j=1:d)

        this%raise(j,i) = dot_product(this%spin_states(j,:),d1_state)
        this%lower(j,i) = dot_product(this%spin_states(j,:),d2_state)
      end do
    end do

    this%x = (this%raise + this%lower)/2
    this%y = (this%raise - this%lower)/(2*(0.0,1.0))
  end subroutine create_spin_matrices
end module arb_spin
