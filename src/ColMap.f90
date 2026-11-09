module m_colmap

  implicit none
  
  ! Data structure to keep track of where data is
  type :: t_dataset_colmap
    integer  :: nsets, maxcol
    integer, allocatable :: class_ids(:)  ! Global ID number for each class (in order) (nsets)
    integer, allocatable :: col_val(:)    ! Class value columns (nsets)
    integer, allocatable :: col_var(:)    ! Class variance columns (nsets)
    integer, allocatable :: col_hsu(:)    ! HSU columns (nlayers)
  end type
  
  ! Related routines
  contains
  
  subroutine append_colmap(colmaps, cm)
    implicit none
    type(t_dataset_colmap), allocatable, intent(inout) :: colmaps(:)
    type(t_dataset_colmap), intent(in)                 :: cm
    type(t_dataset_colmap), allocatable                :: tmp(:)
    integer :: n

    if (.not. allocated(colmaps)) then
      allocate(colmaps(1))
      colmaps(1) = cm
      return
    end if

    n = size(colmaps)
    allocate(tmp(n+1))
    tmp(1:n) = colmaps
    tmp(n+1) = cm
    call move_alloc(tmp, colmaps)
  end subroutine append_colmap
  
end module m_colmap