module m_cellblocks
  implicit none  
!-------------------------------------------------------------------------------------------------!
! Block Kriging Container
  
  type t_cellblocks
    real,    allocatable :: block_size(:,:)       ! (dim, nodes), native libkriging GQ
  end type

!-------------------------------------------------------------------------------------------------!
  
  contains
  
!-------------------------------------------------------------------------------------------------!
  
  subroutine alloc_block_sizes(blck, nnode, dim)
    implicit none
    type(t_cellblocks), intent(inout) :: blck
    integer, intent(in)               :: nnode, dim

    if (allocated(blck%block_size)) deallocate(blck%block_size)
    allocate(blck%block_size(dim, nnode))
  end subroutine alloc_block_sizes
  
!-------------------------------------------------------------------------------------------------!
  
end module m_cellblocks
