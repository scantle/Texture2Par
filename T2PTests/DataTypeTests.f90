module DataTypeTests
  use assert
  use list_module
  implicit none
  
  contains
  !-----------------------------------------------------------------------------------------------!
  ! Main data type testing routine
  !-----------------------------------------------------------------------------------------------!
  subroutine test_data_types(errors)
    implicit none
    integer, intent(inout)  :: errors
    
    call test_stack_type(errors)
    
  end subroutine test_data_types
  !-----------------------------------------------------------------------------------------------!
  
  !-----------------------------------------------------------------------------------------------!
  ! Subroutines called by main routine
  !-----------------------------------------------------------------------------------------------!
  
  !-----------------------------------------------------------------------------------------------!
  ! Test Stack 
  subroutine test_stack_type(errors)
    implicit none
    integer, intent(inout)  :: errors
    integer                 :: int_array(4) = (/ 10, 12, 42, 98 /), kernel
    type(stack)             :: teststack
    type(list_node)         :: corn
    
    ! Test empty
    call assert_true(teststack%size == 0, "Error in empty stack type size from list module", errors)
    
    ! Test adding values
    call teststack%push(int_array(1))
    call assert_true(teststack%size == 1, "Error in incremented stack type size from list module", errors)
    
    call teststack%push(int_array(2))
    call assert_true(teststack%size == 2, "Error in incremented stack type size from list module", errors)
    
    call teststack%push(int_array(3))
    call assert_true(teststack%size == 3, "Error in incremented stack type size from list module", errors)
    
    call teststack%push(int_array(4))
    call assert_true(teststack%size == 4, "Error in incremented stack type size from list module", errors)
    
    ! Test removing values
    corn = teststack%pop()  ! pop = corn, get it?
    call assert_true(corn%value==int_array(4), "Error in stack pop 1st return value from list module", errors)
    call assert_true(teststack%size == 3, "Error in reduced stack type size from list module", errors)
    
    corn = teststack%pop()
    call assert_true(corn%value==int_array(3), "Error in stack pop 2nd return value from list module", errors)
    call assert_true(teststack%size == 2, "Error in reduced stack type size from list module", errors)
    
    corn = teststack%pop()
    call assert_true(corn%value==int_array(2), "Error in stack pop 3rd return value from list module", errors)
    call assert_true(teststack%size == 1, "Error in reduced stack type size from list module", errors)
    
    corn = teststack%pop()
    call assert_true(corn%value==int_array(1), "Error in stack pop 4th return value from list module", errors)
    call assert_true(teststack%size == 0, "Error in emptied stack type size from list module", errors)
  
  end subroutine test_stack_type
  !-----------------------------------------------------------------------------------------------!
  
end module DataTypeTests