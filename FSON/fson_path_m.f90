! Copyright (c) 2012 Joseph A. Levin
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of this
! software and associated documentation files (the "Software"), to deal in the Software
! without restriction, including without limitation the rights to use, copy, modify, merge,
! publish, distribute, sublicense, and/or sell copies of the Software, and to permit 
! persons to whom the Software is furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or 
! substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
! INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
! PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
! LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
! OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! DEALINGS IN THE SOFTWARE.

!     
! File:   fson_path_m.f95
! Author: Joseph A. Levin
!
! Created on March 10, 2012, 11:01 PM
!

module fson_path_m
    
    use fson_value_m 
    use fson_string_m

    private
    
    public :: fson_path_get
    
    interface fson_path_get
        module procedure get_by_path
        module procedure get_integer
        module procedure get_real
        module procedure get_double
        module procedure get_logical
        module procedure get_chars
        module procedure get_array_1d_integer
        module procedure get_array_2d_integer
        module procedure get_array_1d_real
        module procedure get_array_2d_real
        module procedure get_array_1d_double
        module procedure get_array_2d_double
        module procedure get_array_1d_logical
        module procedure get_array_2d_logical
        module procedure get_array_1d_char
        module procedure get_array_2d_char
    end interface fson_path_get

    abstract interface

       subroutine array_callback_1d(element, i, count)
         use fson_value_m
         implicit none
         type(fson_value), pointer,intent(in) :: element
         integer, intent(in) :: i        ! index
         integer, intent(in) :: count    ! size of array
       end subroutine array_callback_1d

       subroutine array_callback_2d(element, i1, i2, count1, count2)
         use fson_value_m
         implicit none
         type(fson_value), pointer,intent(in) :: element
         integer, intent(in) :: i1, i2
         integer, intent(in) :: count1, count2
       end subroutine array_callback_2d

    end interface
        integer,allocatable:: arr_int_1d(:),arr_int_2d(:,:)
        real,allocatable:: arr_real_1d(:),arr_real_2d(:,:)
        double precision,allocatable:: arr_double_1d(:),arr_double_2d(:,:)
        logical ,allocatable:: arr_log_1d(:),arr_log_2d(:,:)
        character ,allocatable:: arr_ch_1d(:),arr_ch_2d(:,:)
contains
    !
    ! GET BY PATH
    !
    ! $     = root 
    ! @     = this
    ! .     = child object member
    ! []    = child array element
    !
    recursive subroutine get_by_path(this, path, p)
        type(fson_value), pointer,intent(in) :: this        
        type(fson_value), pointer,intent(inout) ::  p        
        character(len=*),intent(in) :: path
        integer :: i, length, child_i
        character :: c
        logical :: array        
                
        ! default to assuming relative to this
        p => this
        
        child_i = 1          
        
        array = .false.
        
        length = len_trim(path)
        
        do i=1, length
            c = path(i:i)    
            select case (c)
                case ("$")
                    ! root
                    do while (associated (p % parent))
                        p => p % parent
                    end do
                    child_i = i + 1
                case ("@")
                    ! this                    
                    p => this
                    child_i = i + 1
                case (".", "[")                    
                    ! get child member from p                          
                    if (child_i < i) then                          
                        p => fson_value_get(p, path(child_i:i-1))
                    else
                        child_i = i + 1
                        cycle
                    end if
                    
                    if(.not.associated(p)) then
                        return                                        
                    end if
                    
                    child_i = i+1
                    
                    ! check if this is an array
                    ! if so set the array flag
                    if (c == "[") then
                        ! start looking for the array element index
                        array = .true.
                    end if
                case ("]")
                    if (.not.array) then
                        print *, "ERROR: Unexpected ], not missing preceding ["
                        stop(1)
                    end if
                    array = .false.
                    child_i = parse_integer(path(child_i:i-1))                                                
                    p => fson_value_get(p, child_i)                                                                                                                    
                    
                    child_i= i + 1                                     
            end select            
        end do
                
        ! grab the last child if present in the path
        if (child_i <= length) then            
            p => fson_value_get(p, path(child_i:i-1))                    
            if(.not.associated(p)) then
                return
            else                
            end if
        end if
                
        
    end subroutine get_by_path
    
    !
    ! PARSE INTEGER
    !
    integer function parse_integer(chars) result(integral)
        character(len=*) :: chars
        character :: c
        integer :: tmp, i
                
        integral = 0        
        do i=1, len_trim(chars)
            c = chars(i:i)            
            select case(c)
                case ("0":"9")
                    ! digit        
                    read (c, '(i1)') tmp                                               
                    
                    ! shift
                    if(i > 1) then
                        integral = integral * 10
                    end if
                    ! add
                    integral = integral + tmp
                                                    
                case default                          
                    return
            end select            
        end do
    
    end function parse_integer    
    
    !
    ! GET INTEGER
    !
    subroutine get_integer(this, path, value)
        type(fson_value), pointer,intent(in) :: this
        character(len=*),intent(in), optional :: path
        integer,intent(out) :: value        
        type(fson_value), pointer :: p
        
        
        nullify(p)                
        if(present(path)) then
            call get_by_path(this=this, path=path, p=p)
        else
            p => this
        end if
        
        if(.not.associated(p)) then
            print *, "Unable to resolve path: ", path
            stop(1)
        end if
                
        
        if(p % value_type == TYPE_INTEGER) then            
            value = p % value_integer
        else if (p % value_type == TYPE_REAL) then
            value = p % value_real
        else if (p % value_type == TYPE_LOGICAL) then
            if (p % value_logical) then
                value = 1
            else
                value = 0
            end if
        else
            print *, "Unable to resolve value to integer: ", path
            stop(1)
        end if
        
    end subroutine get_integer
    
    !
    ! GET REAL
    !
    subroutine get_real(this, path, value)
        type(fson_value), pointer,intent(in) :: this 
        character(len=*), optional,intent(in) :: path
        real,intent(out) :: value        
        type(fson_value), pointer :: p
        
        
        nullify(p)                
        
        if(present(path)) then
            call get_by_path(this=this, path=path, p=p)
        else
            p => this
        end if
        
        if(.not.associated(p)) then
            print *, "Unable to resolve path: ", path
            stop(1)
        end if
                
        
        if(p % value_type == TYPE_INTEGER) then            
            value = p % value_integer
        else if (p % value_type == TYPE_REAL) then
            value = p % value_real
        else if (p % value_type == TYPE_LOGICAL) then
            if (p % value_logical) then
                value = 1
            else
                value = 0
            end if
        else
            print *, "Unable to resolve value to real: ", path
            stop(1)
        end if
        
    end subroutine get_real
    
    !
    ! GET DOUBLE
    !
    subroutine get_double(this, path, value)
        type(fson_value), pointer, intent(in) :: this
        character(len=*), optional, intent(in) :: path
        double precision, intent(out) :: value        
        type(fson_value), pointer ::  p
        
        
        nullify(p)                
        
        if(present(path)) then
            call get_by_path(this=this, path=path, p=p)
        else
            p => this
        end if
        
        if(.not.associated(p)) then
            print *, "Unable to resolve path: ", path
            stop(1)
        end if
                
        
        if(p % value_type == TYPE_INTEGER) then            
            value = p % value_integer
        else if (p % value_type == TYPE_REAL) then
            value = p % value_double
        else if (p % value_type == TYPE_LOGICAL) then
            if (p % value_logical) then
                value = 1
            else
                value = 0
            end if
        else
            print *, "Unable to resolve value to double: ", path
            stop(1)
        end if
        
    end subroutine get_double
    
    
    !
    ! GET LOGICAL
    !
    subroutine get_logical(this, path, value)
        type(fson_value), pointer, intent(in) :: this
        character(len=*), optional,intent(in) :: path
        logical,intent(out) :: value        
        type(fson_value), pointer ::  p
        
        
        nullify(p)                
        
        if(present(path)) then
            call get_by_path(this=this, path=path, p=p)
        else
            p => this
        end if
        
        if(.not.associated(p)) then
            print *, "Unable to resolve path: ", path
            stop(1)
        end if
                
        
        if(p % value_type == TYPE_INTEGER) then            
            value = (p % value_integer > 0)       
        else if (p % value_type == TYPE_LOGICAL) then
            value = p % value_logical
        else
            print *, "Unable to resolve value to real: ", path
            stop(1)
        end if
        
    end subroutine get_logical
    
    !
    ! GET CHARS
    !
    subroutine get_chars(this, path, value)
        type(fson_value), pointer,intent(in) :: this
        character(len=*), optional,intent(in) :: path
        character(len=*),intent(inout) :: value  
        type(fson_value), pointer ::  p
        nullify(p)                
        
        if(present(path)) then
            call get_by_path(this=this, path=path, p=p)
        else
            p => this
        end if
        
        if(.not.associated(p)) then
            print *, "Unable to resolve path: ", path
            stop(1)
        end if
                
        
        if(p % value_type == TYPE_STRING) then            
            call fson_string_copy(p % value_string, value)          
        else
            print *, "Unable to resolve value to characters: ", path
            stop(1)
        end if
        
    end subroutine get_chars
    
    !
    ! GET ARRAY 1D
    !
    
    subroutine get_array_1d(this, path, array_callback)
        type(fson_value), pointer :: this
        character(len = *), optional :: path
        procedure(array_callback_1d) :: array_callback

        type(fson_value), pointer :: p, element
        integer :: index, count
                
        nullify(p)                
        
        ! resolve the path to the value
        if(present(path)) then
            call get_by_path(this=this, path=path, p=p)
        else
            p => this
        end if
            
        if(.not.associated(p)) then
            print *, "Unable to resolve path: ", path
            stop(1)
        end if
        
        if(p % value_type == TYPE_ARRAY) then            
            count = fson_value_count(p)
            element => p % children
            do index = 1, count
                call array_callback(element, index, count)
                element => element % next
            end do
        else
            print *, "Resolved value is not an array. ", path
            stop(1)
        end if

        if (associated(p)) nullify(p)

      end subroutine get_array_1d

!
! GET ARRAY INTEGER 1D
!
    subroutine get_array_1d_integer(this, path, arr)

      implicit none
      type(fson_value), pointer, intent(in) :: this
      character(len=*), intent(in), optional :: path   
      integer, allocatable, intent(out) :: arr(:)
        integer::l
      if (allocated(arr)) deallocate(arr)
      if (allocated(arr_int_1d)) deallocate(arr_int_1d)
      call get_array_1d(this, path, array_callback_1d_integer)

     if (allocated(arr_int_1d)) then
              l=SIZE(arr_int_1d)
              allocate(arr(l))
              arr=arr_int_1d
              deallocate(arr_int_1d)
      endif

    end subroutine get_array_1d_integer

      subroutine array_callback_1d_integer(element, i, count)
        implicit none
        type(fson_value), pointer, intent(in) :: element
        integer, intent(in) :: i, count
        if (.not. allocated(arr_int_1d)) allocate(arr_int_1d(count))
        call fson_path_get(element, "", arr_int_1d(i))
      end subroutine array_callback_1d_integer
!
! GET ARRAY REAL 1D
!
    subroutine get_array_1d_real(this, path, arr)

      implicit none
      type(fson_value), pointer, intent(in) :: this
      character(len=*), intent(in), optional :: path   
      real, allocatable, intent(out) :: arr(:)

        integer::l
      if (allocated(arr)) deallocate(arr)
      if (allocated(arr_real_1d)) deallocate(arr_real_1d)
      call get_array_1d(this, path, array_callback_1d_real)

     if (allocated(arr_real_1d)) then
              l=SIZE(arr_real_1d)
              allocate(arr(l))
              arr=arr_real_1d
              deallocate(arr_real_1d)
      endif


    end subroutine get_array_1d_real

      subroutine array_callback_1d_real(element, i, count)
        implicit none
        type(fson_value), pointer, intent(in) :: element
        integer, intent(in) :: i, count
        if (.not. allocated(arr_real_1d)) allocate(arr_real_1d(count))
        call fson_path_get(element, "", arr_real_1d(i))
      end subroutine array_callback_1d_real
!
! GET ARRAY DOUBLE 1D
!
    subroutine get_array_1d_double(this, path, arr)

      implicit none
      type(fson_value), pointer, intent(in) :: this
      character(len=*), intent(in), optional :: path   
      double precision, allocatable, intent(out) :: arr(:)
        integer::l
        
      if (allocated(arr)) deallocate(arr)
      if (allocated(arr_double_1d)) deallocate(arr_double_1d)
      call get_array_1d(this, path, array_callback_1d_double)

     if (allocated(arr_double_1d)) then
              l=SIZE(arr_double_1d)
              allocate(arr(l))
              arr=arr_double_1d
              deallocate(arr_double_1d)
      endif
    end subroutine get_array_1d_double
      subroutine array_callback_1d_double(element, i, count)
        implicit none
        type(fson_value), pointer, intent(in) :: element
        integer, intent(in) :: i, count
        if (.not. allocated(arr_double_1d)) allocate(arr_double_1d(count))
        call fson_path_get(element, "", arr_double_1d(i))
      end subroutine array_callback_1d_double

!
! GET ARRAY LOGICAL 1D
!
    subroutine get_array_1d_logical(this, path, arr)

      implicit none
      type(fson_value), pointer, intent(in) :: this
      character(len=*), intent(in), optional :: path   
      logical, allocatable, intent(out) :: arr(:)
        integer::l
      if (allocated(arr)) deallocate(arr)
      if (allocated(arr_log_1d)) deallocate(arr_log_1d)
      call get_array_1d(this, path, array_callback_1d_logical)

     if (allocated(arr_log_1d)) then
              l=SIZE(arr_log_1d)
              allocate(arr(l))
              arr=arr_log_1d
              deallocate(arr_log_1d)
      endif


    end subroutine get_array_1d_logical
      subroutine array_callback_1d_logical(element, i, count)
        implicit none
        type(fson_value), pointer, intent(in) :: element
        integer, intent(in) :: i, count
        if (.not. allocated(arr_log_1d)) allocate(arr_log_1d(count))
        call fson_path_get(element, "", arr_log_1d(i))
      end subroutine array_callback_1d_logical

!
! GET ARRAY CHAR 1D
!
    subroutine get_array_1d_char(this, path, arr)

      implicit none
      type(fson_value), pointer, intent(in) :: this
      character(len=*), intent(in), optional :: path
      character(len = *), allocatable, intent(out) :: arr(:)
        integer::l
      if (allocated(arr)) deallocate(arr)
      if (allocated(arr_ch_1d)) deallocate(arr_ch_1d)
      call get_array_1d(this, path, array_callback_1d_char)

     if (allocated(arr_ch_1d)) then
              l=SIZE(arr_ch_1d)
              allocate(arr(l))
              arr=arr_ch_1d
              deallocate(arr_ch_1d)
      endif


    end subroutine get_array_1d_char

      subroutine array_callback_1d_char(element, i, count)
        implicit none
        type(fson_value), pointer, intent(in) :: element
        integer, intent(in) :: i, count
        if (.not. allocated(arr_ch_1d)) allocate(arr_ch_1d(count))
        call fson_path_get(element, "", arr_ch_1d(i))
      end subroutine array_callback_1d_char

    !
    ! GET ARRAY 2D
    !
    
    subroutine get_array_2d(this, path, array_callback)
        type(fson_value), pointer :: this
        character(len = *), optional :: path
        procedure(array_callback_2d) :: array_callback

        type(fson_value), pointer :: p, element, item
        integer :: i1, i2, count1, count2, c
                
        nullify(p)                
        
        ! resolve the path to the value
        if(present(path)) then
            call get_by_path(this=this, path=path, p=p)
        else
            p => this
        end if
            
        if(.not.associated(p)) then
            print *, "Unable to resolve path: ", path
            stop(1)
        end if
        
        if(p % value_type == TYPE_ARRAY) then            
            count1 = fson_value_count(p)
            element => p % children
            do i1 = 1, count1
               if (element % value_type == TYPE_ARRAY) then
                  c = fson_value_count(element)
                  if (i1 == 1) then
                     count2 = c
                  else if (c /= count2) then
                     print *, "Resolved value has the wrong number of elements. ", &
                          path, "[", i1, "]"
                     stop(1)
                  end if
                  item => element % children
                  do i2 = 1, count2
                     call array_callback(item, i1, i2, count1, count2)
                     item => item % next
                  end do
                  element => element % next
               else
                  print *, "Resolved value is not an array. ", path, "[", i1, "]"
                  stop(1)
               end if
            end do
        else
            print *, "Resolved value is not an array. ", path
            stop(1)
        end if

        if (associated(p)) nullify(p)

      end subroutine get_array_2d

!
! GET ARRAY INTEGER 2D
!
    subroutine get_array_2d_integer(this, path, arr)

      implicit none
      type(fson_value), pointer, intent(in) :: this
      character(len=*), intent(in), optional :: path   
      integer, allocatable, intent(out) :: arr(:, :)
        integer::sh(2)
      if (allocated(arr)) deallocate(arr)
      if (allocated(arr_int_2d)) deallocate(arr_int_2d)
      call get_array_2d(this, path, array_callback_2d_integer)

     if (allocated(arr_int_2d)) then
              sh=SHAPE(arr_int_2d)
              allocate(arr(sh(1),sh(2)))
              arr=arr_int_2d
              deallocate(arr_int_2d)
      endif
    end subroutine get_array_2d_integer
      subroutine array_callback_2d_integer(element, i1, i2, count1, count2)
        implicit none
        type(fson_value), pointer, intent(in) :: element
        integer, intent(in) :: i1, i2, count1, count2
        if (.not. allocated(arr_int_2d)) allocate(arr_int_2d(count1, count2))
        call fson_path_get(element, "", arr_int_2d(i1, i2))
      end subroutine array_callback_2d_integer

!
! GET ARRAY REAL 2D
!
    subroutine get_array_2d_real(this, path, arr)

      implicit none
      type(fson_value), pointer, intent(in) :: this
      character(len=*), intent(in), optional :: path   
      real, allocatable, intent(out) :: arr(:, :)
        integer::sh(2)
      if (allocated(arr)) deallocate(arr)
      if (allocated(arr_real_2d)) deallocate(arr_real_2d)
      call get_array_2d(this, path, array_callback_2d_real)

      if (allocated(arr_real_2d)) then
              sh=SHAPE(arr_real_2d)
              allocate(arr(sh(1),sh(2)))
              arr=arr_real_2d
              deallocate(arr_real_2d)
      endif


    end subroutine get_array_2d_real
      subroutine array_callback_2d_real(element, i1, i2, count1, count2)
        implicit none
        type(fson_value), pointer, intent(in) :: element
        integer, intent(in) :: i1, i2, count1, count2
        if (.not. allocated(arr_real_2d)) allocate(arr_real_2d(count1, count2))
        call fson_path_get(element, "", arr_real_2d(i1, i2))
      end subroutine array_callback_2d_real

!
! GET ARRAY DOUBLE 2D
!
    subroutine get_array_2d_double(this, path, arr)

      implicit none
      type(fson_value), pointer, intent(in) :: this
      character(len=*), intent(in), optional :: path   
      double precision, allocatable, intent(out) :: arr(:, :)
        integer::sh(2)
      if (allocated(arr)) deallocate(arr)
      if (allocated(arr_double_2d)) deallocate(arr_double_2d)
      call get_array_2d(this, path, array_callback_2d_double)

      if (allocated(arr_double_2d)) then
              sh=SHAPE(arr_double_2d)
              allocate(arr(sh(1),sh(2)))
              arr=arr_double_2d
              deallocate(arr_double_2d)
      endif
    end subroutine get_array_2d_double
      subroutine array_callback_2d_double(element, i1, i2, count1, count2)
        implicit none
        type(fson_value), pointer, intent(in) :: element
        integer, intent(in) :: i1, i2, count1, count2
        if (.not. allocated(arr_double_2d)) allocate(arr_double_2d(count1, count2))
        call fson_path_get(element, "", arr_double_2d(i1, i2))
      end subroutine array_callback_2d_double

!
! GET ARRAY LOGICAL 2D
!
    subroutine get_array_2d_logical(this, path, arr)

      implicit none
      type(fson_value), pointer, intent(in) :: this
      character(len=*), intent(in), optional :: path   
      logical, allocatable, intent(out) :: arr(:, :)
        integer::sh(2)
      if (allocated(arr)) deallocate(arr)
      if (allocated(arr_log_2d)) deallocate(arr_log_2d)
      call get_array_2d(this, path, array_callback_2d_logical)
      
      if (allocated(arr_log_2d)) then
              sh=SHAPE(arr_log_2d)
              allocate(arr(sh(1),sh(2)))
              arr=arr_log_2d
              deallocate(arr_log_2d)
      endif

    end subroutine get_array_2d_logical
      subroutine array_callback_2d_logical(element, i1, i2, count1, count2)
        implicit none
        type(fson_value), pointer, intent(in) :: element
        integer, intent(in) :: i1, i2, count1, count2
        if (.not. allocated(arr_log_2d)) allocate(arr_log_2d(count1, count2))
        call fson_path_get(element, "", arr_log_2d(i1, i2))
      end subroutine array_callback_2d_logical

!
! GET ARRAY CHAR 2D
!
    subroutine get_array_2d_char(this, path, arr)

      implicit none
      type(fson_value), pointer, intent(in) :: this
      character(len=*), intent(in), optional :: path
      character(len = *), allocatable, intent(out) :: arr(:, :)
        integer::sh(2)

      if (allocated(arr)) deallocate(arr)
      if (allocated(arr_ch_2d)) deallocate(arr_ch_2d)
      call get_array_2d(this, path, array_callback_2d_char)

      if (allocated(arr_ch_2d)) then
              sh=SHAPE(arr_ch_2d)
              allocate(arr(sh(1),sh(2)))
              arr=arr_ch_2d
              deallocate(arr_ch_2d)
      endif


    end subroutine get_array_2d_char

      subroutine array_callback_2d_char(element, i1, i2, count1, count2)
        implicit none
        type(fson_value), pointer, intent(in) :: element
        integer, intent(in) :: i1, i2, count1, count2
        if (.not. allocated(arr_ch_2d)) allocate(arr_ch_2d(count1, count2))
        call fson_path_get(element, "", arr_ch_2d(i1, i2))
      end subroutine array_callback_2d_char
end module fson_path_m
