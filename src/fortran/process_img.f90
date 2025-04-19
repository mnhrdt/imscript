program process_img
	use iso_c_binding
        use img
        implicit none

        character(len=100) :: f90_filename
        character(kind=c_char), dimension(:), allocatable :: c_filename
        integer(c_int) :: w, h
        type(c_ptr) :: pixel_ptr
        real(c_float), pointer :: pixels(:)
        real :: average
        integer :: i

        f90_filename = "gbarbara.npy"
        call make_c_string(trim(f90_filename), c_filename)

        pixel_ptr = img_read(c_filename, w, h)

        if (.not. c_associated(pixel_ptr)) then
                print *, "Failed to read image."
                stop 1
        end if

        call c_f_pointer(pixel_ptr, pixels, [w * h])

        average = 0.0
        do i = 1, w * h
                average = average + pixels(i)
        end do

        average = average / (w * h)

        print *, "Image size: ", w, "x", h
        print *, "Average brightness:", average

        call free_c(pixel_ptr)
end program process_img
