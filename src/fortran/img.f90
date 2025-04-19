module img
        use iso_c_binding
        implicit none

        interface
                function img_read(fname, w, h) bind(C, name="iio_read_float")
                        import :: c_char, c_int, c_float, c_ptr
                        character(kind=c_char), dimension(*) :: fname
                        integer(c_int), intent(out) :: w, h
                        type(c_ptr) :: read_image
                end function img_read

                subroutine free_c(ptr) bind(C, name="free")
                        import :: c_ptr
                        type(c_ptr), value :: ptr
                end subroutine free_c
         end interface

contains

        subroutine make_c_string(fstr, cstr)
                character(len=*), intent(in) :: fstr
                character(kind=c_char), dimension(:), allocatable, intent(out) :: cstr
                integer :: lenf

                lenf = len_trim(fstr)
                allocate(cstr(0:lenf))
                cstr(0:lenf-1) = transfer(fstr(1:lenf), cstr(0:lenf-1))
                cstr(lenf) = c_null_char
        end subroutine make_c_string

end module img

