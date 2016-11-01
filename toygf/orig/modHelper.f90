module helper
  use parameters
  !use ode

  implicit none
  character(len=50):: file_name


contains

  !**********************************
  !*                                *
  !* Set Parameters of The System   *
  !*                                *
  !**********************************  
  subroutine set_evn()
    implicit none 
    integer:: i = 0
    character(len=30):: buffer, path, evn


    ! set record file path and name
    path = '../Data/'
    file_name = "tra2D"

    file_name = trim(path) // trim(file_name) // '.txt'

    print *, '*****************************************'
    print *, ''

  end subroutine set_evn

  subroutine to_lower(str)
     character(*), intent(in out) :: str
     integer :: i
 
     do i = 1, len(str)
       select case(str(i:i))
         case("A":"Z")
           str(i:i) = achar(iachar(str(i:i))+32)
       end select
     end do
   end subroutine to_lower


  !=======================!
  !=                     =!
  != Recording Functions =!
  !=                     =!
  !=======================!
  subroutine write_file(te)
    implicit none   
    integer:: i
    integer:: mid = 20
    double precision, intent(in):: te

    open(unit = mid, file = file_name, position = "append", action = "write")

    !do i = 1, SNUM
        !write(unit = mid, fmt = *) te
        write(unit = mid, fmt = *) p
    !end do

    close(mid)


  end subroutine write_file



end module helper


