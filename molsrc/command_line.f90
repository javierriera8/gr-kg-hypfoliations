module command_line

! use f2kcli
use prec
use parameters_mod
use userparameters

!use molparameters
!use userparameters

use mol
use grid
use output

implicit none

contains

    subroutine read_params(parameters, file_names, file_prefix, conv_level, attach)
    
    type(parameters_type),                            intent(out)  :: parameters
    character(len = 120), dimension(number_of_files), intent(out)  :: file_names
    character(len = 120), 			      intent(out)  :: file_prefix
    integer,                                          intent(inout):: conv_level
    character(len=120),                               intent(inout):: attach

    character(len = 200) :: input_file_name, program_name, conv_test
    integer              :: number_of_arguments
    character(len = 1)   :: conv_level_char

    namelist / param / parameters, file_names, file_prefix


!    conv_level = 0

    call set_default_grid_parameters
    call set_default_output_parameters

    parameters = default_parameters

    call Get_Command_Argument(0, program_name)
    write (*, *) ! made with http://patorjk.com/software/taag/#p=display&f=Graffiti&t=Type%20Something%20
    write (*, *)'     /$$   /$$                      /$$   /$$$$$$$ '
    write (*, *)'    | $$  | $$                    /$$$$  | $$__  $$'
    write (*, *)'    | $$  | $$ /$$   /$$  /$$$$$$|_  $$  | $$  \ $$'
    write (*, *)'    | $$$$$$$$| $$  | $$ /$$__  $$ | $$  | $$  | $$'
    write (*, *)'    | $$__  $$| $$  | $$| $$  \ $$ | $$  | $$  | $$'
    write (*, *)'    | $$  | $$| $$  | $$| $$  | $$ | $$  | $$  | $$'
    write (*, *)'    | $$  | $$|  $$$$$$$| $$$$$$$//$$$$$$| $$$$$$$/'
    write (*, *)'    |__/  |__/ \____  $$| $$____/|______/|_______/ '
    write (*, *)'               /$$  | $$| $$                       '
    write (*, *)'              |  $$$$$$/| $$                       '
    write (*, *)'               \______/ |__/                       '
    write (*, *) "Hello! This program has been called as: ", trim(program_name)

    number_of_arguments = Command_Argument_Count()
    if(number_of_arguments == 0 .or. number_of_arguments > 3) then
       call command_line_help(trim(program_name))
    else
       call Get_Command_Argument(1, input_file_name)
       if(input_file_name == '?' .or. input_file_name == '-?' .or. &
            input_file_name == 'h' .or. input_file_name == '-h' .or. &
            input_file_name == 'help' .or. input_file_name == '-help') then
          call command_line_help(trim(program_name))
       endif
       if(number_of_arguments == 1) then
          conv_level = 1
       else
          call Get_Command_Argument(2, conv_test)
          if(conv_test == '-conv_test') then
             attach = 'convres'
             if(number_of_arguments == 2) then
                conv_level = 3
             elseif(number_of_arguments == 3) then
                call Get_Command_Argument(3, conv_level_char)
                conv_level = iachar(conv_level_char) - iachar("0")
                if (conv_level < 0 .or. conv_level > 8 .or. &
                     & len(trim(conv_level_char)) > 1 ) &
                     & call command_line_help(trim(program_name))
             endif
          else
             call command_line_help(trim(program_name))
          endif
       endif
    endif

    write(*, *) "Trying to open the input file: ", trim(input_file_name), " ..."
    open(4, err = 120, file = trim(input_file_name), status = 'old', action = 'read')
    write(*, *) " successful"
!    read(4, nml = param, end = 110, err = 130)
    read(4, nml = param)
    close(4, status = 'keep')

    return
    110 call end_of_file
    120 call input_file_error
    130 call parameter_file_error
    
    end subroutine read_params

    subroutine end_of_file
        write (*, *) "The file doesn't contain all parameters!"
        stop 'File error'
    end subroutine end_of_file

    subroutine input_file_error
        write (*, *) "File error: probably the file name is wrong."
        stop 'File does not exist!'
    end subroutine input_file_error

    subroutine parameter_file_error
        write (*, *) "File error:"
        stop 'There is a problem in the parameter list'
      end subroutine parameter_file_error


   subroutine command_line_help(programname)
      character(len=*), intent(in) :: programname
      write(*,*)
      write(*,*) "Proper usage: ", programname, " [filename or switch] [-conv_test] [number] "
      write(*,*)
      write(*,*) "switch: ?, -?, h, -h, help, -help "
      write(*,*) "        --- this calls the help you see now."
      write(*,*)
      write(*,*) "filename"
      write(*,*) "        --- the name of a file containing parameters for the program."
      write(*,*)
      write(*,*) "-conv_test"
      write(*,*) "         --- tells me to make a convergence test."
      write(*,*)
      write(*,*) "number" 
      write(*,*) "        --- gives the level of convergence test, has to be between 1 and 9. If no number given, the convergence test will be done with 3 levels"
      stop
  end subroutine command_line_help
end module command_line
