!=================================================!
! this subroutine is used to initalize the        !
!   pes from bas                                  !
!=================================================!
subroutine pes_init_3b(dname)
  use pes
  character (len=*), intent(in) :: dname

  call pes0_init (dir=dname)
  call pes1_init ('x6y3')

  return
end subroutine pes_init_3b
!=========================================================!
! this is the main function to get the potential          !
! energy from fitting code                                !
!=========================================================!
function fpes(x) 
  use pes
  double precision,dimension(0:2,0:8),intent(in)::x
  double precision::fpes
  ! ::::::::::::::::::::

  fpes=pes_x6y3_pot(x)

  return
end function fpes
