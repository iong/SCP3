program main
  use pes_shell
  implicit none

  real,dimension(:),allocatable::xx,gd
  character(len=2)::symb
  character(len=32)::filename
  integer::i,natm
  real::pot,dip(3)

  call getarg(1,filename)
  open(21,status='old',file=filename)
  read(21,*) natm
  read(21,*)
  allocate(xx(3*natm))
  allocate(gd(3*natm))

  do i=1,natm
     read(21,*) symb,xx(i*3-2:i*3)
  end do
  xx=xx/auang
  call pes_init(natm/3)

  pot=f(xx)
!  call dipole(xx,dip)
!  gd=grad(xx,1.e-3)

  write(*,'(A,E)') "POT (hartree) =",pot
!  write(*,'(A,3F9.5)') "dip (a.u.)    =",dip
!  write(*,'(A)') "gradient      ="
!  print *,gd



end program main

