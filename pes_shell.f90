module pes_shell
  use pot_monomer_mod
  implicit none
  ! constants
  real,parameter::auang=0.5291772083
  real,parameter::aucm=219474.6313710
  real,parameter::aukcal=627.51

  ! parameters for switching function
  real,parameter::r2i=6.5/auang
  real,parameter::r2f=7.5/auang
  real,parameter::r3i=5.0/auang
  real,parameter::r3f=6.0/auang

  ! some global variables
  integer,dimension(:,:),allocatable::idx_3b,idx_2b

contains
  include 'nbody.f90'
  !=================================================!
  ! Initializing HBB water potential                !
  !=================================================!
  subroutine pes_init(nw)
    integer,intent(in)::nw
    !::::::::::::::::::::

    ! 3-body init
    call pes_init_3b('./data/WHBB_MP2_3b5')
    allocate(idx_3b(9,nw*(nw-1)*(nw-2)/6))
    call map_3b(nw,idx_3b)

    ! 2-body init
    call prepot('./data/fit0000.out.pes4X1-T1.090915')
    call predip('./data/h4o2.dms4.coeff.dat')
    allocate(idx_2b(6,nw*(nw-1)/2))
    call map_2b(nw,idx_2b)

    ! monomer init
    call monomer_init()

    return
  end subroutine pes_init

  !==================================================!
  ! water dipole moment                              !
  !==================================================!
  subroutine dipole(x,dp)
    real,dimension(:),intent(in)::x
    real,dimension(3),intent(out)::dp
    ! ::::::::::::::::::::
    real,dimension(1:size(x,1),3)::xn
    integer::natm

    natm=size(x)/3
    xn=reshape(x,(/3,natm/))

    call dip12bhbb(natm,xn,dp)

    return
  end subroutine dipole

  !==================================================!
  ! water potential                                  !
  !   x: cartesian coordinates in bohr               !
  !  im: optional augument should only be used for   !
  !      efficient numerical gradient calculation    !
  ! * To obtain the full potential of water cluster, !
  !   one should call this function with ONLY ONE    ! 
  !   argument, simply as f(x)                       !
  !==================================================!
  function f(x,xc,im) result(pot)
    real,dimension(:),intent(in)::x
    real,dimension(:),optional,intent(in)::xc
    integer,optional,intent(in)::im
    real::pot
    ! ::::::::::::::::::::
    real,dimension(3,1:size(x)/3)::xn,xnc
    real::p1,p2,p3,pmb
    integer::natm

    natm=size(x)/3
    xn=reshape(x,(/3,natm/))
    if (present(xc)) then
       xnc=reshape(xc,(/3,natm/))
    end if

    p1=0.d0; p2=0.d0; p3=0.d0; pmb=0.d0

    if (present(im))then
       call pot1b_g(natm,xn,p1,im)
       call potc2bt1_g(natm,xn,p2,im)
       call potc3b_gs(natm,xn,xnc,p3,im)
    else 
       call pot1b(natm,xn,p1)
       call potc2bt1(natm,xn,p2)
       call potc3b(natm,xn,p3)
    end if
    call potmb(natm,xn,pmb)
    pot=p1+p2+p3+pmb

    return
  end function f
  !==================================================!
  ! numerical gradient of the water potential        !
  !==================================================!  
  function grad(x,eps) result(gd)
    real,dimension(:),intent(in)::x
    real,intent(in)::eps
    !::::::::::::::::::::
    real,dimension(1:size(x))::xt,gd
    real::fa,fb,fc,fd
    integer::dim,i,iacc,im

    dim=size(x)
!$omp parallel do private(im,xt,fa,fb,fc,fd)
    do i=1,dim
       ! im: index of the atom 
       im=(i+2)/3
       xt=x; xt(i)=xt(i)-eps
       fa=f(xt,x,im)
       xt=x; xt(i)=xt(i)+eps
       fb=f(xt,x,im)
       gd(i)=0.5d0*(fb-fa)/eps
    end do
!$omp end parallel do

    return
  end function grad
  
end module pes_shell
