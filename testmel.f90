program testmel
    use water0, only: TIP4P
    use bsm
    use sobol
    implicit none

    character(256) :: arg

    integer :: i, j, info, melout = 34

    integer :: N, NO, Nmodes

    integer(8) :: sobol_skip=2**30, skip2, NSobol = 2**20

    double precision, allocatable :: r0(:,:), r(:,:), mass(:), Vr(:,:), H(:,:), &
            omegasq(:), omega(:), alpha(:), MUa(:,:), q(:), Vq(:)
    double precision, allocatable :: work(:)

    double precision :: s, t, V, V0, VHO,rn

    double precision, parameter :: bohr=0.52917721092, autocm=2.194746313e5, qM = 1.1128

    call get_command_argument(1, arg)
    call load_xyz_and_hessian(arg, r0, mass, H)

    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read(arg, *) NSobol

        sobol_skip = 2**ceiling(log(real(NSobol))/log(2d0), 8)
        sobol_skip = NSobol
        print *, 'sobol_skip =', sobol_skip
    end if

    N = size(r0, 2)
    Nmodes = 3*N - 6

    allocate(r(3, N), Vr(3, N), omegasq(3*N), omega(Nmodes), alpha(Nmodes), MUa(3*N, Nmodes), q(Nmodes),  Vq(Nmodes) )

    allocate (work(3*3*N))
    call dsyev('V', 'U', 3*N, H, 3*N, omegasq, work, size(work), info)

    omega = sqrt(omegasq(7:3*N))
    alpha = sqrt(omega)

    do i=1,Nmodes
        MUa(:,i) = H(:,i+6) / (sqrt(mass) * alpha(i) )
    end do

    open(melout, file='testmel.dat')
    s = 0d0
    do i=1,NSobol
        call random_number(rn)
        skip2 = sobol_skip+rn*Nsobol
        call sobol_stdnormal(sobol_skip, q)
        q = q / sqrt(2d0)
        r = r0
        call dgemv('N', 3*N, Nmodes, 1d0, MUa, 3*N, q, 1, 1d0, r, 1)
        call TIP4P(N/3, r, V0, Vr)

        Vq = omega * q
        call dgemv('T', 3*N, Nmodes, -1d0, MUa, 3*N, Vr, 1, -1d0, Vq, 1)

        VHO = 0.5d0*sum(omega*q**2)
        V = V0 - VHO

        t =  m2_n1_f_m1_l1_j1(1, 18, 11, 16)
        s = s + t
        write(melout, '(I10,2(" ",ES15.8))') i, s/real(i), t
    end do

    close(melout)
contains

subroutine load_xyz_and_hessian(fname, r, mass, H)
    implicit none
    real*8, intent(out), allocatable :: r(:,:), mass(:), H(:,:)
    character(LEN=*), intent(in) :: fname
    character(4):: species
    integer :: i, j, N

    real(8) :: M
    real(8), parameter :: Hmass=1837.15137, Cmass=21891.6543

    if (allocated(r)) then
        deallocate (r)
    end if

    if (allocated(mass)) then
        deallocate (mass)
    end if

    if (allocated(H)) then
        deallocate (H)
    end if



    open(33,file=fname,STATUS='OLD')
    read(33,*) N
    read(33, *)

    allocate( r(3,N), mass(3*N), H(3*N, 3*N) )

    do j=1,N
        read(33, *) species, r(:, j)
        if (species == 'H') then
            M = Hmass
        else if (species == 'C') then
            M = Cmass
        else if (species == 'O') then
            M = Cmass*15.9994d0/12.0107d0
        else
            M = -1d0
        end if
        mass(3*j - 2 : 3*j) = M
    end do

    r = r / bohr

    read(33,*) ((H(i,j),i=1,j),j=1,3*N)
    
    close(33)
end subroutine

function m2_n1_f_m1_l1_j1(m, n, l, j) result(res)
    implicit none
    integer, intent(in) :: m, n, l,  j
    double precision :: res
    res = q(l)*q(m)*q(n)*Vq(j) + q(j)*( q(m)*q(n)*Vq(l) + q(l)*q(n)*Vq(m) + q(l)*q(m)*Vq(n) )
    res = res * (-1.0 + 2.0*q(m)**2)
    res = res + 4.0*V*q(j)*q(l)*q(m)*q(n)
    res = res / (2d0*sqrt(2d0))
end function

end program
