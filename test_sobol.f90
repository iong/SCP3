program test_sobol
    use sobol
    implicit none

    character(256) :: arg

    integer :: seqlen=200000
    integer :: i, melout = 34
    integer(8) :: lp2, dim=18, sobol_skip = 100000

    double precision, allocatable :: q(:)

    if (command_argument_count() == 2) then
        call get_command_argument(1, arg)
        read(arg, *) sobol_skip

        call get_command_argument(2, arg)
        read(arg, *) seqlen
    end if

    lp2 = 2**int(log(real(sobol_skip))/log(2d0))
    print *, 'distance from lowest power of 2 =', sobol_skip - lp2

    allocate(q(dim))

    open(melout, file='test_sobol.dat')
    do i=1,seqlen    
        q = i8_sobol(dim, sobol_skip)
        !call random_number(q)
        write(melout, *) q(1)*q(11)*q(16)*q(18)
    end do
    close(melout)
end program
