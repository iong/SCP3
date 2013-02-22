MODULE water_module

  !DOUBLE PRECISION, PARAMETER :: epsilon = 0.1852D0            ! kcal/mol
  DOUBLE PRECISION, PARAMETER :: bohr = 0.52917721092           ! atomic units
  DOUBLE PRECISION, PARAMETER :: autokJmol =2625.5002
  DOUBLE PRECISION, PARAMETER :: epsilon = 0.0002951349           ! atomic units
  DOUBLE PRECISION, PARAMETER :: sigma = 3.1589D0/bohr             ! bohr
  DOUBLE PRECISION, PARAMETER :: qM = 1.1128D0                 ! |e|
  DOUBLE PRECISION, PARAMETER :: gamma = 0.73612D0             ! (unit-less)
  !DOUBLE PRECISION, PARAMETER :: Dr = 116.09D0                 ! kcal/mol
  DOUBLE PRECISION, PARAMETER :: Dr = 0.1850012                 ! atomic units
  DOUBLE PRECISION, PARAMETER :: alphar = 2.287*bohr             ! bohr**(-1)
  !DOUBLE PRECISION, PARAMETER :: alphar = 2.287D0             ! angstrom**(-1)
  DOUBLE PRECISION, PARAMETER :: req =0.9419D0 /bohr                ! bohr
  !DOUBLE PRECISION, PARAMETER :: ktheta = 87.85D0              ! kcal/(mol*rad**2)
  DOUBLE PRECISION, PARAMETER :: ktheta = 0.139998              ! atomic units/(rad**2)
  !DOUBLE PRECISION, PARAMETER :: thetaeq = 107.4D0            ! degrees
  DOUBLE PRECISION, PARAMETER :: thetaeq = 1.87448361664D0      ! radians; 107.4D0*pi/180

CONTAINS


  SUBROUTINE TIP4P_UF(q, Utot, force_out)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)  :: q(:,:)          ! coordinates; (Ox, Oy, Oz, H1x, H1y, H1z, H2x, H2y, H2z)
    DOUBLE PRECISION, INTENT(OUT) :: Utot, force_out(:,:) ! energy, force

    DOUBLE PRECISION, DIMENSION(size(q, 2), 3) :: q_block
    DOUBLE PRECISION, DIMENSION(size(q, 2), 3) :: F

    INTEGER :: i

    q_block = TRANSPOSE(q)
    Utot = 0.0d0
    F = 0d0

    ! inter
    CALL Coulomb(q_block, Utot, F)
    CALL LJ_UF(q_block, Utot, F)


    ! intra
    CALL TIP4P_UF_intra(size(q, 2)/3, q_block, Utot, F)

    force_out = TRANSPOSE(F)
  END SUBROUTINE


  FUNCTION TIP4P_U(q) result(Utot)
    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(:, :), INTENT(IN) :: q     ! coordinates; (Ox, Oy, Oz, H1x, H1y, H1z, H2x, H2y, H2z)
    DOUBLE PRECISION :: Utot

    DOUBLE PRECISION, allocatable :: q_block(:,:)

    allocate ( q_block(size(q, 2) , 3) )
    q_block = TRANSPOSE(q)
    Utot = 0.0d0

    ! inter
    CALL Coulomb(q_block, Utot)
    Utot = Utot +  LJ(q_block)

    ! intra
    Utot = Utot + TIP4P_U_intra(q_block)

    deallocate(q_block)

  END FUNCTION


  SUBROUTINE Coulomb(r, Utot, F)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(INOUT) :: r(:,:)      ! coordinates
    DOUBLE PRECISION, INTENT(INOUT) :: Utot                    ! intermolecular potential
    DOUBLE PRECISION, INTENT(INOUT), optional :: F(:,:)      ! intermolecular forces

    INTEGER :: i, Natoms
    DOUBLE PRECISION, ALLOCATABLE :: q(:), rO(:,:), FCoulomb(:,:)      ! charges, O coordinates

    Natoms = size(r, 1)
    ALLOCATE(q(Natoms), rO(Natoms/3 , 3))

    q(1::3) = -qM
    q(2::3) = 0.5d0*qM
    q(3::3) = 0.5d0*qM

    rO = r(1::3,:)
    r(1::3,:) = gamma*r(1::3,:) + 0.5*(1.0-gamma)*(r(2::3,:) + r(3::3,:))

    if ( present(F) ) then

        allocate( FCoulomb(Natoms, 3) )

        FCoulomb = 0d0

        CALL Coulomb_UF(q, r, Utot, FCoulomb)

        FCoulomb(2::3,:) = FCoulomb(2::3,:) + 0.5*(1.0-gamma) * FCoulomb(1::3,:)
        FCoulomb(3::3,:) = FCoulomb(3::3,:) + 0.5*(1.0-gamma) * FCoulomb(1::3,:)
        FCoulomb(1::3,:) = gamma * FCoulomb(1::3,:)

        F = F + Fcoulomb

        DEALLOCATE(FCoulomb)
    else

        Utot = Utot + Coulomb_U(q, r)

    end if

    r(1::3,:) = rO

    DEALLOCATE(q, rO)
  END SUBROUTINE


  FUNCTION LJ(r) result(Utot)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: r(:,:)
    DOUBLE PRECISION :: Utot

    DOUBLE PRECISION, DIMENSION(size(r, 1)/3) :: rsq, ir2, sir6, sir12
    DOUBLE PRECISION, DIMENSION(size(r, 1)/3, 3) ::  dr
    INTEGER :: i, k,j, NO
    
    NO = size(r, 1)/3
    Utot = 0d0
    DO i=1,NO-1
       DO k=1,3
          dr(i+1:NO,k) = r(3*i+1::3,k) - r(3*i-2,k)
       END DO
       rsq(i+1:NO) = dr(i+1:NO,1)**2 + dr(i+1:NO,2)**2 + dr(i+1:NO,3)**2
       ir2(i+1:NO) = 1.0/rsq(i+1:NO)


       sir6(i+1 : NO) = (sigma**2 * ir2(i+1:NO))**3
       sir12(i+1:NO) = sir6(i+1:NO)**2
       Utot  = Utot + 4.0 *epsilon* SUM(sir12(i+1:NO) - sir6(i+1:NO))
    END DO
  END FUNCTION


  SUBROUTINE LJ_UF(r, Utot, F)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: r(:,:)
    DOUBLE PRECISION, INTENT(inout) :: F(:,:), Utot

    DOUBLE PRECISION, DIMENSION(size(r, 1)/3) :: rsq, ir2, sir6, sir12, Cij
    DOUBLE PRECISION, DIMENSION(size(r, 1)/3, 3) ::  dr, Fij
    INTEGER :: i, k,j, NO

    NO = size(r, 1)/3
    DO i=1,NO-1
       DO k=1,3
          dr(i+1:NO,k) = r(3*i+1::3,k) - r(3*i-2,k)
       END DO
       rsq(i+1:NO) = dr(i+1:NO,1)**2 + dr(i+1:NO,2)**2 + dr(i+1:NO,3)**2
       ir2(i+1:NO) = 1.0/rsq(i+1:NO)


       sir6(i+1 : NO) = (sigma**2 * ir2(i+1:NO))**3
       sir12(i+1:NO) = sir6(i+1:NO)**2
       Utot  = Utot + 4.0 *epsilon* SUM(sir12(i+1:NO) - sir6(i+1:NO))

       Cij(i+1:NO) = 4.0*epsilon*(12.0 * sir12(i+1:NO)  - 6.0 * sir6(i+1:NO)) * ir2(i+1:NO)
       DO k=1,3
          Fij(i+1:NO,k) = dr(i+1:NO,k) * Cij(i+1:NO)
       END DO


       F(3*i-2,:) = F(3*i-2,:) - SUM(Fij(i+1:NO,:), 1)
       F(3*i+1::3,:) = F(3*i+1::3,:) +  Fij(i+1:NO,:)
    END DO
  END SUBROUTINE

  FUNCTION Coulomb_U(q, r) RESULT(Utot)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: q(:), r(:,:)
    DOUBLE PRECISION :: Utot

    INTEGER :: Natoms
    DOUBLE PRECISION, DIMENSION(size(q)) :: rsq, ir2, Cij
    DOUBLE PRECISION, DIMENSION(size(q), 3) ::  dr
    INTEGER :: i, k

    Utot = 0d0
    Natoms = size(q)
    DO i=1,Natoms-3
       DO k=1,3
          dr(i+1:Natoms,k) = r(i+1:Natoms,k) - r(i,k)
       END DO

       rsq(i+1:Natoms) = dr(i+1:Natoms,1)**2 + dr(i+1:Natoms,2)**2 + dr(i+1:Natoms,3)**2
       ir2(i+1:Natoms) = 1.0/rsq(i+1:Natoms)

       Cij(i+1 : Natoms) = q(i) * q(i+1:Natoms) * SQRT(ir2(i+1:Natoms))

       Cij(i : 3*((i + 2)/3)) = 0.0

       Utot  = Utot + SUM(Cij(i+1:Natoms))
    END DO
  END FUNCTION


  SUBROUTINE Coulomb_UF(q, r, Utot, F)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: q(:), r(:,:)
    DOUBLE PRECISION, INTENT(out) :: F(:,:), Utot

    INTEGER :: Natoms
    DOUBLE PRECISION, DIMENSION( size(q) ) :: rsq, ir2, Cij
    DOUBLE PRECISION, DIMENSION( size(q) , 3) ::  dr, Fij
    INTEGER :: i, k

    Natoms = size(q)
    DO i=1,Natoms-3
       DO k=1,3
          dr(i+1:Natoms,k) = r(i+1:Natoms,k) - r(i,k)
       END DO

       rsq(i+1:Natoms) = dr(i+1:Natoms,1)**2 + dr(i+1:Natoms,2)**2 + dr(i+1:Natoms,3)**2
       ir2(i+1:Natoms) = 1.0/rsq(i+1:Natoms)

       Cij(i+1 : Natoms) = q(i) * q(i+1:Natoms) * SQRT(ir2(i+1:Natoms))

       Cij(i : 3*((i + 2)/3)) = 0.0

       Utot  = Utot + SUM(Cij(i+1:Natoms))

       Cij(i+1:Natoms) = Cij(i+1:Natoms) * ir2(i+1:Natoms)
       DO k=1,3
          Fij(i+1:Natoms,k) = dr(i+1:Natoms,k) *Cij(i+1:Natoms)
       END DO

       F(i,:) = F(i,:) - SUM(Fij(i+1:Natoms,:), 1)
       F(i+1:Natoms,:) = F(i+1:Natoms,:) +  Fij(i+1:Natoms,:)
    END DO
  END SUBROUTINE


  function TIP4P_U_intra(coords) result(U_intra)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: coords(:,:)           ! coordinates
    DOUBLE PRECISION :: U_intra

    INTEGER :: i, k, Natoms
    DOUBLE PRECISION :: U1, U2, Uth
    DOUBLE PRECISION :: r1, r2, r1dotr2, rtilde, theta, denom
    DOUBLE PRECISION, DIMENSION(3) :: r1vec, r2vec
    DOUBLE PRECISION, DIMENSION(3,3) :: q


    U_intra = 0.0d0
    Natoms = size(coords, 1)
    DO i = 1, Natoms, 3


       q = (coords(i:i+2, :))
       
       r1vec = q(2,:) - q(1,:)
       r2vec = q(3,:) - q(1,:)

       r1dotr2= r1vec(1)*r2vec(1) + r1vec(2)*r2vec(2) + r1vec(3)*r2vec(3)
       r1 = SQRT(r1vec(1)**2 + r1vec(2)**2 + r1vec(3)**2) 
       r2 = SQRT(r2vec(1)**2 + r2vec(2)**2 + r2vec(3)**2)

       rtilde = r1dotr2/(r1*r2)
       theta = ACOS(rtilde)
       denom = SQRT(1.0d0-rtilde**2)


       ! intramolecular potential energy
       U1 = (Dr * alphar**2 * (r1 - req)**2) * (1.0d0 - alphar*(r1-req)*(1.0d0 - (7.0d0/12.0d0)*alphar*(r1-req)))
       U2 = (Dr * alphar**2 * (r2 - req)**2) * (1.0d0 - alphar*(r2-req)*(1.0d0 - (7.0d0/12.0d0)*alphar*(r2-req)))
       Uth = (ktheta/2.0D0) * (theta - thetaeq)**2

       U_intra = U_intra + (U1 + U2 + Uth)
    END DO
  END function


  SUBROUTINE TIP4P_UF_intra(NO, coords, U_intra, F)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: NO                                            ! number of water molecules
    DOUBLE PRECISION, DIMENSION(3*NO, 3), INTENT(IN) :: coords           ! coordinates
    DOUBLE PRECISION, DIMENSION(3*NO, 3), INTENT(INOUT) :: F       ! intramolecular forces
    DOUBLE PRECISION, INTENT(INOUT) :: U_intra              

    INTEGER :: i, k
    DOUBLE PRECISION :: U1, U2, Uth
    DOUBLE PRECISION :: r1, r2, r1dotr2, rtilde, theta, denom
    DOUBLE PRECISION, DIMENSION(3) :: r1vec, r2vec
    DOUBLE PRECISION, DIMENSION(3,3) :: q, F1, F2, Fth


    U1 = 0.0d0
    U2 = 0.0d0
    Uth = 0.0d0
    F1 = 0.0d0
    F2 = 0.0d0
    Fth = 0.0d0

    DO i = 1, 3*(NO-1)+1, 3


       q = (coords(i:i+2, :))
       
       r1vec = q(2,:) - q(1,:)
       r2vec = q(3,:) - q(1,:)

       r1dotr2= r1vec(1)*r2vec(1) + r1vec(2)*r2vec(2) + r1vec(3)*r2vec(3)
       r1 = SQRT(r1vec(1)**2 + r1vec(2)**2 + r1vec(3)**2) 
       r2 = SQRT(r2vec(1)**2 + r2vec(2)**2 + r2vec(3)**2)

       rtilde = r1dotr2/(r1*r2)
       theta = ACOS(rtilde)
       denom = SQRT(1.0d0-rtilde**2)


       ! intramolecular potential energy
       U1 = (Dr * alphar**2 * (r1 - req)**2) * (1.0d0 - alphar*(r1-req)*(1.0d0 - (7.0d0/12.0d0)*alphar*(r1-req)))
       U2 = (Dr * alphar**2 * (r2 - req)**2) * (1.0d0 - alphar*(r2-req)*(1.0d0 - (7.0d0/12.0d0)*alphar*(r2-req)))
       Uth = (ktheta/2.0D0) * (theta - thetaeq)**2

       U_intra = U_intra + (U1 + U2 + Uth)


       !Gradient associated with O-H1 stretch
       F1(1,:) = (Dr*alphar**2*(r1-req)  *  (2.0d0 - 3.0d0*alphar*(r1-req) + (7.0d0/3.0d0) * alphar**2 * (r1-req)**2) )/r1
       F1(2,:) = F1(1,:)

       DO k = 1, 3
          F1(1, k) = F1(1, k) * (-(q(2,k)-q(1,k)))
          F1(2, k) = F1(2, k) * (q(2,k)-q(1,k))
       END DO


       !Gradient associated with O-H2 stretch
       F2(1,:) = (Dr*alphar**2*(r2-req)  *  (2.0d0 - 3.0d0*alphar*(r2-req) + (7.0d0/3.0d0) * alphar**2 * (r2-req)**2) )/r2
       F2(3,:) = F2(1,:)      

       DO k = 1, 3
          F2(1, k) = F2(1, k) * (-(q(3,k)-q(1,k)))
          F2(3, k) = F2(3, k) * (q(3,k)-q(1,k))
       END DO


       ! Gradient associated with H1-O-H2 bond angle
       DO k = 1, 3
          Fth(1, k) = (2.0d0*q(1,k)-q(2,k)-q(3,k))/(r1 * r2) &
               + (q(3,k) - q(1,k))*r1dotr2/(r1 * r2**3) &
               + (q(2,k) - q(1,k))*r1dotr2/(r1**3 * r2)

          Fth(2, k) = (q(3,k)-q(1,k))/(r1*r2) &
               - (q(2,k) - q(1,k))*r1dotr2/(r1**3 * r2)

          Fth(3, k) = (q(2,k)-q(1,k))/(r1*r2) &
               - (q(3,k) - q(1,k))*r1dotr2/(r1 * r2**3)
              
       END DO

       Fth = (-ktheta*(theta - thetaeq)/denom) * Fth

       ! intramolecular force
       F(i:i+2, :) = F(i:i+2, :) -(F1 + F2 + Fth)

    END DO
  END SUBROUTINE


  SUBROUTINE test_grad(NO,q)
    ! gives the force using finite differences

    IMPLICIT NONE
    INTEGER :: i, NO, k
    REAL(8) :: s,Ener,Ener0, qold
    REAL(8), DIMENSION(3,3*NO) :: F, F0, q

    s=1d-4

    CALL TIP4P_UF(q, Ener, F0)
    DO i=1,3*NO
       do k=1,3
            qold = q(k,i)
            q(k,i)=qold - s

            CALL TIP4P_UF(q, Ener0, F)

            q(k,i)=qold + s
            CALL TIP4P_UF(q, Ener, F)

            F(k, i)=(Ener0-Ener)/(2*s)
            q(k, i)=qold
            WRITE(*,*) F(k, i), F0(k, i)
        end do
    ENDDO


  END SUBROUTINE test_grad
END MODULE water_module


!! Interoperability with C/C++
subroutine TIP4P_UF(N_H2O, r, U, UX) bind(C, name="TIP4P_UF")
    use water_module, only: fortran_TIP4P_UF => TIP4P_UF
    use iso_c_binding
    integer(c_int), value :: N_H2O
    real(c_double) :: r(3, 3*N_H2O), U, UX(3, 3*N_H2O)
    
    call fortran_TIP4P_UF(r, U, UX)
end subroutine

function TIP4P_U(N_H2O, r) result(Utot) bind(C, name="TIP4P_U")
    use water_module, only: fortran_TIP4P_U => TIP4P_U
    use iso_c_binding
    integer(c_int), value :: N_H2O
    real(c_double) :: r(3, 3*N_H2O), Utot
    
    Utot = fortran_TIP4P_U(r)
end function
