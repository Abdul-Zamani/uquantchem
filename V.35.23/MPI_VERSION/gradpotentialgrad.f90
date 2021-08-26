SUBROUTINE gradpotentialgrad(BAS,NATOMS,ATOMS,GVG)
      ! This subroutine calculates the the following 
      ! matrix elements <grad(PHI_i) | V | grad(PHI_j) > . 
      ! See for example A. Wolf, M. Reiner, B. A. Hess Chapter 11, 
      ! page 651, eqn (113) in Chapter 11 of Relativistic Electronic Structure Theory,
      ! part 1: Fundamentals, Theoretical and Computational Chemistry, Vol 11.
      ! Or J. Chem. Phys., 117, 9215 (2002), equation (A5).
      ! The implementation is that of my blue folder, "Jonisering + DKH + Atomic Sphere", 
      ! page 10. Ultimately this matrix will be used in the Douglas-Kroll transformation
      ! of the Dirac Hamiltonian into the approximately block diagonalized DKH hamitonian 
      ! to be used in scalar relativistic calculations. See J. Chem. Phys., 117, 9215 (2002)
      USE datatypemodule
      IMPLICIT NONE
      INCLUDE "mpif.h"
      DOUBLE PRECISION, EXTERNAL :: primpotential
      INTEGER, INTENT(IN) :: NATOMS
      TYPE(BASIS), INTENT(IN) :: BAS
      TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
      DOUBLE PRECISION, INTENT(OUT) :: GVG(BAS%NBAS,BAS%NBAS)
      INTEGER :: I,J,N,K,M,L1,M1,N1,L2,M2,N2
      INTEGER :: ierr
      DOUBLE PRECISION :: NO1,NO2,NP1,NP2
      DOUBLE PRECISION :: A(3),B(3),alpha1,alpha2,coeff1,coeff2

      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      GVG(:,:) = 0.0d0

      DO I=1,BAS%NBAS
                
                L1= BAS%PSI(I)%L(1)
                M1= BAS%PSI(I)%L(2)
                N1= BAS%PSI(I)%L(3)
                A = BAS%PSI(I)%R
                NO1 = BAS%PSI(I)%NORM 
                
                DO J=I,BAS%NBAS

                        GVG(I,J) = 0.0d0
                        
                        L2= BAS%PSI(J)%L(1)
                        M2= BAS%PSI(J)%L(2)
                        N2= BAS%PSI(J)%L(3)
                        B = BAS%PSI(J)%R
                        NO2 = BAS%PSI(J)%NORM
                        
                        DO N=1,BAS%PSI(I)%NPRIM
                                alpha1 = BAS%PSI(I)%EXPON(N)
                                coeff1 = BAS%PSI(I)%CONTRCOEFF(N)
                                NP1 = BAS%PSI(I)%PRIMNORM(N)
                                DO K=1,BAS%PSI(J)%NPRIM
                                        alpha2 = BAS%PSI(J)%EXPON(K)
                                        coeff2 = BAS%PSI(J)%CONTRCOEFF(K)
                                        NP2 = BAS%PSI(J)%PRIMNORM(K)
                                        DO M=1,NATOMS
                                                !-----------------------------------------------
                                                ! Here we calculate the gradpotentialgrad matrix
                                                !-----------------------------------------------
                                                
                                                ! x-coordinate:
                                                GVG(I,J) = GVG(I,J) - &
                                                & 4.0d0*alpha1*alpha2*ATOMS(M)%Z*NO1*NO2*coeff1*coeff2*NP1*NP2* &
                                                & primpotential(L1+1,M1,N1,A,alpha1,ATOMS(M)%R,L2+1,M2,N2,B,alpha2)
                                                
                                                ! y-coordinate:
                                                GVG(I,J) = GVG(I,J) - &
                                                & 4.0d0*alpha1*alpha2*ATOMS(M)%Z*NO1*NO2*coeff1*coeff2*NP1*NP2* &
                                                & primpotential(L1,M1+1,N1,A,alpha1,ATOMS(M)%R,L2,M2+1,N2,B,alpha2)
                                                
                                                ! z-coordinate:
                                                GVG(I,J) = GVG(I,J) - &
                                                & 4.0d0*alpha1*alpha2*ATOMS(M)%Z*NO1*NO2*coeff1*coeff2*NP1*NP2* &
                                                & primpotential(L1,M1,N1+1,A,alpha1,ATOMS(M)%R,L2,M2,N2+1,B,alpha2)
                                              

                                                ! x-coordinate:
                                                IF ( L1 .GT. 0 ) THEN
                                                        GVG(I,J) = GVG(I,J) + &
                                                        & 2.0d0*alpha2*L1*ATOMS(M)%Z*NO1*NO2*coeff1*coeff2*NP1*NP2* &
                                                        & primpotential(L1-1,M1,N1,A,alpha1,ATOMS(M)%R,L2+1,M2,N2,B,alpha2)
                                                ENDIF
                                                
                                                IF ( L2 .GT. 0 ) THEN
                                                        GVG(I,J) = GVG(I,J) + &
                                                        & 2.0d0*alpha1*L2*ATOMS(M)%Z*NO1*NO2*coeff1*coeff2*NP1*NP2* &
                                                        & primpotential(L1+1,M1,N1,A,alpha1,ATOMS(M)%R,L2-1,M2,N2,B,alpha2)
                                                ENDIF

                                                IF ( L1 .GT. 0 .AND. L2 .GT. 0 ) THEN
                                                        GVG(I,J) = GVG(I,J) - &
                                                        & L1*L2*ATOMS(M)%Z*NO1*NO2*coeff1*coeff2*NP1*NP2* &
                                                        & primpotential(L1-1,M1,N1,A,alpha1,ATOMS(M)%R,L2-1,M2,N2,B,alpha2)
                                                ENDIF

                                                ! y-coordinate:
                                                IF ( M1 .GT. 0 ) THEN
                                                        GVG(I,J) = GVG(I,J) + &
                                                        & 2.0d0*alpha2*M1*ATOMS(M)%Z*NO1*NO2*coeff1*coeff2*NP1*NP2* &
                                                        & primpotential(L1,M1-1,N1,A,alpha1,ATOMS(M)%R,L2,M2+1,N2,B,alpha2)
                                                ENDIF
                                                
                                                IF ( M2 .GT. 0 ) THEN
                                                        GVG(I,J) = GVG(I,J) + &
                                                        & 2.0d0*alpha1*M2*ATOMS(M)%Z*NO1*NO2*coeff1*coeff2*NP1*NP2* &
                                                        & primpotential(L1,M1+1,N1,A,alpha1,ATOMS(M)%R,L2,M2-1,N2,B,alpha2)
                                                ENDIF

                                                IF ( M1 .GT. 0 .AND. M2 .GT. 0 ) THEN
                                                        GVG(I,J) = GVG(I,J) - &
                                                        & M1*M2*ATOMS(M)%Z*NO1*NO2*coeff1*coeff2*NP1*NP2* &
                                                        & primpotential(L1,M1-1,N1,A,alpha1,ATOMS(M)%R,L2,M2-1,N2,B,alpha2)
                                                ENDIF
                                                
                                                ! z-coordinate:
                                                IF ( N1 .GT. 0 ) THEN
                                                        GVG(I,J) = GVG(I,J) + &
                                                        & 2.0d0*alpha2*N1*ATOMS(M)%Z*NO1*NO2*coeff1*coeff2*NP1*NP2* &
                                                        & primpotential(L1,M1,N1-1,A,alpha1,ATOMS(M)%R,L2,M2,N2+1,B,alpha2)
                                                ENDIF
                                                
                                                IF ( N2 .GT. 0 ) THEN
                                                        GVG(I,J) = GVG(I,J) + &
                                                        & 2.0d0*alpha1*N2*ATOMS(M)%Z*NO1*NO2*coeff1*coeff2*NP1*NP2* &
                                                        & primpotential(L1,M1,N1+1,A,alpha1,ATOMS(M)%R,L2,M2,N2-1,B,alpha2)
                                                ENDIF

                                                IF ( N1 .GT. 0 .AND. N2 .GT. 0 ) THEN
                                                        GVG(I,J) = GVG(I,J) - &
                                                        & N1*N2*ATOMS(M)%Z*NO1*NO2*coeff1*coeff2*NP1*NP2* &
                                                        & primpotential(L1,M1,N1-1,A,alpha1,ATOMS(M)%R,L2,M2,N2-1,B,alpha2)
                                                ENDIF


                                        ENDDO
                                ENDDO
                        ENDDO

                        IF ( I .NE. J ) THEN 
                                GVG(J,I) = GVG(I,J)
                        ENDIF
                ENDDO

      ENDDO
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
END SUBROUTINE gradpotentialgrad
