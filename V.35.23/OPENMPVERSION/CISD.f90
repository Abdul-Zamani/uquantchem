SUBROUTINE CISD(Cup,Cdown,Ints,h,S,NB,Ne,E0,nuce,ECISD,WRITECICOEF,EHFeigenup,EHFeigendown,LEXCSP,NEEXC,ENEXCM,SPINCONSERVE)
      ! This subroutine constructs the many-electron basis set in the form
      ! of slater determinants, these slater determinants are in turn
      ! constructed from the single and double excitations of the
      ! unrestricted solutions to the Hartree-Fock equations. With this
      ! basis-set the hamiltonian matrix in the CI-basis is calculated
      ! from the slater rules (Cook Book page.71-74), and diagonalized.
      ! The electron repulsion tensor, nuclear attraction matrix and
      ! kinetic energy matrix are expressed in terms of primitive
      ! gaussian orbitals. The matrices are preconstructed together with
      ! the corresponding gaussian basis set expansion coefficients by
      ! the python program CISD.py, and consequently is needed as input
      ! to this program. By Petros Souvatzis Thu Nov 24 10:25:47 CET
      ! 2011

      IMPLICIT NONE
      EXTERNAL :: readmat,makedens,getJ,getK,diagh
      INTEGER, INTENT(IN) :: NB,Ne,NEEXC
      DOUBLE PRECISION, INTENT(IN) :: Cup(NB,NB),Cdown(NB,NB),Ints(NB,NB,NB,NB),h(NB,NB),S(NB,NB),E0,nuce
      DOUBLE PRECISION, INTENT(IN) :: EHFeigenup(NB),EHFeigendown(NB),ENEXCM
      LOGICAL, INTENT(IN) :: WRITECICOEF,LEXCSP,SPINCONSERVE
      DOUBLE PRECISION, INTENT(OUT) :: ECISD
      INTEGER :: I,J,N,M,NCI,K,L,K1,II,JJ,DIFFER,LU,LD,Q,A(2),B(2),Nel,NBl
      INTEGER, ALLOCATABLE :: EXC(:,:),V1(:)
      DOUBLE PRECISION :: TEMP1,TEMP2,TEMP3
      DOUBLE PRECISION, ALLOCATABLE :: PT(:,:),Pup(:,:),Pdown(:,:),HAMILTONIAN(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: COEFFUP(:,:),COEFFDOWN(:,:),Jup(:,:),Jdown(:,:),Kup(:,:),Kdown(:,:),Fup(:,:),Fdown(:,:),F(:,:),JT(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: Calpha(:),Ca(:),Cbeta(:),Cb(:),TEMPMAT(:,:),EIGENVAL(:),EIGENVECT(:,:),KT(:,:)
      DOUBLE PRECISION :: TOTALTIME
      LOGICAL :: EXCITE
      LOGICAL, EXTERNAL :: exciteornot
      
      ! Checking if the space of excitatons have been limited by the user or 
      ! if a full CISD calculation is to be performed
      IF ( LEXCSP ) THEN
              Nel = NEEXC
              
              IF ( Nel .GT. Ne ) Nel = Ne
              
              NBl = 0
              DO WHILE ( NBL .LT. NB .AND. ( EHFeigenup(NBl+1) .LT. ENEXCM .OR. EHFeigendown(NBl+1) .LT. ENEXCM) ) 
                        NBl = NBl + 1
              ENDDO
              IF ( NBl .LE. (Ne+MOD(Ne,2))/2 ) THEN
                      print*,'----------------------------------------------------------------------------'
                      print*,'     WARNING! The number HF-energy levels accesible for excitations = 0     '
                      print*,'ABORTING CISD calculation. Try increasing ENEXCM, or setting LEXCSP =.FALSE.'
                      print*,'----------------------------------------------------------------------------'
                      STOP
              ENDIF
              IF ( NEEXC .EQ. 0 ) THEN
                      print*,'------------------------------------------------------------------'
                      print*,'WARNING! Number of ground state lavels allowed for excitation = 0 '
                      print*,'   NEEXC must be greater than 1!   ABORTING CISD calculation.  '  
                      print*,'------------------------------------------------------------------'
                      STOP
              ENDIF
              IF ( NEEXC .EQ. 1 ) THEN
                      print*,'--------------------------------------------------------------------------'
                      print*,'WARNING! The number of accesible HF-groundstate levels for exitation is 1!'
                      print*,' Hence this will be only a SINGLES CI-CALCULATION!  ABORTING!             '
                      print*,'--------------------------------------------------------------------------'
                      STOP
              ENDIF

      ELSE
              Nel = Ne
              NBl = NB
      ENDIF
      
      ! Calculating the size of the CI-basis set (Singles and Doubles only )
      ! See page 12 in red folder
      IF ( SPINCONSERVE ) THEN
        NCI = (Nel-MOD(Nel,2))*(NBl-(Ne-MOD(Ne,2))/2)/2 + (Nel+MOD(Nel,2))*(NBl-(Ne+MOD(Ne,2))/2)/2 + 1
        NCI = NCI + (Nel-MOD(Nel,2))*(NBl-(Ne-MOD(Ne,2))/2)*(Nel+MOD(Nel,2))*(NBl-(Ne+MOD(Ne,2))/2)/4
        NCI = NCI + (Nel+MOD(Nel,2))*((Nel+MOD(Nel,2))/2-1)*(NBl - (Ne+MOD(Ne,2))/2 )*(NBl - (Ne+MOD(Ne,2))/2 - 1 )/8
        NCI = NCI + (Nel-MOD(Nel,2))*((Nel-MOD(Nel,2))/2-1)*(NBl - (Ne-MOD(Ne,2))/2 )*(NBl - (Ne-MOD(Ne,2))/2 - 1 )/8
      ELSE
        !NCI = Ne*(2*NB-Ne)+Ne*(Ne-1)*(2*NB-Ne)*(2*NB-Ne-1)/4 + 1
        NCI = 1 + Nel*(2*NBl-Ne)
        NCI = NCI + ((Nel-MOD(Nel,2))/2)*((Nel-MOD(Nel,2))/2 - 1)*(2*NBl-Ne)*(2*NBl-Ne-1)/4
        NCI = NCI + ((Nel+MOD(Nel,2))/2)*((Nel+MOD(Nel,2))/2 - 1)*(2*NBl-Ne)*(2*NBl-Ne-1)/4
        NCI = NCI + ((Nel-MOD(Nel,2))/2)*((Nel+MOD(Nel,2))/2)*(2*NBl-Ne)*(2*NBl-Ne-1)/2
     ENDIF
      
      print*,' ==========================================================='
      print*,'          NUMBER OF CI-SLATER-DETERMINANTS=',NCI
      print*,' ==========================================================='

      
      ALLOCATE(EXC(NCI,2*NBl),HAMILTONIAN(NCI,NCI),EIGENVAL(NCI),EIGENVECT(NCI,NCI))
      
      ! Setting the HF-ground state in the excitations matrix EXC:
      
      EXC(:,:) = 0
      
      M = (Ne-MOD(Ne,2))/2

      DO I=1,M
        EXC(1,I) = 1
        EXC(1,I+NBl) = 1
      ENDDO
      EXC(1,M+NBl+MOD(Ne,2)) = 1

      !---------------------------------------------------------------------
      ! (1) Storing all the different single excitations in the matrix EXC:
      !---------------------------------------------------------------------
      
      K1 = 1
      DO I=  M - (Nel-MOD(Nel,2))/2 + 1,M
        DO J=1,2*NBl
          IF ( J .GT. M .AND. ( J .LE. NBl .OR. J .GT. NBl+M+MOD(Ne,2) ) ) THEN
                  EXCITE = .FALSE.
                  IF ( .NOT. SPINCONSERVE ) THEN
                          EXCITE = .TRUE.
                  ELSE
                          IF ( J .LE. NBl ) EXCITE = .TRUE.
                  ENDIF
                  IF ( EXCITE ) THEN
                        K1=K1+1
                        DO K=1,M
                                EXC(K1,K) = 1
                        ENDDO
                        EXC(K1,I) = 0
                        EXC(K1,J) = 1
                        DO K=1,M+MOD(Ne,2)
                                EXC(K1,K+NBl) = 1
                        ENDDO
                        !print*,K1,EXC(K1,:)
                        !print*,'#####################################'
                ENDIF
           ENDIF
         ENDDO
       ENDDO

      DO I= M+MOD(Ne,2) - (Nel+MOD(Nel,2))/2 + 1,M + MOD(Ne,2)
        DO J=1,2*NBl
          IF ( J .GT. M .AND. ( J .LE. NBl .OR. J .GT. NBl+M+MOD(Ne,2) ) ) THEN
                EXCITE = .FALSE.
                IF ( .NOT. SPINCONSERVE ) THEN
                        EXCITE = .TRUE.
                ELSE
                        IF ( J .GT. NBl ) EXCITE = .TRUE.
                ENDIF
                IF ( EXCITE ) THEN
                        K1=K1+1
                        DO K=1,M+MOD(Ne,2)
                                EXC(K1,K+NBl) = 1
                        ENDDO
                        EXC(K1,NBl+I) = 0
                        EXC(K1,J) = 1
                        DO K=1,M
                                EXC(K1,K) = 1
                        ENDDO
                        !print*,K1,EXC(K1,:)
                        !print*,'#####################################'
                ENDIF
           ENDIF
         ENDDO
       ENDDO

       !---------------------------------------------------------------------
       ! (2) Storing all the different double excitations in the matrix EXC:
       !---------------------------------------------------------------------

       DO I=1,Ne
        DO J=I+1,Ne
          DO K=1,2*NBl
           DO L=K+1,2*NBl
            IF ( K .GT. M  .AND.  ( K .LE. NBl .OR. K .GT. NBl+M+MOD(Ne,2)) .AND. L .GT. M .AND. ( L .LE. NBl .OR. L .GT. NBl+M+MOD(Ne,2)) ) THEN
                EXCITE = .FALSE.
                II = I
                JJ = J
                IF ( I .GT. (Ne-MOD(Ne,2))/2 ) II = I - (Ne-MOD(Ne,2))/2 + NBl
                IF ( J .GT. (Ne-MOD(Ne,2))/2 ) JJ = J - (Ne-MOD(Ne,2))/2 + NBl
                IF ( .NOT. SPINCONSERVE ) THEN
                        EXCITE = exciteornot(Ne,Nel,I,J)
                ELSE
                        IF ( (II .LE. NBl .AND. K .LE. NBl ) .OR. (II .GT. NBl .AND. K .GT. NBl ) ) THEN
                                IF ( (JJ .LE. NBl .AND. L .LE. NBl ) .OR. (JJ .GT. NBl .AND. L .GT. NBl ) ) THEN
                                        EXCITE = exciteornot(Ne,Nel,I,J)
                                ENDIF
                        ENDIF
                        IF ( (II .LE. NBl .AND. L .LE. NBl ) .OR. (II .GT. NBl .AND. L .GT. NBl ) ) THEN
                                IF ( (JJ .LE. NBl .AND. K .LE. NBl ) .OR. (JJ .GT. NBl .AND. K .GT. NBl ) ) THEN
                                        EXCITE = exciteornot(Ne,Nel,I,J)
                                ENDIF
                        ENDIF
                ENDIF
                IF ( EXCITE ) THEN
                        K1 = K1+1
                        DO N=1,M
                                EXC(K1,N) = 1
                                EXC(K1,NBl+N) =1
                        ENDDO
                        EXC(K1,NBl+M+MOD(Ne,2)) = 1
                        EXC(K1,II) = 0
                        EXC(K1,JJ) = 0
                        EXC(K1,K) = 1
                        EXC(K1,L) = 1
                        !print*,K1,EXC(K1,:)
                        !print*,'#####################################'
                ENDIF
            ENDIF
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        
        !------------------------------------------------
        ! (3) Calculating the Hamiltonian matrix elements
        !------------------------------------------------
        
        !print*,'==========================================================='
        !print*,'        CONSTRUCTING THE CI-HAMILTONIAN MATRIX:            '
        !print*,'==========================================================='
        !$OMP PARALLEL SHARED(HAMILTONIAN,Ints,Cup,Cdown,h,E0,NB,NBl,NCI,EXC) &
        !$OMP & PRIVATE(K,L,A,B,Q,LD,LU,DIFFER,COEFFUP,COEFFDOWN,Pup,Pdown,PT,Jup,Jdown,JT,Kup,Kdown,KT,Fup,Fdown,F,V1,Ca,Calpha,Cb,Cbeta,TEMPMAT,J)
        ALLOCATE(PT(NB,NB),Pup(NB,NB),Pdown(NB,NB),Calpha(NB),Ca(NB),Cbeta(NB),Cb(NB),TEMPMAT(NB,NB),F(NB,NB),JT(NB,NB),KT(NB,NB),V1(2*NBl))
        ALLOCATE(COEFFUP(NB,NB),COEFFDOWN(NB,NB),Jup(NB,NB),Jdown(NB,NB),Kup(NB,NB),Kdown(NB,NB),Fup(NB,NB),Fdown(NB,NB))
        !$OMP DO 
        DO I=1,NCI
           DO J=I,NCI
            IF ( I .EQ. 1 .AND. J .EQ. 1 ) HAMILTONIAN(I,J) = E0 

            ! Calculating the number of orbitals which the slater determinant
            ! Phi_i and Ph_j differ
            
            V1 = EXC(I,:)-EXC(J,:)
            DIFFER = Int(DOT_PRODUCT(V1,V1))/2
            
            !----------------------------------
            ! See red folder p. 16 Eqn. 41-42
            !----------------------------------

            IF ( DIFFER .EQ. 0 ) THEN
                    !----------------------------------
                    ! Determining the density marices:
                    !----------------------------------
                    LU = 0
                    LD = 0
                    COEFFUP(:,:) = 0.0d0
                    COEFFDOWN(:,:) = 0.0d0
                    DO K=1,2*NBl
                        IF ( EXC(I,K) .EQ. 1 ) THEN
                                IF ( K <= NBl) THEN
                                        LU = LU+1
                                        COEFFUP(:,LU) = Cup(:,K)
                                ENDIF
                                IF ( K > NBl) THEN
                                        LD = LD+1
                                        COEFFDOWN(:,LD) = Cdown(:,K-NBl)
                                ENDIF
                        ENDIF
                    ENDDO
                    CALL makedens(COEFFUP,NB,Pup)
                    CALL makedens(COEFFDOWN,NB,Pdown)
                    PT = Pup+Pdown
                    !--------------------------------------------
                    ! The density matrices have been determined
                    !--------------------------------------------
                    CALL getJ(Pup,NB,Ints,Jup)
                    CALL getJ(Pdown,NB,Ints,Jdown)
                    CALL getK(Pup,NB,Ints,Kup)
                    CALL getK(Pdown,NB,Ints,Kdown)
                    Fup = h + Jup - Kup + Jdown
                    Fdown = h + Jdown - Kdown + Jup
                    HAMILTONIAN(I,J) = 0.50d0*( SUM(h*PT) + SUM(Fup*Pup) + SUM(Fdown*Pdown)  )
             ENDIF

             ! Brillouins theorem is applied here:
             IF ( DIFFER .EQ. 1 .AND. ( I .EQ. 1 .OR. J .EQ. 1) ) HAMILTONIAN(I,J)  = 0.0d0
                
             IF ( DIFFER .EQ. 1 .AND. I .NE. 1 .AND. J .NE. 1 ) THEN
                     L = 1
                     A(:) = 0
                     ! determining which orbitals differ: a[0] <--> a[1]
                     DO K=1,2*NBl
                        IF ( V1(K) .NE.  0 ) THEN
                                A(L) =  K
                                L = L+1
                        ENDIF
                     ENDDO

                     !-----------------------------------------------
                     ! See red folder p. 13-14, Eqn 36,37:
                     !-----------------------------------------------
                     
                     ! if the differing orbitals have oposite spin, i.e if S(a[0]) != S(a[1]):
                     IF (( A(1) .LE. NBl .AND.  A(2) .GT. NBl ) .OR. ( A(1) .GT. NBl .AND. A(2) .LE. NBl )) HAMILTONIAN(I,J) = 0.0d0
                     
                     ! if the spins differing spins are parallell, i.e if S(a[0]) = S(a[1]):
                     IF (( A(1) .LE. NBl .AND.  A(2) .LE. NBl ) .OR. ( A(1) .GT. NBl .AND. A(2) .GT. NBl )) THEN
                        !----------------------------------
                        ! Determining the density marices:
                        !----------------------------------
                        LU = 0
                        LD = 0
                        COEFFUP(:,:) = 0.0d0
                        COEFFDOWN(:,:) = 0.0d0
                        DO K=1,2*NBl
                                IF ( EXC(I,K) .EQ. 1 .AND. K .NE. A(1) .AND. K .NE. A(2) ) THEN
                                        IF ( K <= NBl ) THEN
                                                LU = LU+1
                                                COEFFUP(:,LU) = Cup(:,K)
                                        ENDIF
                                        IF ( K > NBl ) THEN
                                                LD = LD+1
                                                COEFFDOWN(:,LD) = Cdown(:,K-NBl)
                                        ENDIF
                                ENDIF
                        ENDDO
                        CALL makedens(COEFFUP,NB,Pup)
                        CALL makedens(COEFFDOWN,NB,Pdown)
                        PT = Pup+Pdown
                        !--------------------------------------------
                        ! The density matrices have been determined
                        !--------------------------------------------
                        CALL getJ(PT,NB,Ints,JT)
                        CALL getK(Pup,NB,Ints,Kup)
                        CALL getK(Pdown,NB,Ints,Kdown)

                        IF ( A(1) .LE. NBl ) THEN
                                F = h + JT - Kup
                                HAMILTONIAN(I,J) = DOT_PRODUCT(MATMUL(Cup(:,A(1)),F),Cup(:,A(2)))
                        ENDIF

                        IF ( A(1) .GT. NBl ) THEN
                                F = h + JT - Kdown
                                HAMILTONIAN(I,J) = DOT_PRODUCT(MATMUL(Cdown(:,A(1)-NBl),F),Cdown(:,A(2)-NBl))
                        ENDIF
                     ENDIF
             ENDIF

             !------------------------------------------------------------------
             !Here we calculate matrix elements between slater determinants that
             !differ with 2 orbitals See red folder p. 15 Eqn. 39,40
             !------------------------------------------------------------------
             
             IF ( DIFFER .EQ. 2 ) THEN
                     L = 1
                     Q = 1
                     A(:) = 0
                     B(:) = 0
                     DO K=1,2*NBl
                        IF  ( V1(K) .NE. 0 .AND. EXC(I,K) .NE. 0 ) THEN
                                A(L) = K
                                L = L+1
                        ENDIF
                        IF  ( V1(K) .NE. 0 .AND. EXC(J,K) .NE. 0 ) THEN
                                B(Q) = K
                                Q = Q+1
                        ENDIF
                     ENDDO
                     ! if the spins of the differing orbitals are parallell, i.e S(a(1)) = S(b(1))
                     ! and S(a(2)) = S(b(2)):
                     IF ( ( A(1) .LE. NBl .AND. B(1) .LE. NBl ) .OR. ( A(1) .GT. NBl .AND. B(1) .GT. NBl ) ) THEN
                             IF ( ( A(2) .LE. NBl .AND. B(2) .LE. NBl ) .OR. ( A(2) .GT. NBl .AND. B(2) .GT. NBl ) ) THEN
                                     JT = 0.0d0
                                     IF ( A(1) .LE. NBl ) THEN
                                             Calpha = Cup(:,B(1))
                                             Ca = Cup(:,A(1))
                                     ENDIF
                                     IF ( A(1) .GT. NBl ) THEN
                                             Calpha = Cdown(:,B(1)-NBl)
                                             Ca = Cdown(:,A(1)-NBl)
                                     ENDIF
                                     IF ( A(2) .LE. NBl ) THEN
                                             Cbeta = Cup(:,B(2))
                                             Cb = Cup(:,A(2))
                                     ENDIF
                                     IF ( A(2) .GT. NBl ) THEN
                                             Cbeta = Cdown(:,B(2)-NBl)
                                             Cb = Cdown(:,A(2)-NBl)
                                     ENDIF
                                     DO M=1,NB
                                        DO N=M,NB
                                                TEMPMAT(:,:) = Ints(M,N,:,:)
                                                JT(M,N) = DOT_PRODUCT(Cbeta,MATMUL(TEMPMAT,Cb))
                                                IF ( M .NE. N ) JT(N,M) = JT(M,N)
                                        ENDDO
                                     ENDDO
                                     HAMILTONIAN(I,J) = DOT_PRODUCT(Calpha,MATMUL(JT,Ca))
                               ENDIF
                      ENDIF
                      ! if the spinns of the differing orbitals are parallell, i.e S(a[0]) == S(b[1])
                      ! and S(a[1]) = S(b[0])
                      IF ( ( A(1) .LE. NBl .AND. B(2) .LE. NBl ) .OR. ( A(1) .GT. NBl .AND. B(2) .GT. NBl ) ) THEN
                             IF ( ( A(2) .LE. NBl .AND. B(1) .LE. NBl ) .OR. ( A(2) .GT. NBl .AND. B(1) .GT. NBl ) ) THEN
                                     KT = 0.0d0
                                     IF ( A(1) .LE. NBl ) THEN
                                             Cbeta = Cup(:,B(2))
                                             Ca = Cup(:,A(1))
                                     ENDIF
                                     IF ( A(1) .GT. NBl ) THEN
                                             Cbeta = Cdown(:,B(2)-NBl)
                                             Ca = Cdown(:,A(1)-NBl)
                                     ENDIF
                                     IF ( A(2) .LE. NBl ) THEN
                                             Calpha = Cup(:,B(1))
                                             Cb = Cup(:,A(2))
                                     ENDIF
                                     IF ( A(2) .GT. NBl ) THEN
                                             Calpha = Cdown(:,B(1)-NBl)
                                             Cb = Cdown(:,A(2)-NBl)
                                     ENDIF
                                     DO M=1,NB
                                        DO N=M,NB
                                                TEMPMAT(:,:) = Ints(:,:,M,N)
                                                KT(M,N) = DOT_PRODUCT(Calpha,MATMUL(TEMPMAT,Cb))
                                                IF ( M .NE. N ) KT(N,M) =KT(M,N)
                                        ENDDO
                                     ENDDO
                                     HAMILTONIAN(I,J) = HAMILTONIAN(I,J) - DOT_PRODUCT(Cbeta,MATMUL(KT,Ca))
                              ENDIF
                      ENDIF  
             ENDIF
             IF ( I .NE. J ) HAMILTONIAN(J,I) = HAMILTONIAN(I,J)
          ENDDO
        ENDDO
        !$OMP END DO
        !$OMP END PARALLEL
        
        ! Diagonalization of Hamiltonian
        !print*,'==========================================================='
        !print*,'            DIAGONALIZING THE HAMILTONIAN:                 '
        !print*,'==========================================================='
        CALL diagh( HAMILTONIAN,NCI,EIGENVAL,EIGENVECT)
        print*,' '
        print*,' '
        print*,'==========================================================='
        print*,' Results from the Configuration Interaction calculation:   '
        print*,'==========================================================='
        print*,'                (Singles and doubles only)                 '
        print*,'==========================================================='
        print*,' '
        print*,' '
        !-----------------------------
        ! HERE THE OUTPUT IS GENERATED
        !-----------------------------
        WRITE(*,'(A27,F30.20,A3)')'      Correlation Energy =',EIGENVAL(1)-E0,' au'
        WRITE(*,'(A27,F30.20,A3)')'     Hartree-Fock Energy =',E0,' au'
        IF ( nuce .NE. 0) WRITE(*,'(A28,F30.20,A3)')' Nuclear repulsion Energy =',nuce,' au'
        IF ( nuce .NE. 0) WRITE(*,'(A28,F30.20,A3)')'             Total Energy =',EIGENVAL(1)+nuce,' au'
        IF ( nuce .EQ. 0) WRITE(*,'(A28,F30.20,A3)')'             Total Energy =',EIGENVAL(1),' au'
        print*,' '
        print*,' '
        1 FORMAT(100(F30.20))
        
        ECISD = EIGENVAL(1)+nuce

        IF ( WRITECICOEF ) OPEN(7,FILE='CIEXPANSIONCOEFF.dat',ACTION='WRITE')
        
        ! END OF FORMAT CONSTRUCTION
        IF ( WRITECICOEF ) THEN
                DO I=1,NCI
                        WRITE(7,FMT=1)EIGENVECT(I,:)
                ENDDO
                close(7)
        ENDIF
        OPEN(8,FILE='CIEIGENVALUES.dat',ACTION='WRITE')
        WRITE(8,FMT=1)EIGENVAL
        close(8)
       
        END SUBROUTINE CISD
