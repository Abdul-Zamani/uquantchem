SUBROUTINE moleculardynamicssoft(gradS,gradT,gradV,gradIntsv,S,H0,Intsv,NB,NRED,Ne,MULTIPLICITY,BSURHF,nucE,Tol,MIX,DIISORD,DIISSTART,NATOMS,NTIMESTEPS,DT,BAS,ATOMS,APPROXEE,CORRLEVEL,PRYSR,PRYSW, &
           & WRITEONFLY,MOVIE,SAMPLERATE,TEMPERATURE,ZEROSCF,XLBOMD,kappa,alpha,CN,PULAY,FIXNSCF,DORDER,NTOTALQUAD,Q1,Q2,Q3,LORDER,CGORDER,LQ,CGQ,EETOL, &
           & CNSTART,alphastart,kappastart,ETEMP,OMEGA,EDIR,NEPERIOD,EPROFILE,EFIELDMAX,ADEF)
        USE datatypemodule
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NTOTALQUAD,LORDER,CGORDER,Q1(NTOTALQUAD),Q2(NTOTALQUAD),Q3(NTOTALQUAD)
        INTEGER, INTENT(IN) :: NB,Ne,multiplicity,NATOMS,NTIMESTEPS,DIISORD,DIISSTART,SAMPLERATE,PULAY,FIXNSCF,DORDER,EDIR,NEPERIOD
        INTEGER*8, INTENT(IN) :: NRED
        LOGICAL, INTENT(IN) :: APPROXEE,WRITEONFLY,MOVIE,ZEROSCF,XLBOMD,ADEF
        CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL,EPROFILE
        DOUBLE PRECISION, INTENT(INOUT) :: S(NB,NB), H0(NB,NB), Intsv(NRED)
        TYPE(ATOM), INTENT(INOUT) :: ATOMS(NATOMS)
        TYPE(BASIS), INTENT(INOUT)  :: BAS
        DOUBLE PRECISION, INTENT(IN) :: MIX(2),Tol,PRYSR(25,25),PRYSW(25,25),kappa,alpha,CN(DORDER),EETOL,CNSTART(4),alphastart,kappastart,ETEMP,OMEGA
        DOUBLE PRECISION, INTENT(IN) :: DT,TEMPERATURE,LQ(LORDER,3),CGQ(CGORDER,2),EFIELDMAX
        DOUBLE PRECISION :: ETOT, EHFeigen(NB),EIGENVECT(NB,NB),EHFeigenup(NB),EHFeigendown(NB),Cup(NB,NB),Cdown(NB,NB),DE,EOLD,mu,ENTROPY
        DOUBLE PRECISION, INTENT(INOUT) :: gradS(NATOMS,3,NB,NB),gradT(NATOMS,3,NB,NB),gradV(NATOMS,3,NB,NB),gradIntsv(NATOMS,3,NRED),NucE
        DOUBLE PRECISION :: leng, f(3), g(3),Rn(3),force(NATOMS,3),T(NB,NB),V(NB,NB),R0(NATOMS,3),R1(NATOMS,3),R2(NATOMS,3)
        DOUBLE PRECISION :: forceold(NATOMS,3),gradold(NATOMS,3),grad(NATOMS,3),E0,E1,E2,E3,E4,const1,const2,DRM,nom,normgrad,DRMOLD,VCM(3),MTOT
        DOUBLE PRECISION :: Pup(NB,NB),Pdown(NB,NB),P(NB,NB),PNu(DORDER,NB,NB),PNd(DORDER,NB,NB),PT(NB,NB),EFIELD,TIME
        DOUBLE PRECISION :: Jup(NB,NB),Jdown(NB,NB),Kup(NB,NB),Kdown(NB,NB),Fup(NB,NB),Fdown(NB,NB),DPTENSOR(NB,NB),DPTENSORT(3,NB,NB)
        LOGICAL :: CFORCE,SCRATCH,LSEARCH,ST,RESTART,VELOINIT,ZEROSCFF,bsurhf
        INTEGER :: I,J,II,JJ,KK,NPOINTS,JSTART,JEND,MM,clock,I0,NSCF,I1,I2,FIXNSCFF
        DOUBLE PRECISION :: ENERGY(NTIMESTEPS),VEL(NTIMESTEPS,NATOMS,3),R(NTIMESTEPS,NATOMS,3),EKIN(NTIMESTEPS),EPOT(NTIMESTEPS),TEMP(NTIMESTEPS)
        INTEGER, ALLOCATABLE :: seed(:)
        REAL :: RAND
        DOUBLE PRECISION, PARAMETER :: KB = 0.000003166811524300030d0 !  Boltzmanns constant in au/K
        DOUBLE PRECISION, EXTERNAL :: exc
        DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0
        
        IF ( EPROFILE .NE. 'ST' ) THEN
                EFIELD = 0.0d0
        ELSE
                EFIELD = EFIELDMAX
        ENDIF
        ENTROPY = 0.0d0

        OPEN(100,FILE='MOLDYNENERGY.dat',ACTION='WRITE')
        
        IF ( MOVIE ) OPEN(200,FILE='MOLDYNMOVIE.xsf',ACTION='WRITE')
        
        !---------------------------------------
        ! Seeding the random number generator
        !---------------------------------------
        ! only used for selecting the directions
        ! of the velocities
        !---------------------------------------
        call random_seed(size=MM)
        ALLOCATE(seed(MM))
        CALL SYSTEM_CLOCK(COUNT=clock)
        seed = clock + 37 * (/ (i - 1, i = 1, MM) /)
        CALL RANDOM_SEED(PUT = seed )
        !-------------------------------------

        

        I = 1
        I0 = 0
        EOLD = 0.0d0
        SCRATCH = .TRUE.
        ZEROSCFF = .FALSE.
        FIXNSCFF = -1
        ST = .TRUE.
        E0 = 0.0d0
        KK = 0
        VCM = 0.0d0
        MTOT = 0.0d0
        RESTART = .FALSE.
        VELOINIT = .FALSE.

        !=================================================================
        ! Are we going to continue the calculation from a previous MD-run?
        !=================================================================
        inquire(file='MOLDYNRESTART.dat',exist=RESTART)
        IF ( RESTART ) THEN
                OPEN(300,FILE='MOLDYNRESTART.dat',STATUS='OLD',ACTION='READ')
                READ(300,'(I10)')I0
                DO J=1,NATOMS
                        READ(300,'(3(E30.20))')ATOMS(J)%R
                ENDDO
                DO J=1,NATOMS
                        READ(300,'(3(E30.20))')VEL(1,J,:)
                ENDDO
                DO J=1,NATOMS
                        READ(300,'(3(E30.20))')force(J,:)
                ENDDO
                ! reading the regular density matrices
                DO I1=1,NB
                     DO I2=1,NB
                          READ(300,'(E30.20)')Pup(I1,I2)
                     ENDDO
                ENDDO
                DO I1=1,NB
                     DO I2=1,NB
                          READ(300,'(E30.20)')Pdown(I1,I2)
                     ENDDO
                ENDDO
                IF ( XLBOMD ) THEN
                        ! Reading the auxiliary density matrices
                        DO J=1,DORDER
                                DO I1=1,NB
                                        DO I2=1,NB
                                                READ(300,'(E30.20)')PNu(J,I1,I2)
                                        ENDDO
                                ENDDO
                        ENDDO   
                        DO J=1,DORDER
                                DO I1=1,NB
                                        DO I2=1,NB
                                                READ(300,'(E30.20)')PNd(J,I1,I2)
                                        ENDDO
                                ENDDO
                        ENDDO  
                ENDIF
                CLOSE(300)
                I0 = I0-1
        ENDIF
        !=====================================================================
        ! If this is not a continuation, have initial velocities been provided?
        !=====================================================================
        IF ( .not. RESTART ) inquire(file='INITVELO.dat',exist=VELOINIT)
        
        IF ( VELOINIT ) THEN
                OPEN(400,FILE='INITVELO.dat',STATUS='OLD',ACTION='READ')
                DO J=1,NATOMS
                        READ(400,*)VEL(1,J,1),VEL(1,J,2),VEL(1,J,3)
                        VCM = VCM + VEL(1,J,:)*ATOMS(J)%M
                        MTOT = MTOT + ATOMS(J)%M
                ENDDO
                CLOSE(400)
                VCM = VCM/MTOT
        ENDIF
        
        ! Initialization of positions:
        DO J=1,NATOMS
                R(1,J,:) = ATOMS(J)%R
        ENDDO

        IF ( .not. RESTART .AND. .not. VELOINIT ) THEN
                ! Initialization of velocities:
                DO J=1,NATOMS
                        CALL RANDOM_NUMBER(RAND)
                        f(1) = RAND
                        CALL RANDOM_NUMBER(RAND)
                        f(2) = RAND
                        CALL RANDOM_NUMBER(RAND)
                        f(3) = RAND
                        IF ( MOD(J,2) .NE. 0 ) THEN
                                VEL(1,J,:) = (f/sqrt(DOT_PRODUCT(f,f)))*sqrt(6.0d0*KB*TEMPERATURE/ATOMS(J)%M)
                        ELSE
                                VEL(1,J,:) = -(VEL(1,J-1,:)/sqrt(DOT_PRODUCT(VEL(1,J-1,:),VEL(1,J-1,:))))*sqrt(6.0d0*KB*TEMPERATURE/ATOMS(J)%M)
                        ENDIF
                        ! Calculating center of mass velocity (to be subtracted later)
                        VCM = VCM + VEL(1,J,:)*ATOMS(J)%M
                        MTOT = MTOT + ATOMS(J)%M
                ENDDO
                VCM = VCM/MTOT
        ENDIF

        ! subtracting the center of mass velocity from all the atoms:
        DO J=1,NATOMS
                VEL(1,J,:)  = VEL(1,J,:) - VCM
        ENDDO
        
        DO WHILE ( I .LT. NTIMESTEPS )
                TIME = (I-1)*DT
                IF ( XLBOMD ) THEN
                        IF ( I .GT. 50 .AND. I .LE. 100 ) ZEROSCFF = .TRUE.
                        IF ( I .EQ. 101 ) ZEROSCFF = .FALSE.
                ENDIF
                !===================================================
                ! Making sure that the first 10 time steps are fully 
                ! self consistent, eaven though ZEROSCF = .TRUE.
                !====================================================
                IF ( ZEROSCF ) THEN
                        IF ( RESTART ) ZEROSCFF = .TRUE.
                        IF ( .not. RESTART .AND. I .GT. 100 ) ZEROSCFF = .TRUE.
                ENDIF
                !===================================================
                ! Making sure that the first 10 time steps are fully 
                ! self consistent, eaven though FIXNSCF > 0 
                !====================================================
                IF ( FIXNSCF .GT. 0 .AND. .not. ZEROSCF ) THEN
                        IF ( RESTART ) FIXNSCFF = FIXNSCF
                        IF ( .not. RESTART .AND. I .GT. 100 ) FIXNSCFF = FIXNSCF
                ENDIF
         !=======================================================================================================================
                
                IF ( I+I0 .EQ. 1  ) THEN
                        IF ( CORRLEVEL .EQ. 'RHF' ) THEN
                                CALL RHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol,EHFeigen,ETOT,EIGENVECT,P,MIX,DIISORD,DIISSTART,NSCF,-1,.FALSE., SCRATCH,.FALSE.)
                                EHFeigenup = EHFeigen
                                EHFeigendown = EHFeigen
                                Cup = EIGENVECT
                                Cdown = EIGENVECT
                                Pup = P/2.0d0
                                Pdown = P/2.0d0
                        ENDIF
      
                        IF ( CORRLEVEL .EQ. 'URHF' ) THEN
                                CALL URHF(S,H0,Intsv,NB,NRED,Ne,MULTIPLICITY,BSURHF,nucE,Tol,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1,.FALSE.,SCRATCH,.FALSE.,ETEMP,ENTROPY)
                        ENDIF

                        IF ( CORRLEVEL .EQ. 'PBE' .OR. CORRLEVEL .EQ. 'LDA' .OR.  CORRLEVEL .EQ. 'B3LYP' ) THEN
                               CALL DFT(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,S,gradS,H0,Intsv,NB,NRED,Ne,MULTIPLICITY,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown, &
                                & ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1,.FALSE.,SCRATCH,.FALSE.,ETEMP,mu,ENTROPY)
                        ENDIF
     
                        ! Calculating forces on atoms:
                        CALL forces(CORRLEVEL,NATOMS,Ne,NB,NRED,Cup,Cdown,Pup,Pdown,EHFeigenup,EHFeigendown,ATOMS,BAS,gradS,gradT,gradV,gradIntsv,S,H0,Intsv,PULAY, &
                        & NTOTALQUAD,Q1,Q2,Q3,LORDER,CGORDER,LQ,CGQ,force,EDIR,EPROFILE,EFIELD,OMEGA,TIME,ADEF)

                        ! In the case of XLBOMD we here initialize the auxiliary density matrices
                        IF ( XLBOMD .AND. .not. RESTART ) THEN
                                DO J=1,DORDER
                                        PNu(J,:,:) = Pup
                                        PNd(J,:,:) = Pdown
                                ENDDO
                        ENDIF
                 ENDIF
                 
                 ! Saveing the force of the current time-step:
                 forceold = force

                 !=======================================================================================================================
                 ! The velocity Verlet algorithm as presented in J.M Thijssen's
                 ! "Computational Physics" P.181. Here we update the positions:
                 !=======================================================================================================================
                  DO J=1,NATOMS
                          R(I+1,J,:) = R(I,J,:) + DT*VEL(I,J,:) + (DT**2)*force(J,:)/(2.0d0*ATOMS(J)%M)
                  ENDDO

                 ! Uppdating the atomic positions to the positions of the
                 ! current time-step:
                 DO J=1,NATOMS
                        ATOMS(J)%R = R(I+1,J,:)
                 ENDDO

                 !-----------------------------------------------------------
                 ! Calculating the nuclear-nuclear repulsion energy of the
                 ! updated nuclear configuration
                 !-----------------------------------------------------------
                 NucE = 0.0d0
                 DO II=1,NATOMS
                      DO JJ=II+1,NATOMS
                           Rn = ATOMS(II)%R - ATOMS(JJ)%R
                           NucE = NucE + ATOMS(II)%Z*ATOMS(JJ)%Z/sqrt(DOT_PRODUCT(Rn,Rn))
                      ENDDO
                 ENDDO
                 
                 !=================================================================================================
                 ! In the case of XLBOMD, we here integrate the electronic degrees of freedoom ala A.M Niklasson
                 !=================================================================================================
                 IF ( XLBOMD ) THEN
                         IF ( I .LE. 100 ) THEN
                                Pup   = 2.0d0*PNu(1,:,:) - PNu(2,:,:) + Kappastart*(Pup   - PNu(1,:,:) )
                                Pdown = 2.0d0*PNd(1,:,:) - PNd(2,:,:) + Kappastart*(Pdown - PNd(1,:,:) )
                                DO J=1,4
                                        Pup   = Pup   + alphastart*CNSTART(J)*PNu(J,:,:)
                                        Pdown = Pdown + alphastart*CNSTART(J)*PNd(J,:,:)
                                ENDDO
                        ELSE IF ( I .GT. 100 ) THEN
                                Pup   = 2.0d0*PNu(1,:,:) - PNu(2,:,:) + Kappa*(Pup   - PNu(1,:,:) )
                                Pdown = 2.0d0*PNd(1,:,:) - PNd(2,:,:) + Kappa*(Pdown - PNd(1,:,:) )
                                DO J=1,DORDER
                                        Pup   = Pup   + alpha*CN(J)*PNu(J,:,:)
                                        Pdown = Pdown + alpha*CN(J)*PNd(J,:,:)
                                ENDDO
                        ENDIF
                         IF ( CORRLEVEL .EQ. 'RHF' ) P = Pup + Pdown
                         ! Here we uppgrade the auxiliary density matrices:
                         DO J=1,DORDER-1
                                PNu(DORDER+1-J,:,:) = PNu(DORDER-J,:,:)
                                PNd(DORDER+1-J,:,:) = PNd(DORDER-J,:,:)
                         ENDDO
                         PNu(1,:,:) = Pup
                         PNd(1,:,:) = Pdown
                         !==================================================
                 ENDIF
                 !=================================================================================================

                 ! Updating the basis set:
                 CALL gettotalbasis(NATOMS,ATOMS,NB,BAS,.FALSE.)
                 ! Here the normalization of the basis functions is performed
                 CALL normalize(BAS)
                 ! Calculating the matrices in the basis of the updated positions:
                 CALL overlap(NATOMS,BAS,S,gradS)
                 CALL kinetic(NATOMS,BAS,T,gradT)
                 CALL potential(BAS,NATOMS,ATOMS,V,gradV)
                
                 H0 = T + V
                
                 ! In the case we have an electric field present:
                 IF ( ADEF ) THEN
                         IF ( EPROFILE .EQ. 'AC' ) THEN
                                 CALL dipoletensorinhom(BAS,EDIR,OMEGA,TIME,DPTENSOR)
                         ELSE
                                 CALL dipoletensor(BAS,DPTENSORT)
                                 IF ( EDIR .EQ. 1 ) DPTENSOR = DPTENSORT(2,:,:)
                                 IF ( EDIR .EQ. 2 ) DPTENSOR = DPTENSORT(3,:,:)
                                 IF ( EDIR .EQ. 3 ) DPTENSOR = DPTENSORT(1,:,:)
                         ENDIF

                         IF ( EPROFILE .EQ. 'AC' .OR. EPROFILE .EQ. 'HO' ) THEN
                                IF ( TIME .GE. 0.0d0 .AND. TIME .LE. 2.0d0*pi/OMEGA ) THEN
                                        EFIELD = EFIELDMAX*(OMEGA*TIME/(2.0d0*pi))
                                ENDIF

                                IF ( TIME .GT. 2.0d0*pi/OMEGA .AND. TIME .LE. (NEPERIOD + 1.0d0)*(2.0d0*pi/OMEGA) ) THEN
                                        EFIELD = EFIELDMAX
                                ENDIF

                                IF ( TIME .GT. (NEPERIOD + 1.0d0)*(2.0d0*pi/OMEGA) .AND.  TIME .LT. (NEPERIOD + 2.0d0)*(2.0d0*pi/OMEGA) ) THEN
                                        EFIELD = EFIELDMAX*( (NEPERIOD + 2.0d0) - OMEGA*TIME/(2.0d0*pi) )
                                ENDIF

                                IF ( TIME .GT. (NEPERIOD + 2.0d0)*(2.0d0*pi/OMEGA) ) EFIELD = 0.0d0

                                IF ( EPROFILE .EQ. 'AC' ) THEN
                                        DPTENSOR = DPTENSOR*EFIELD
                                ELSE IF ( EPROFILE .EQ. 'HO' ) THEN
                                        EFIELD = EFIELD*sin(OMEGA*TIME)
                                        DPTENSOR = DPTENSOR*EFIELD
                                ENDIF
                        ELSE IF ( EPROFILE .EQ. 'ST' ) THEN
                                EFIELD = EFIELDMAX
                                DPTENSOR = DPTENSOR*EFIELD
                        ENDIF

                        H0 = H0 - DPTENSOR
                 ENDIF

                 ! Calculation of electron repulsion tensor (ij|kl) in the basis of the updated positions:
                 CALL eeints(NATOMS,NATOMS,BAS,Intsv,gradIntsv,NRED,NRED,PRYSR,PRYSW,APPROXEE,EETOL,.TRUE.,Pup+Pdown)
                
                 !=============================================================================
                 ! Here we have that: P --> RHF/URHF --> P = D  ( or Pup = Dup, Pdown = Ddown )
                 !=============================================================================
                 IF ( CORRLEVEL .EQ. 'RHF' ) THEN
                         CALL RHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol,EHFeigen,ETOT,EIGENVECT,P,MIX,DIISORD,DIISSTART,NSCF,FIXNSCFF,.FALSE., SCRATCH,ZEROSCFF)
                         EHFeigenup = EHFeigen
                         EHFeigendown = EHFeigen
                         Cup = EIGENVECT
                         Cdown = EIGENVECT
                         Pup = P/2.0d0
                         Pdown = P/2.0d0
                 ENDIF
      
                 IF ( CORRLEVEL .EQ. 'URHF' ) THEN
                        CALL URHF(S,H0,Intsv,NB,NRED,Ne,MULTIPLICITY,BSURHF,nucE,Tol,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,FIXNSCFF,.FALSE.,SCRATCH,ZEROSCFF,ETEMP,ENTROPY)
                 ENDIF
                
                 IF ( CORRLEVEL .EQ. 'PBE' .OR. CORRLEVEL .EQ. 'LDA' .OR.  CORRLEVEL .EQ. 'B3LYP' ) THEN
                        CALL DFT(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,S,gradS,H0,Intsv,NB,NRED,Ne,MULTIPLICITY,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown, &
                                & ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,FIXNSCFF,.FALSE.,SCRATCH,ZEROSCFF,ETEMP,mu,ENTROPY)

                 ENDIF
                 
                 ! Calculating forces on atoms:
                 CALL forces(CORRLEVEL,NATOMS,Ne,NB,NRED,Cup,Cdown,Pup,Pdown,EHFeigenup,EHFeigendown,ATOMS,BAS,gradS,gradT,gradV,gradIntsv,S,H0,Intsv,PULAY,& 
                 & NTOTALQUAD,Q1,Q2,Q3,LORDER,CGORDER,LQ,CGQ,force,EDIR,EPROFILE,EFIELD,OMEGA,TIME,ADEF)
                 
                 !==========================================================================
                 ! Here we calculate the potential energy using the density
                 ! matrix D ( or Dup and Ddown) that was previously produced by RHF or URHF:
                 ! ETOT = trace(F[D]*D) - trace( (J[D] - 0.50d0*K[D])*D )
                 !==========================================================================
                 IF ( XLBOMD ) THEN
                         PT = Pup + Pdown
                         CALL getJv(Pup,NB,NRED,Intsv,Jup)
                         CALL getJv(Pdown,NB,NRED,Intsv,Jdown)
                         CALL getKv(Pup,NB,NRED,Intsv,Kup)
                         CALL getKv(Pdown,NB,NRED,Intsv,Kdown)
                         IF ( CORRLEVEL .EQ. 'RHF' .OR. CORRLEVEL .EQ. 'URHF' ) THEN
                                Fup   = H0 + Jdown - Kup   + Jup
                                Fdown = H0 + Jdown - Kdown + Jup
                                ETOT = 0.50d0*(SUM(H0*PT) + SUM(Fup*Pup) + SUM(Fdown*Pdown) ) + nucE
                        ELSE IF  ( CORRLEVEL .EQ. 'PBE' .OR. CORRLEVEL .EQ.  'LDA' .OR.  CORRLEVEL .EQ. 'B3LYP' ) THEN
                                Fup   = H0 +  Jdown + Jup
                                Fdown = H0 +  Jdown + Jup
                                IF ( CORRLEVEL .EQ. 'B3LYP'  ) THEN
                                        Fup = Fup - 0.20d0*Kup
                                        Fdown = Fdown - 0.20d0*Kdown
                                ENDIF
                                ETOT = 0.50d0*(SUM(H0*PT) + SUM(Fup*Pup) + SUM(Fdown*Pdown) ) + &
                                & exc(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,Pup,Pdown,gradS,LORDER,CGORDER,LQ,CGQ) + nucE
                        ENDIF
                 ENDIF

                 
                 !=======================================================================================================================
                 ! The velocity Verlet algorithm as presented in J.M Thijssen's
                 ! "Computational Physics" P.181. Here we update the velocities:
                 !=======================================================================================================================
                  DO J=1,NATOMS
                          VEL(I+1,J,:) = VEL(I,J,:) + (DT/2.0d0)*(force(J,:) + forceold(J,:) )/ATOMS(J)%M
                  ENDDO
                 
                 !===============================================
                 ! Calculating the diferent energy contributions:
                 !===============================================

                 ! The kinetic energy:
                 EKIN(I+1) = 0.0d0
                 DO J=1,NATOMS
                        EKIN(I+1) = EKIN(I+1) + (ATOMS(J)%M/2.0d0)*DOT_PRODUCT(VEL(I+1,J,:),VEL(I+1,J,:))
                 ENDDO

                 ! The potential energy
                 EPOT(I+1) = ETOT
                 
                 ! The total energy
                 ENERGY(I+1) = EPOT(I+1) + EKIN(I+1)

                 !=================================================
                 ! Calculating the temperature:
                 !=================================================

                 TEMP(I+1) = 2.0d0*EKIN(I+1)/(3.0d0*NATOMS*KB)
                
                 IF ( I .EQ. 1 ) THEN
                         IF ( ETEMP .LT. 0.0d0 ) THEN
                          WRITE(100,'(A138)')'        ITER           E [au]                        EK [au]                       EP [au]                        T [K]              NSCF '
                          IF (WRITEONFLY) WRITE(*,'(A138)')'        ITER           E [au]                        EK [au]                       EP [au]                        T [K]              NSCF '
                         ELSE
                          WRITE(100,'(A12,A17,A30,A30,A30,A30,A18)')'ITER','E [au]','F [au]','EK [au]','EP [au]','T [K]','NSCF'
                          IF (WRITEONFLY) WRITE(*,'(A12,A17,A30,A30,A30,A30,A18)')'ITER','E [au]','F [au]','EK [au]','EP [au]','T [K]','NSCF'
                         ENDIF
                         IF ( MOVIE ) THEN
                                 WRITE(200,'(A9,I10)')'ANIMSTEPS',INT(NTIMESTEPS*1.0/(SAMPLERATE*1.0))-1
                                 IF (WRITEONFLY) WRITE(220,'(A9,I10)')'ANIMSTEPS',INT(NTIMESTEPS*1.0/(SAMPLERATE*1.0))-1
                         ENDIF
                 ENDIF 

                 ! increasing the time index one step
                 I = I + 1
                 
                 !=============================================================
                 ! Enabaling the use of the density matrix of the previous step 
                 ! to be used as a starting guess for the next time step
                 !=============================================================
                 SCRATCH = .FALSE.

                 !==============================
                 ! Saving the energies to file:
                 !==============================
                 IF (WRITEONFLY)  THEN
                         IF ( ETEMP .LT. 0.0d0 ) THEN
                                WRITE(100,'(I10,E30.20,E30.20,E30.20,E30.20,I6)')I+I0,ENERGY(I),EKIN(I),EPOT(I),TEMP(I),NSCF
                                WRITE(*,'(I10,E30.20,E30.20,E30.20,E30.20,I6)'  )I+I0,ENERGY(I),EKIN(I),EPOT(I),TEMP(I),NSCF
                         ELSE
                                WRITE(100,'(I10,E30.20,E30.20,E30.20,E30.20,E30.20,I6)')I+I0,ENERGY(I),ENERGY(I)-ETEMP*ENTROPY,EKIN(I),EPOT(I),TEMP(I),NSCF
                                WRITE(*,'(I10,E30.20,E30.20,E30.20,E30.20,E30.20,I6)'  )I+I0,ENERGY(I),ENERGY(I)-ETEMP*ENTROPY,EKIN(I),EPOT(I),TEMP(I),NSCF
                        ENDIF
                 ENDIF    
                
              
                 !==================================================
                 ! saving a xsf file used by xcrysden to make movies
                 !==================================================
                 IF ( MOVIE .AND. MOD(I,SAMPLERATE) .EQ. 0 ) THEN
                        KK = KK + 1
                        IF (WRITEONFLY) THEN
                                WRITE(200,'(A5,I10)')'ATOMS',KK
                                DO J=1,NATOMS
                                        f = R(I,J,:)*0.52917720859 ! positions saved in angstroem
                                        g = VEL(I,J,:)*0.52917720859 ! velocities in angstroem/au
                                        WRITE(200,'(I4,6(E30.20))')ATOMS(J)%Z,f,g
                                ENDDO
                       ENDIF
                 ENDIF
                 !===============================
                 ! saving the restart file 
                 !===============================
                 IF ( MOD(I,SAMPLERATE) .EQ. 0 ) THEN
                        OPEN(300,FILE='MOLDYNRESTART.dat',ACTION='WRITE')
                        WRITE(300,'(I10)')I+I0
                        DO J=1,NATOMS
                                WRITE(300,'(3(E30.20))')ATOMS(J)%R
                        ENDDO
                        DO J=1,NATOMS
                                WRITE(300,'(3(E30.20))')VEL(I,J,:)
                        ENDDO
                        DO J=1,NATOMS
                                WRITE(300,'(3(E30.20))')force(J,:)
                        ENDDO
                        ! saving the regular density matrices
                        DO I1=1,NB
                              DO I2=1,NB
                                 WRITE(300,'(E30.20)')Pup(I1,I2)
                              ENDDO
                        ENDDO
                        DO I1=1,NB
                               DO I2=1,NB
                                    WRITE(300,'(E30.20)')Pdown(I1,I2)
                               ENDDO
                        ENDDO
                        IF ( XLBOMD ) THEN
                                ! Saving the auxiliary density matrices
                                DO J=1,DORDER
                                        DO I1=1,NB
                                                DO I2=1,NB
                                                        WRITE(300,'(E30.20)')PNu(J,I1,I2)
                                                ENDDO
                                         ENDDO
                                ENDDO   
                                DO J=1,DORDER
                                        DO I1=1,NB
                                                DO I2=1,NB
                                                        WRITE(300,'(E30.20)')PNd(J,I1,I2)
                                                ENDDO
                                         ENDDO
                                ENDDO  
                        ENDIF

                        CLOSE(300)
                ENDIF
      ENDDO
      
      !===============================
      ! saving the restart file
      !===============================
      OPEN(300,FILE='MOLDYNRESTART.dat',ACTION='WRITE')
      WRITE(300,'(I10)')I+I0
      DO J=1,NATOMS
              WRITE(300,'(3(E30.20))')ATOMS(J)%R
      ENDDO
      DO J=1,NATOMS
              WRITE(300,'(3(E30.20))')VEL(I,J,:)
      ENDDO
      DO J=1,NATOMS
              WRITE(300,'(3(E30.20))')force(J,:)
      ENDDO
      ! saving the regular density matrices
      DO I1=1,NB
           DO I2=1,NB
                WRITE(300,'(E30.20)')Pup(I1,I2)
         ENDDO
      ENDDO
      DO I1=1,NB
           DO I2=1,NB
                WRITE(300,'(E30.20)')Pdown(I1,I2)
           ENDDO
      ENDDO
      IF ( XLBOMD ) THEN
             ! Saving the auxiliary density matrices
             DO J=1,DORDER
                  DO I1=1,NB
                       DO I2=1,NB
                            WRITE(300,'(E30.20)')PNu(J,I1,I2)
                       ENDDO
                  ENDDO
             ENDDO   
             DO J=1,DORDER
                 DO I1=1,NB
                      DO I2=1,NB
                           WRITE(300,'(E30.20)')PNd(J,I1,I2)
                      ENDDO
                 ENDDO
             ENDDO  
      ENDIF
      CLOSE(300)
      !=================================================

      !===================================================
      ! If specified not to save on-the-fly, we save after
      ! the calculation has been finished 
      !===================================================
      KK = 0
      IF ( .not. WRITEONFLY) THEN
                DO J=2,NTIMESTEPS
                        WRITE(100,'(I10,E30.20,E30.20,E30.20,E30.20)')J+I0,ENERGY(J),EKIN(J),EPOT(J),TEMP(J)
                        IF ( MOVIE .AND. MOD(J,SAMPLERATE) .EQ. 0 ) THEN
                                KK = KK+1
                                WRITE(200,'(A5,I10)')'ATOMS',KK
                                DO JJ=1,NATOMS
                                       f = R(J,JJ,:)*0.52917720859 ! positions saved in angstroem
                                       g = VEL(J,JJ,:)*0.52917720859 ! velocities in angstroem/au
                                       WRITE(200,'(I4,6(E30.20))')ATOMS(JJ)%Z,f,g
                                ENDDO
                        ENDIF
               ENDDO
      ENDIF
      close(100)
      IF ( MOVIE ) close(200)
      STOP
      
      END SUBROUTINE moleculardynamicssoft
