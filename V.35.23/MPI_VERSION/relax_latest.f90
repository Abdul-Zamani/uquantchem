SUBROUTINE relax(gradS,gradT,gradV,gradIntsvR,S,H0,IntsvR,NB,NRED,Ne,nucE,Tol,FTol,MIX,DIISORD,DIISSTART,NATOMS,NSTEPS,DR,BAS,ATOMS,APPROXEE,CORRLEVEL,PRYSR,PRYSW, &
& IND1,IND2,IND3,IND4,Istart,Iend,numprocessors,id)
        USE datatypemodule
        IMPLICIT NONE
        INCLUDE "mpif.h"
        INTEGER, INTENT(IN) :: NB,NRED,Ne,NATOMS,NSTEPS,DIISORD,DIISSTART,Istart,Iend,numprocessors,id
        INTEGER, INTENT(IN) :: IND1(Istart:Iend),IND2(Istart:Iend),IND3(Istart:Iend),IND4(Istart:Iend)
        LOGICAL, INTENT(IN) :: APPROXEE
        CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL
        DOUBLE PRECISION, INTENT(INOUT) :: S(NB,NB), H0(NB,NB), IntsvR(Istart:Iend)
        TYPE(ATOM), INTENT(INOUT) :: ATOMS(NATOMS)
        TYPE(BASIS), INTENT(INOUT)  :: BAS
        DOUBLE PRECISION, INTENT(IN) :: MIX(2),Tol,PRYSR(25,25),PRYSW(25,25),FTol
        DOUBLE PRECISION, INTENT(INOUT)  :: DR
        DOUBLE PRECISION :: ETOT, EHFeigen(NB),EIGENVECT(NB,NB),EHFeigenup(NB),EHFeigendown(NB),Cup(NB,NB),Cdown(NB,NB),DE,EOLD
        DOUBLE PRECISION, INTENT(INOUT) :: gradS(NATOMS,3,NB,NB),gradT(NATOMS,3,NB,NB),gradV(NATOMS,3,NB,NB),gradIntsvR(NATOMS,3,Istart:Iend),NucE
        DOUBLE PRECISION :: leng, f(3),g(3),Rn(3),force(NATOMS,3),T(NB,NB),V(NB,NB),DR0,DR1,DR2,R0(NATOMS,3),R1(NATOMS,3),R2(NATOMS,3)
        DOUBLE PRECISION :: forceold(NATOMS,3),gradold(NATOMS,3),grad(NATOMS,3),E0,E1,E2,E3,E4,const1,const2,DRM,DRMO,Ebis,Eprim,Etris,nom,normgrad
        LOGICAL :: CFORCE,SCRATCH,LSEARCH,ST
        INTEGER :: I,J,II,JJ,KK,NLITER,JSTART,JEND

        CFORCE = .TRUE.
        I = 0
        DE = 2*Tol
        leng = 20*FTol
        EOLD = 0.0d0
        SCRATCH = .TRUE.
        NLITER = 10
        ST = .TRUE.

        DO WHILE ( DABS(leng) .GT. FTol .AND. I .LT. NSTEPS ) 

                !Switching to conjugate gradient:
                IF ( DABS(leng) .LT. FTol*10.0d0 ) ST = .FALSE.
     
        !=======================================================================================================================
                
                DO J=1,NATOMS
                        R0(J,:) = ATOMS(J)%R
                ENDDO
  
                IF ( CORRLEVEL .EQ. 'RHF' ) THEN
                        CALL RHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol/100,EHFeigen,ETOT,EIGENVECT,MIX,DIISORD,DIISSTART,numprocessors,id,.FALSE.,SCRATCH)
                        EHFeigenup = EHFeigen
                        EHFeigendown = EHFeigen
                        Cup = EIGENVECT
                        Cdown = EIGENVECT
                ENDIF
                
                IF ( CORRLEVEL .EQ. 'URHF' ) THEN
                        CALL URHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol/100,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,MIX,DIISORD,DIISSTART,numprocessors,id,.FALSE.,SCRATCH)
                ENDIF
                
                ! Calculating forces on atoms:
                CALL forces(IND1,IND2,IND3,IND4,Istart,Iend,NATOMS,Ne,NB,NRED,Cup,Cdown,EHFeigenup,EHFeigendown,ATOMS,gradS,gradT,gradV,gradIntsvR,force,numprocessors,id)
        !=======================================================================================================================

                E0 = ETOT

                IF ( I .EQ. 0 ) THEN
                        grad = -force
                ELSE
                        const1 = 0.0d0
                        const2 = 0.0d0
                        DO J=1,NATOMS
                                const1 = const1 + DOT_PRODUCT(force(J,:)-forceold(J,:),force(J,:))
                                const2 = const2 + DOT_PRODUCT(forceold(J,:),forceold(J,:))
                        ENDDO
                        ! Uppdating the gradient according to Polack and Ribiere
                        DO J=1,NATOMS
                                IF ( ST ) THEN
                                        ! steepest deecent
                                        grad(J,:) = force(J,:) 
                                ELSE
                                        ! Conjugate gradient
                                        grad(J,:) = force(J,:) + (const1/const2)*gradold(J,:)
                                ENDIF
                        ENDDO
                ENDIF
                ! Calculating the norm of the gradient:
                normgrad = 0.0d0
                DO J=1,NATOMS
                        normgrad = normgrad + DOT_PRODUCT(grad(J,:),grad(J,:))
                ENDDO
                normgrad = (1.0d0/(1.0d0*NATOMS))*sqrt(normgrad)

                forceold = force
                gradold  = grad

        !=========================================================================================================================================
        ! Here we search for the energy minimum along the gradient grad.
        !=========================================================================================================================================
                DRMO    = 10.0d0
                Eprim   = 10.0d0
                DRM     = 0.0d0
                KK      = 0
                LSEARCH = .TRUE.
                DO WHILE ( LSEARCH )
                        IF ( DABS(Eprim) .LT. DR**2 .OR. KK .EQ. NLITER ) LSEARCH = .FALSE.

                        JSTART = 1
                        JEND   = 5

                        IF ( KK .EQ. 0 ) JSTART = 2

                        IF ( .not. LSEARCH ) THEN
                                JSTART = 1
                                JEND   = 1
                                DRM = DABS(DRM) ! After all the step should be in the direction of the gradient
                                IF ( DR .GE. DRM/10.0d0) DR = DRM/10.0d0
                        ENDIF

                        DO J=JSTART,JEND
                                ! Calculating new positions on atoms:
                                DO JJ=1,NATOMS
                                        g = grad(JJ,:)/normgrad
                                        IF ( J .EQ. 1 ) ATOMS(JJ)%R = R0(JJ,:) + g*DRM
                                        IF ( J .EQ. 2 ) ATOMS(JJ)%R = R0(JJ,:) + g*(DRM-DR)
                                        IF ( J .EQ. 3 ) ATOMS(JJ)%R = R0(JJ,:) + g*(DRM+DR)
                                        IF ( J .EQ. 4 ) ATOMS(JJ)%R = R0(JJ,:) + g*(DRM-2*DR)
                                        IF ( J .EQ. 5 ) ATOMS(JJ)%R = R0(JJ,:) + g*(DRM+2*DR)
                                ENDDO


                                ! Calculating the nuclear-nuclear repulsion energy of the
                                ! updated nuclear configuration
                                NucE = 0.0d0
                                DO II=1,NATOMS
                                        DO JJ=II+1,NATOMS
                                                Rn = ATOMS(II)%R - ATOMS(JJ)%R
                                                NucE = NucE + ATOMS(II)%Z*ATOMS(JJ)%Z/sqrt(DOT_PRODUCT(Rn,Rn))
                                        ENDDO
                                ENDDO

                                ! Updating the basis set:
                                CALL gettotalbasis(NATOMS,ATOMS,NB,BAS,.FALSE.)

                                ! Here the normalization of the basis functions is performed
                                CALL normalize(BAS)
                
                                ! Calculating the matrices in the basis of the updated positions:
                                CALL overlap(NATOMS,BAS,S,gradS)
                                CALL kinetic(NATOMS,BAS,T,gradT)
                                CALL potential(BAS,NATOMS,ATOMS,V,gradV)
                
                                H0 = T + V

                                ! Calculation of electron repulsion tensor (ij|kl) in the basis of the updated positions:
                                IF ( LSEARCH ) THEN
                                        CALL eeints(NATOMS,NATOMS,BAS,IntsvR,gradIntsvR,IND1,IND2,IND3,IND4,NRED,Istart,Iend,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,Tol/100,.FALSE.,id)
                                ELSE
                                        CALL eeints(NATOMS,NATOMS,BAS,IntsvR,gradIntsvR,IND1,IND2,IND3,IND4,NRED,Istart,Iend,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,Tol/100,.TRUE.,id)
                                ENDIF

                                IF ( CORRLEVEL .EQ. 'RHF' .AND. LSEARCH ) THEN
                                        CALL RHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol/100,EHFeigen,ETOT,EIGENVECT,MIX,DIISORD,DIISSTART,numprocessors,id,.FALSE.,SCRATCH)
                                        EHFeigenup = EHFeigen
                                        EHFeigendown = EHFeigen
                                        Cup = EIGENVECT
                                        Cdown = EIGENVECT
                                ENDIF
                
                                IF ( CORRLEVEL .EQ. 'URHF' .AND. LSEARCH ) THEN
                                        CALL URHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol/100,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,MIX,DIISORD,DIISSTART,& 
                                        & numprocessors,id,.FALSE.,SCRATCH)
                                ENDIF
                        
                                IF ( J .EQ. 1 ) E0 = ETOT
                                IF ( J .EQ. 2 ) E1 = ETOT
                                IF ( J .EQ. 3 ) E2 = ETOT
                                IF ( J .EQ. 4 ) E3 = ETOT
                                IF ( J .EQ. 5 ) E4 = ETOT
                       
                                ! Calculating the new position in the line-searc for the minimum: 
                                IF ( J .EQ. 5 ) THEN
                                       ! Here we use the 5-point finite difference
                                       ! formulas presented in Koonin's book
                                       ! Computaional Physics. page 6, table 1.2
                                       Eprim = (1.0d0/12.0d0)*(E3 - 8*E1 + 8*E2 - E4 )/DR
                                       Ebis  = (1.0d0/12.0d0)*(-E3 + 16*E1 - 30*E0 + 16*E2 - E4 )/(DR**2)
                                       Etris = (1.0d0/2.00d0)*(-E3 + 2*E1 -2*E2 + E4 )/(DR**3)
                                       
                                       IF ( Etris .NE. 0.0d0 ) THEN
                                               IF ( (Ebis/Etris)**2 .GE.  2*Eprim/Etris ) THEN
                                                       IF ( Etris .GT. 0.0d0 ) DRM = DRM - (Ebis/Etris) + sqrt( (Ebis/Etris)**2 - 2*Eprim/Etris )
                                                       IF ( Etris .LT. 0.0d0 ) DRM = DRM - (Ebis/Etris) - sqrt( (Ebis/Etris)**2 - 2*Eprim/Etris )
                                               ELSE
                                                       IF ( Ebis .NE. 0.0d0 ) DRM = DRM - Eprim/DABS(Ebis)
                                               ENDIF
                                       ELSE
                                               IF ( Ebis .NE. 0.0d0 ) DRM = DRM - Eprim/DABS(Ebis)
                                       ENDIF
                                       
                                       ! If the energy landscape is verry flatt, i.e, Eprim = Ebis = 0
                                       IF ( DRM .EQ. 0.0d0 .AND. leng .GT. FTol ) DRM = DRM + 0.10d0

                                       ! Safety cutoff:
                                       IF ( DABS(DRM) .GT. 10.0d0 ) DRM = DRM/DABS(DRM)
                                ENDIF

                          ENDDO
                          KK = KK + 1
                  ENDDO
          !=========================================================================================================================================

               IF ( I .EQ. 0 .AND. id .EQ. 0 ) THEN 
                        WRITE(*,'(A77)')'            ================================================================='
                        WRITE(*,'(A79)')'                             Relaxing the nuclear positions:                   '
                        WRITE(*,'(A77)')'            ================================================================='
                        WRITE(*,*)' '
                        WRITE(*,'(A114)')'  ITER           E [au]                        DE [au]                     <|F|> [au]                      DR [au]'
                ENDIF
                
                DE = (ETOT-EOLD)
                EOLD = ETOT
                leng = 0.0d0
                DO II=1,NATOMS
                        leng = sqrt(DOT_PRODUCT(force(II,:),force(II,:)))
                ENDDO
                
                leng = leng/NATOMS
                
                IF ( id .EQ. 0 ) WRITE(*,'(I4,E30.20,E30.20,E30.20,E30.20,I4)')I,ETOT,DE,leng,DRM,KK

                I = I + 1
                SCRATCH = .FALSE.

      ENDDO
        
      IF ( leng .LT. FTol .AND. id .EQ. 0 ) THEN
              WRITE(*,*)' '
              WRITE(*,'(A69)')'                 The nuclear positions have been relaxed with respect'
              WRITE(*,'(A68)')'                 to the interatomic forces. An force convergence of '
              WRITE(*,'(A27,E12.5,A24)')' F  < ',FTol,' [au] has been obtained.'
              WRITE(*,*)' ' 
              WRITE(*,'(A70)')'                The relaxed positions of the nuclea are the following:'
              WRITE(*,*)' '
              WRITE(*,'(A84)')'  ATOM           Rx [au]                       Ry [au]                       Rz [au]'
              DO J=1,NATOMS
                WRITE(*,'(I4,E30.20,E30.20,E30.20)')J,ATOMS(J)%R
              ENDDO
      ELSE
              IF ( id .EQ. 0 ) THEN
                        WRITE(*,*)' '
                        WRITE(*,'(A84)')'          Failure to relax the nuclear positions within prescribed force tolerance '
                        WRITE(*,*)' ' 
                        WRITE(*,'(A71)')'                     The last positions of the nuclea are the following:'
                        WRITE(*,*)' '
                        WRITE(*,'(A84)')'  ATOM           Rx [au]                       Ry [au]                       Rz [au]'
                        DO J=1,NATOMS
                                WRITE(*,'(I4,E30.20,E30.20,E30.20)')J,ATOMS(J)%R
                        ENDDO
                ENDIF
      ENDIF

      STOP

      END SUBROUTINE relax
