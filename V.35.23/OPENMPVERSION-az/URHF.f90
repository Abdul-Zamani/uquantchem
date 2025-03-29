SUBROUTINE URHF(S,H0,Intsv,NB,NRED,Ne,MULTIPLICITY,BSURHF,nucE,Tol,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,&
                & NSCF,FIXNSCF,POUT,SCRATCH,ZEROSCF,ETEMP,ENTROPY,NBAUX,VRI,WRI,RIAPPROX)
      ! This subroutine calculates the self consistent UN-Restricted
      ! Hartree-Fock solution
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: POUT,SCRATCH,ZEROSCF,RIAPPROX
      INTEGER, INTENT(IN) :: NB,Ne,FIXNSCF,NBAUX,MULTIPLICITY
      INTEGER*8, INTENT(IN) :: NRED
      INTEGER, INTENT(INOUT) :: DIISORD,DIISSTART
      INTEGER, INTENT(OUT) :: NSCF
      DOUBLE PRECISION, INTENT(IN) :: S(NB,NB),H0(NB,NB),Intsv(NRED),Tol,nucE,MIX(2),ETEMP,VRI(NBAUX,NBAUX),WRI(NB,NB,NBAUX)
      DOUBLE PRECISION, INTENT(OUT) :: EHFeigenup(NB),EHFeigendown(NB),ETOT,ENTROPY
      DOUBLE PRECISION, INTENT(INOUT) :: Cup(NB,NB),Cdown(NB,NB),Pup(NB,NB),Pdown(NB,NB)
      DOUBLE PRECISION :: PT(NB,NB),Jup(NB,NB),Jdown(NB,NB),Kup(NB,NB),Kdown(NB,NB),FTOT,FOLD,mu,TOLDNe
      DOUBLE PRECISION :: Fup(NB,NB),Fdown(NB,NB),G(NB,NB),C1(NB,NB),C2(NB,NB),C3(NB,NB),C4(NB,NB),DE,EOLD,DELTAP,LAMDAu,LAMDAd,MIXING
      DOUBLE PRECISION :: PTold(NB,NB),Pupold(NB,NB),Pdownold(NB,NB),Pups(50,NB,NB),Pdowns(50,NB,NB),Pupt(NB,NB),Pdownt(NB,NB)
      DOUBLE PRECISION :: ERRSU(50,NB,NB),ERRSD(50,NB,NB),ERRU(NB,NB),ERRD(NB,NB),SH(NB,NB),SL(NB,NB),LAM(NB),EIGENVECT(NB,NB),Sz,twoSP1,Nalpha,Nbeta
      DOUBLE PRECISION :: LSHIFTU(NB,NB),LSHIFTD(NB,NB),shift,CSHu(NB,NB),CSHd(NB,NB),FSHu(NB,NB),FSHd(NB,NB)
      INTEGER :: I,J,II,III,L,M,N,Neup,Nedown,INFO1,INFO2
      INTEGER :: MAXITER
      INTEGER, EXTERNAL :: ijkl
      LOGICAL :: STARTPRINTDIISIFO
!AZ 8/23
      DOUBLE PRECISION, ALLOCATABLE :: SABMO(:,:)
      DOUBLE PRECISION :: SSIJ,S2
      DOUBLE PRECISION :: CupOld(NB,NB), CdownOld(NB,NB) !for MOM 3/28/25
      LOGICAL :: BSURHF !momflag 
!AZ 8/23      

      MIXING  = MIX(1)
      shift = MIX(2)
      IF ( FIXNSCF .GT. 0 ) THEN
              MAXITER = FIXNSCF-1
      ELSE
              MAXITER = 1000
      ENDIF
      ! The tolerance used when calculating 
      ! the chemical potential 
      TOLDNe = 1.0E-8

      ! Calculating the lowdin S^(-1/2) matrix
      CALL diagh( S,NB,LAM,EIGENVECT)
      SL = 0.0d0
      DO I=1,NB
        SL(I,I) = 1.0d0/sqrt(LAM(I))
      ENDDO
      SH = MATMUL(EIGENVECT,MATMUL(SL,TRANSPOSE(EIGENVECT)))

      IF ( SCRATCH ) THEN
        CALL diaghHF( H0,S,NB,EHFeigenup,C1,INFO1)
        C2 = C1
      ELSE
        C1 = Cup
        C2 = Cdown
      ENDIF
      
      IF ( DIISORD .GT. 25 ) THEN
              WRITE(*,*)'The order of DIIS mixing set =',DIISORD,' Is too Big.'
              WRITE(*,*)'Maximum allowed value is 25. This value has been used instead.'
              DIISORD = 25
      ENDIF

      !AZ 
      !fix double precision print
      !2sp1 = 2s+1 
      twoSP1 = Multiplicity
      Sz = ((twoSP1) - 1) / 2
      if(Sz.eq.0) then
        Neup = (Ne / 2)
        Nedown = (Ne / 2)
      elseif(Sz.ne.0) then
        Nbeta  = (Ne-(Sz/0.5)) / 2
        Nalpha = NBeta + (Sz/0.5)
        Neup = Nalpha
        Nedown=Nbeta
      endif

      write(*,*)
      print*,'Sz',Sz
      print*,'Multiplicity',MULTIPLICITY
      print*,'NAlpha',Neup
      print*,'NBeta',Nedown
      print*,'NBas',NB
      write(*,*)

!      Neup   = ( Ne - MOD(Ne,2) )/2
!      Nedown = ( Ne + MOD(Ne,2) )/2

      Pups = 0.0d0
      Pdowns = 0.0d0

      DE = 2.0d0*Tol
      I = 0
      II = 0
      PTold = 0.0d0
      STARTPRINTDIISIFO = .FALSE.
        
      !=======================================================
      ! This is used for the zero scf cycle (FAST-MD)
      !=======================================================
      IF ( ZEROSCF ) THEN
                PT = Pup + Pdown
                CALL getJv(Pup,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Jup)
                CALL getJv(Pdown,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Jdown)
                CALL getKv(Pup,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Kup)
                CALL getKv(Pdown,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Kdown)
                Fup   = H0 + Jdown - Kup   + Jup
                Fdown = H0 + Jdown - Kdown + Jup
                CALL diaghHF( Fup,S,NB,EHFeigenup,C1,INFO1)
                CALL diaghHF( Fdown,S,NB,EHFeigendown,C2,INFO2)
                !ETOT = 0.50d0*(SUM(H0*PT) + SUM(Fup*Pup) + SUM(Fdown*Pdown) ) + nucE
                IF ( ETEMP .LT. 0.0d0 ) THEN
                        Cup(:,:) = 0.0d0
                        DO M=1,Neup
                                Cup(:,M) = C1(:,M)
                        ENDDO
                        Cdown(:,:) = 0.0d0
                        DO M=1,Nedown
                                Cdown(:,M) = C2(:,M)
                        ENDDO
                        CALL makedens(Cup,NB,Pup)
                        CALL makedens(Cdown,NB,Pdown)
                ELSE
                        CALL makedensT(TOLDNe,Cup,Cdown,EHFeigenup,EHFeigendown,ETEMP,NB,Ne,Pup,Pdown,mu,ENTROPY)
                ENDIF
                NSCF = 0
                RETURN
      ENDIF
      !=====================================================================
      
      DO WHILE ( (DABS(DE) .GT. Tol .AND. I .LE. MAXITER) .OR. ( DELTAP .GT. sqrt(Tol) .AND. I .LE. MAXITER) )
                IF ( ETEMP .LE. 0.0d0 ) THEN
                        Cup(:,:) = 0.0d0
                        DO M=1,Neup
                                Cup(:,M) = C1(:,M)
                        ENDDO

                        Cdown(:,:) = 0.0d0
                        DO M=1,Nedown
                                Cdown(:,M) = C2(:,M)
                        ENDDO
                ELSE
                        Cup = C1
                        Cdown = C2
                ENDIF
               
                !==============================================================
                ! In the case of XL-BOMD, it is crusial that the scf starts 
                ! with the density matrices Pup and Pdown provided by the input,
                ! since it is these matrices that have been propagated by the 
                ! XL-BOMD scheme and not the coefficient matrices Cup and Cdown
                ! Therefor we have the if-construct here.
                !==============================================================

                IF ( I .GT. 0 ) THEN
                        IF ( ETEMP .LE. 0.0d0 ) THEN
                                if(BSURHF) then
                                  if(I.eq.1) then
                                  !AZ 3/29 hard code imom
                                  cupold = cup
                                  !cdownold = cdown 
                                  cupOld(:,1) =  cupOld(:,1)*0.5d0
                                  !cdownOld(:,1) =  cdownOld(:,1)*0.5d0 
                                  cup=cup*.9
                                  !cdown=cdown*.9 
                                  endif 
                                  if(I.ge.1) then 
                                  call mom(Cup,CupOld,S,NB,Neup,Pup,1)
                                  !call mom(Cdown,CdownOld,S,NB,Nedown,Pdown,2)
                                  CALL makedens(Cdown,NB,Pdown)
                                  endif 
                                else 
                                CALL makedens(Cup,NB,Pup)
                                CALL makedens(Cdown,NB,Pdown)
                                endif 
 
                        ELSE
                                CALL makedensT(TOLDNe,Cup,Cdown,EHFeigenup,EHFeigendown,ETEMP,NB,Ne,Pup,Pdown,mu,ENTROPY)
                        ENDIF
                ENDIF

                PT = Pup + Pdown
            
                ! Calculating the change of the total density matrix:
                DELTAP = sqrt(DOT_PRODUCT(reshape(PT-PTold,(/NB**2/)),reshape(PT-PTold,(/NB**2/))))
                
                ! Linear mixing:
                IF ( I .GT. 0 ) THEN
                        Pup = Pup*(1.0d0-MIXING) + Pupold*MIXING
                        Pdown = Pdown*(1.0d0-MIXING) + Pdownold*MIXING
                ELSE
                        Pup = Pup
                        Pdown = Pdown
               ENDIF

                Pupold = Pup
                Pdownold = Pdown

                PTold = Pupold + Pdownold

                ! Saving the density matrices from the previous
                ! iterations, to be used with the DIIS-method.
                IF ( II .LT. 2*DIISORD ) THEN
                        II = II + 1
                        Pups(II,:,:) = Pup
                        Pdowns(II,:,:) = Pdown
                        ERRSU(II,:,:) = ERRU
                        ERRSD(II,:,:) = ERRD
                ELSE
                        ! refreshing the density matrices used by the DIIS
                        IF ( DIISORD  .NE. 0 ) THEN
                                DO III=1,II-1
                                        Pups(III,:,:) = Pups(III+1,:,:)
                                        Pdowns(III,:,:) = Pdowns(III+1,:,:)
                                        ERRSU(III,:,:) = ERRSU(III+1,:,:)
                                        ERRSD(III,:,:) = ERRSD(III+1,:,:)
                                ENDDO
                                Pups(II,:,:) = Pup
                                Pdowns(II,:,:)  = Pdown
                                ERRSU(II,:,:) = ERRU
                                ERRSD(II,:,:) = ERRD
                        ENDIF
                ENDIF

                ! Here we use the DIIS method ( CHEM. Phys. Lett. 73, 393 (1980) )
                ! in order to estimate the self-consistent density matrices:
                IF ( II .GE. 4 ) THEN

                        CALL DIIS(NB,II,ERRSU,Pups,Pupt,LAMDAu,INFO1)
                        CALL DIIS(NB,II,ERRSD,Pdowns,Pdownt,LAMDAd,INFO2)
                        
                        IF ( (INFO1 .EQ. 0 .AND. INFO2 .EQ. 0 .AND. LAMDAu+LAMDAd .LT. 1.0d0) .OR. (INFO1 .EQ. 0 .AND. INFO2 .EQ. 0 .AND. I .GT. DIISSTART )  ) THEN
                                IF ( I .GE. DIISSTART ) THEN
                                        Pup = Pupt
                                        Pdown = Pdownt
                                        STARTPRINTDIISIFO = .TRUE.
                                        MIXING = 0.0d0
                                ENDIF
                        ENDIF
                ENDIF

                PT = Pup + Pdown

                CALL getJv(Pup,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Jup)

                CALL getJv(Pdown,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Jdown)
                
                CALL getKv(Pup,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Kup)
                
                CALL getKv(Pdown,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Kdown)
                
                Fup   = H0 + Jdown - Kup   + Jup
                Fdown = H0 + Jdown - Kdown + Jup
               
                ! Here we calculate the error vector/matrix according to 
                ! P. Pulay, J. Comp. Chem. 3, 556 (1982) (Eqn 4)
                ! to be used in the DIIS algorithm:
                IF ( DIISORD .NE. 0 ) THEN
                        ERRU = MATMUL(Fup,MATMUL(Pup,S)) - MATMUL(S,MATMUL(Pup,Fup))
                        ERRD = MATMUL(Fdown,MATMUL(Pdown,S)) - MATMUL(S,MATMUL(Pdown,Fdown))
                        ERRU = MATMUL(TRANSPOSE(SH),MATMUL(ERRU,SH))
                        ERRD = MATMUL(TRANSPOSE(SH),MATMUL(ERRD,SH))
                ENDIF

                CALL diaghHF( Fup,S,NB,EHFeigenup,C1,INFO1)
                CALL diaghHF( Fdown,S,NB,EHFeigendown,C2,INFO2)

                !----------------------------------------------------------------------------------
                ! Here we do a level shift of the LUMO -states in order to try and accelerate the
                ! convergence of the scf calculation. See Theoret. Chim. Acta. 46, 325-329 (1977).
                ! WARNING! when using this method with FERMI-DIRAC occupation (ETEMP > 0 ) The
                ! calculated entropy,internal energy and free energy  loose there physical meaning.
                !----------------------------------------------------------------------------------
                IF ( shift .GT. 0.0d0 ) THEN
                    LSHIFTU(:,:) = 0.0d0
                    LSHIFTD(:,:) = 0.0d0

                    DO M=1,Neup + 1,NB
                       LSHIFTU(M,M) = shift
                    ENDDO

                    DO M=1,Nedown + 1,NB
                       LSHIFTD(M,M) = shift
                    ENDDO

                    IF ( I .EQ. 0 ) THEN
                       CSHu(:,:) = 0.0d0
                       CSHd(:,:) = 0.0d0
                    ENDIF

                    FSHu = MATMUL(SH,MATMUL(Fup,SH)) + MATMUL(CSHu,MATMUL(LSHIFTU,TRANSPOSE(CSHu)))
                    FSHd = MATMUL(SH,MATMUL(Fdown,SH)) + MATMUL(CSHd,MATMUL(LSHIFTD,TRANSPOSE(CSHd)))

                    CALL diagh(FSHu,NB,EHFeigenup,CSHu,INFO1)
                    CALL diagh(FSHd,NB,EHFeigendown,CSHd,INFO2)

                    C1 = MATMUL(SH,CSHu)
                    C2 = MATMUL(SH,CSHd)
                ENDIF

                IF ( FIXNSCF .LT. 0 ) THEN
                        ETOT = 0.50d0*(SUM(H0*PT) + SUM(Fup*Pup) + SUM(Fdown*Pdown) ) + nucE
                        IF ( ETEMP .GT. 0.0d0  ) FTOT = ETOT - ETEMP*ENTROPY
                ENDIF

                IF ( I .EQ. 0 .AND. POUT ) THEN
                        print*,'   =========================================================='
                        print*,'         Entering the scf unrestricted Hartree-Fock loop      '
                        print*,'   =========================================================='
                        print*,' '
                        IF ( ETEMP .LT. 0.0d0  ) WRITE(*,'(A4,A22,A27,A27,A33)')'N','E [au]','DE [au]',' DP','DIIS'
                        IF ( ETEMP .GT. 0.0d0  ) WRITE(*,'(A4,A20,A29,A31,A27,A31)')'N','E [au]','F [au]','DF [au]',' DP','DIIS'
                ENDIF
        
                IF ( I .GT. 0 .AND. FIXNSCF .LE. 0 .AND. ETEMP .LT. 0.0d0 ) DE = ETOT-EOLD
                IF ( I .GT. 0 .AND. FIXNSCF .LE. 0 .AND. ETEMP .GT. 0.0d0 ) DE = FTOT-FOLD
                
                IF ( STARTPRINTDIISIFO ) THEN
                        IF ( POUT .AND. ETEMP .GT. 0.0d0 ) WRITE(*,'(I4,E30.20,E30.20,E30.20,E30.20,E30.20)'),I,ETOT,FTOT,DE,DELTAP,LAMDAd+LAMDAu
                        IF ( POUT .AND. ETEMP .LT. 0.0d0 ) WRITE(*,'(I4,E30.20,E30.20,E30.20,E30.20)'),I,ETOT,DE,DELTAP,LAMDAd+LAMDAu
                ELSE
                        IF ( POUT .AND. ETEMP .GT. 0.0d0 ) WRITE(*,'(I4,E30.20,E30.20,E30.20,E30.20)'),I,ETOT,FTOT,DE,DELTAP
                        IF ( POUT .AND. ETEMP .LT. 0.0d0 ) WRITE(*,'(I4,E30.20,E30.20,E30.20)'),I,ETOT,DE,DELTAP
                ENDIF
        
                EOLD = ETOT
                FOLD = FTOT
                I = I+1

     ENDDO
     Cup = C1
     Cdown = C2
     IF ( ETEMP .LT. 0.0d0 ) THEN
        C3 = 0.0d0
        DO M=1,Neup
           C3(:,M) = C1(:,M)
        ENDDO
        C4 = 0.0d0
        DO M=1,Nedown
           C4(:,M) = C2(:,M)
        ENDDO
        CALL makedens(C3,NB,Pup)
        CALL makedens(C4,NB,Pdown)
      ELSE
        CALL makedensT(TOLDNe,Cup,Cdown,EHFeigenup,EHFeigendown,ETEMP,NB,Ne,Pup,Pdown,mu,ENTROPY)
      ENDIF
     NSCF = I

     IF ( I .LE. MAXITER .OR. FIXNSCF .GT. 0 ) THEN
             IF ( POUT ) THEN
                print*,' '
                WRITE(*,'(A54)')'                 Convergence reached within tolerance:'
                WRITE(*,'(A22,E8.1,A23)')'Tol=',Tol,' au .Aborting scf loop.'
                print*,' '
                IF ( ETOT .GT. -1.0E03) WRITE(*,'(A33,E27.20,A3)'),'      Hartree-Fock energy:   E = ',ETOT,' au'
                IF ( ETOT .LT. -1.0E03) WRITE(*,'(A33,E30.20,A3)'),'      Hartree-Fock energy:   E = ',ETOT,' au'
                IF ( FTOT .GT. -1.0E03 .AND. ETEMP .GT. 0.0d0 ) WRITE(*,'(A7,A6,A20,E27.20,A3)'),'       ',' URHF ',' free-energy:   F = ',FTOT,' au'
                IF ( FTOT .LT. -1.0E03 .AND. ETEMP .GT. 0.0d0 ) WRITE(*,'(A7,A6,A20,E30.20,A3)'),'       ',' URHF ',' free-energy:   F = ',FTOT,' au'
                print*,' '
!AZ 8/23
                        ALLOCATE(SABMO(NB,NB))
                        !SABMO  = Ca*.S.Cb 
                        SABMO = matmul(matmul(transpose(Cup),S),Cdown)
                        SSIJ=0.0d0
                        do i=1, neup !nalpha
                          do j=1, nedown !nbeta
                             SSIJ = SSIJ + SABMO(i,j)**(2.0d0)
                          enddo
                        enddo
                            !S2  = sz(sz+1) + nB - |Sij|^2 
                            S2  = Sz*(Sz+1.0d0) + (Nedown*1.0d0) - SSIJ
                            print*,'Nalpha',Neup
                            print*,'Nbeta',Nedown
                            print*,'SSIJ',SSIJ
                            print*,'Sz(Sz+1)=', Sz*(Sz+1.0d0)
                            print*,'S2=',S2


                        DEALLOCATE(SABMO)
!AZ 8/23
             ENDIF
     ELSE
             print*,'---------------------------------------------------------'
             print*,'CALCULATION FAILED TO CONVERGE WITHIN ,',Tol,'au'
             print*,'           UQUANTCHEM ABORTS!                    '
             print*,'---------------------------------------------------------'
             STOP
     ENDIF
     IF ( POUT ) THEN
        OPEN(21,FILE='UHFEIGENVALUES.dat',ACTION='WRITE')
        DO I=1,NB
                WRITE(21,*)I,EHFeigenup(I),EHFeigendown(I)
        ENDDO
        CLOSE(21)
     ENDIF
!AZ 3/28 write mo coeff
     !WRITE(*,*) Cup(:,1)
     IF ( POUT ) THEN
        OPEN(24,FILE='alphaMO.dat',ACTION='WRITE')
        DO I=1,NB
          DO J=1,NB
                WRITE(24,*) Cup(I,J)!Cup(J,I) 
          ENDDO 
        ENDDO
        CLOSE(24)
     ENDIF
!AZ 3/28 write mo coeff 
     IF ( ETEMP .GT. 0.0d0 ) ETOT = FTOT
END SUBROUTINE URHF

