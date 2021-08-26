SUBROUTINE DKH(S,T,V,PVP,NB,DKHORDER,HDKH)
      ! This subroutine calculates the scalar relatevistic one-electron hamiltonian
      ! upp to DKHORDER <= 4 and returns the corresponding one-electron scalar relativistic
      ! hamiltonian in HDKH. We here follow strictly the recepie found in 
      ! A. Wolf, M. Reiher and B. A. Hess in "The generalized Douglas-Kroll transformation",
      ! J. Chem. Phys. 117, 9215 (2002).
      ! Here the variables are the following:
      !---------------------------------------------------------------------------------------------------
      ! id = mpi processor id. Used to avoid multiple output.
      ! NB = Basis set size.
      ! DKHORDER = Order of DKH treatment to be used.
      ! T = kinetic energy matrix
      ! V = Potential energy matrix
      ! PVP = <grad(PHI_I) | V | grad(PHI_J) > see Eqn (A5) in J. Chem. Phys. 117, 9215 (2002).
      ! E0, E1, .., E6 are the DKH iscalar relatevistic version of energy terms found 
      ! in Eqn(9,10,44,46,57,58) in J. Chem. Phys. 117, 9215 (2002).
      ! For example the scalar relatevistic version of E2 (Eq 61, J. Chem. Phys. 117, 9215 (2002) )
      ! can be found in Eqn (112) of Chap 11, Relatevistic Electronic Structure Theory, Part 1 
      ! fundamentals Theoretical and Computational Chemistry, Vol 11, by A. Wolf, M. Reiher and B. A. Hess
      !----------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NB
      INTEGER, INTENT(INOUT) :: DKHORDER
      DOUBLE PRECISION, INTENT(IN) :: S(NB,NB),T(NB,NB),V(NB,NB),PVP(NB,NB)
      DOUBLE PRECISION, INTENT(OUT) :: HDKH(NB,NB) ! The DKH skalar relativistic correspondens of H0 = T + V
      DOUBLE PRECISION :: p2(NB),OMEGA(NB,NB),AP(NB,NB),RP(NB,NB),Vp2(NB,NB),PVPp2(NB,NB),Vt(NB,NB),PVPt(NB,NB)
      DOUBLE PRECISION :: IPP(NB,NB), PP(NB,NB), TESTOMEGA(NB,NB),SHI(NB,NB),SH(NB,NB),LAM(NB),DIAG(NB,NB),DIAGI(NB,NB)
      DOUBLE PRECISION :: E0(NB,NB), E1(NB,NB), E2(NB,NB), E3(NB,NB), E4(NB,NB), E5(NB,NB), EP(NB,NB), Et(NB,NB)
      DOUBLE PRECISION :: BE3(NB,NB), CE3(NB,NB), DE3(NB,NB), EE3(NB,NB), RE1R(NB,NB) ! Auxilliary variables for calculation of E3
      DOUBLE PRECISION :: APVPA(NB,NB), APVtPA(NB,NB), AVA(NB,NB), AVtA(NB,NB), TEMP(NB,NB),TERM(NB,NB), DE(NB,NB)
      DOUBLE PRECISION, PARAMETER :: C = 137.0359895 ! Speed of light.
      INTEGER :: I, J, INFO
      INTEGER :: ierr
      LOGICAL :: TEST, ALT

      TEST = .FALSE.
      ALT = .TRUE.

      DKHORDER = ABS(DKHORDER)

      IF ( DKHORDER > 4 ) THEN
        WRITE(*,*)'WARNING! Maximum order (DKHORDER) of scalar relativistic calculations is currently 4'
        WRITE(*,*)'Setting DKHORDER to 4'
        DKHORDER = 4
      ENDIF

      IF ( DKHORDER > 0 ) THEN
             print*,'   ----------------------------------------------------------'
             WRITE(*,'(A59,I2)') '    PERFORMING DKH SCALAR RELATIVISTIC CALCULATION OF ORDER', DKHORDER
             print*,'      See J. Chem. Phys. 117, 9215 (2002) for details.'
             print*,'   ----------------------------------------------------------'
      ENDIF

      HDKH(:,:) = 0.0d0
      EP(:,:) = 0.0d0
      RP(:,:) = 0.0d0
      AP(:,:) = 0.0d0

      !----------------------------------------
      ! Preparing the transformation matrix for
      ! transforming back and forth to the p**2
      ! representation. See Eqn (A2) in 
      ! J. Chem. Phys. 117, 9215 (2002)
      !-----------------------------------------
      
      CALL diagh(S,NB,LAM,OMEGA,INFO)

      !---------------------------------------------------------------------
      ! Calculating the Lovin-matrices S^(1/2) and S(-1/2)
      ! to transform the matrices V,T,PVP into an orthogonal representation
      !---------------------------------------------------------------------

      DIAG(:,:) = 0.0d0
      DIAGI(:,:) = 0.0d0
      DO I = 1,NB
                DIAG(I,I) = sqrt(LAM(I))
                DIAGI(I,I) = 1.0d0/sqrt(LAM(I))
      ENDDO

      SH = MATMUL(OMEGA,MATMUL(DIAG,TRANSPOSE(OMEGA)))
      SHI = MATMUL(OMEGA,MATMUL(DIAGI,TRANSPOSE(OMEGA)))

      CALL diagh(MATMUL(SHI,MATMUL(T,SHI)),NB,p2,OMEGA,INFO)

      !-----------------------------------------------
      ! Transforming V and pVp to p**2 representation:
      !-----------------------------------------------

      Vp2 = MATMUL(TRANSPOSE(OMEGA),MATMUL(MATMUL(SHI,MATMUL(V,SHI)),OMEGA))
      PVPp2 = MATMUL(TRANSPOSE(OMEGA),MATMUL(MATMUL(SHI,MATMUL(PVP,SHI)),OMEGA))
      
     
      !----------------------------------------------------------
      ! Testing Orthonormality of the vectors contained in OMEGA:
      !----------------------------------------------------------
      IF ( TEST ) THEN
                print*,' '
                print*,'-------------------------------------------------------'
                print*,'TESTING THE ORTHONORMALITY OF THE p**2 REPRESENTATION'
                print*,'-------------------------------------------------------'
                print*,' '
                TESTOMEGA = MATMUL(TRANSPOSE(OMEGA),OMEGA)
                DO I=1,NB
                        DO J=1,NB
                                IF (DABS(TESTOMEGA(I,J)) .GT. 1.0E-10 ) print*,I,J,':','<P(i)|P(j)> = ',TESTOMEGA(I,J)
                        ENDDO
                ENDDO
                print*,'--------------------------------------------------------'
                print*,' '
      ENDIF

      !----------------------------------------------------------
      ! Calculating the EP, Rp, E0 and AP -matrices
      ! See Eqn (6,7,9) in J. Chem. Phys. 117, 9215 (2002)
      !----------------------------------------------------------

      E0(:,:) = 0.0d0
      PP(:,:) = 0.0d0
      IPP(:,:) = 0.0d0

      DO I=1,NB
        EP(I,I) = sqrt( 2.0d0*p2(I)*C**2 + C**4 )
        AP(I,I) = sqrt( (EP(I,I) + C**2 )/(2.0d0*EP(I,I)) )
        RP(I,I) = C/(EP(I,I) + C**2 )
        E0(I,I) = EP(I,I) - C**2
        PP(I,I) = (2.0d0*p2(I)*C**2)/( (EP(I,I) + C**2 )**2 )
        IPP(I,I) = 1.0d0/PP(I,I)
      ENDDO
      
      !--------------------------------------------------------------------------
      ! Calculating V-tilde defined by Eqn(63) in J. Chem. Phys. 117, 9215 (2002)
      ! And the corresponding tilde matrix for PVP.
      !--------------------------------------------------------------------------
      DO I=1,NB
         DO J=1,NB
                Vt(I,J) = Vp2(I,J)/( EP(I,I) + EP(J,J) )
                PVPt(I,J) = PVPp2(I,J)/( EP(I,I) + EP(J,J) )
                !E1(I,J) = AP(I,I)*VP2(I,J)*AP(J,J) + AP(I,I)*RP(I,I)*PVPp2(I,J)*RP(J,J)*AP(J,J)
        ENDDO
      ENDDO
      
      !-------------------------------------------------------------------------
      ! Here we calculate the E1 contribution to the DKH -Hamiltonian
      ! according to equation (10) in  J. Chem. Phys. 117, 9215 (2002), 
      ! or from an implementational more close perspective according 
      ! to eqn (9) page 12 in the blue folder "Jonisering + DKH + Atomic-Sphere"
      !--------------------------------------------------------------------------

      IF ( DKHORDER .GT. 0 ) THEN

        IF (DKHORDER .GE. 1 ) THEN
                E1 = MATMUL(MATMUL(AP,Vp2 + MATMUL(RP,MATMUL(PVPp2,RP))),AP)
                HDKH = E0 + E1
        ENDIF

        IF (DKHORDER .GE. 2 ) THEN
                E2 = 0.50d0*( -MATMUL(AP,MATMUL(RP,MATMUL(PVPt,MATMUL(RP,MATMUL(AP,MATMUL(AP,MATMUL(Vp2,AP))))))) & 
                & + MATMUL(AP,MATMUL(RP,MATMUL(PVPt,MATMUL(RP,MATMUL(AP,MATMUL(IPP,MATMUL(AP,MATMUL(RP,MATMUL(PVPp2,MATMUL(RP,AP)))))))))) &
                & + MATMUL(AP,MATMUL(Vt,MATMUL(AP,MATMUL(PP,MATMUL(AP,MATMUL(Vp2,AP)))))) &
                & - MATMUL(AP,MATMUL(Vt,MATMUL(AP,MATMUL(AP,MATMUL(RP,MATMUL(PVPp2,MATMUL(RP,AP))))))) &
                & - MATMUL(AP,MATMUL(RP,MATMUL(PVPp2,MATMUL(RP,MATMUL(AP,MATMUL(AP,MATMUL(Vt,AP))))))) &
                & + MATMUL(AP,MATMUL(RP,MATMUL(PVPp2,MATMUL(RP,MATMUL(AP,MATMUL(IPP,MATMUL(AP,MATMUL(RP,MATMUL(PVPt,MATMUL(RP,AP)))))))))) &
                & + MATMUL(AP,MATMUL(Vp2,MATMUL(AP,MATMUL(PP,MATMUL(AP,MATMUL(Vt,AP)))))) &
                & - MATMUL(AP,MATMUL(Vp2,MATMUL(AP,MATMUL(AP,MATMUL(RP,MATMUL(PVPt,MATMUL(RP,AP))))))) )

                HDKH = HDKH + E2
        ENDIF
 
        IF (DKHORDER .GE. 3 ) THEN
                !-----------------------------------------------------------------------------------------------------------
                ! Initialization of auxilliary variables. Here we use the expressions for the third order DKH correction as
                ! has been deduced in my notes page 15-18, in the blue folder: "Jonisation + DKH + Atomic Sphere".
                !-----------------------------------------------------------------------------------------------------------
                Et(:,:) = 0.0d0
                BE3(:,:) = 0.0d0
                CE3(:,:) = 0.0d0 
                DE3(:,:) = 0.0d0
                EE3(:,:) = 0.0d0 
                RE1R(:,:) = 0.0d0

                
                PVPp2 = MATMUL(RP,MATMUL(PVPp2,RP))
                PVPt = MATMUL(RP,MATMUL(PVPt,RP))                

                IF ( .not. ALT ) THEN 
                        
                        Et = MATMUL(AP,MATMUL(PVPt,MATMUL(AP,MATMUL(AP,MATMUL(Vt,AP))))) & 
                        &  + MATMUL(AP,MATMUL(Vt,MATMUL(AP,MATMUL(AP,MATMUL(PVPt,AP))))) &
                        &  - MATMUL(AP,MATMUL(PVPt,MATMUL(AP,MATMUL(IPP,MATMUL(AP,MATMUL(PVPt,AP)))))) &
                        &  - MATMUL(AP,MATMUL(Vt,MATMUL(AP,MATMUL(PP,MATMUL(AP,MATMUL(Vt,AP)))))) 

                        !-----------------------------------------------------------------------------------------------------
                        ! Here we use the expression for W1*E1*W1 found in J. Chem. Phys. 113, 7786 (2000) Eqn (24) - Eqn(26)
                        !-----------------------------------------------------------------------------------------------------
                        
                        !------------------------------------
                        ! W1*E1*W1 = ( BE3 + CE3 - DE3 - EE3)
                        !------------------------------------
                        
                        RE1R = MATMUL(AP,MATMUL(PVPp2,AP)) + MATMUL(AP,MATMUL(PP,MATMUL(Vp2,MATMUL(PP,AP))))

                        BE3 = MATMUL(AP,MATMUL(PVPt,MATMUL(AP,MATMUL(IPP,MATMUL(RE1R,MATMUL(AP,MATMUL(Vt,AP)))))))

                        DE3 = MATMUL(AP,MATMUL(PVPt,MATMUL(AP,MATMUL(IPP,MATMUL(RE1R,MATMUL(IPP,MATMUL(AP,MATMUL(PVPt,AP))))))))

                        EE3 = MATMUL(AP,MATMUL(Vt,MATMUL(AP,MATMUL(RE1R,MATMUL(AP,MATMUL(Vt,AP))))))

                        CE3 = MATMUL(AP,MATMUL(Vt,MATMUL(AP,MATMUL(RE1R,MATMUL(IPP,MATMUL(AP,MATMUL(PVPt,AP)))))))
                
                        !BE3 = MATMUL(AP,MATMUL(PVPt,MATMUL(AP,MATMUL(IPP,MATMUL(AP,MATMUL(PVPp2,MATMUL(AP,MATMUL(AP,MATMUL(Vt,AP))))))))) &
                        !&   + MATMUL(AP,MATMUL(PVPt,MATMUL(AP,MATMUL(AP,MATMUL(Vp2,MATMUL(AP,MATMUL(PP,MATMUL(AP,MATMUL(Vt,AP))))))))) 
                
                        !CE3 = MATMUL(AP,MATMUL(Vt,MATMUL(AP,MATMUL(AP,MATMUL(PVPp2,MATMUL(AP,MATMUL(IPP,MATMUL(AP,MATMUL(PVPt,AP))))))))) &
                        !&   + MATMUL(AP,MATMUL(Vt,MATMUL(AP,MATMUL(PP,MATMUL(AP,MATMUL(Vp2,MATMUL(AP,MATMUL(AP,MATMUL(PVPt,AP)))))))))

                        !DE3 = MATMUL(AP,MATMUL(PVPt,MATMUL(AP,MATMUL(IPP,MATMUL(AP,MATMUL(PVPp2,MATMUL(AP,MATMUL(IPP,MATMUL(AP,MATMUL(PVPt,AP)))))))))) &
                        !&   + MATMUL(AP,MATMUL(PVPt,MATMUL(AP,MATMUL(AP,MATMUL(Vp2,MATMUL(AP,MATMUL(AP,MATMUL(PVPt,AP))))))))

                        !EE3 = MATMUL(AP,MATMUL(Vt,MATMUL(AP,MATMUL(AP,MATMUL(PVPp2,MATMUL(AP,MATMUL(AP,MATMUL(Vt,AP))))))))  &
                        !&   + MATMUL(AP,MATMUL(Vt,MATMUL(AP,MATMUL(PP,MATMUL(AP,MATMUL(Vp2,MATMUL(AP,MATMUL(PP,MATMUL(AP,MATMUL(Vt,AP))))))))))

                        E3 = 0.50d0*( MATMUL(Et,E1) + MATMUL(E1,Et) - 2.0d0*( BE3 + CE3 - DE3 - EE3) )
                ELSE
                        !-----------------------------------------------------------------------------
                        ! Here we use the explicit expression for the third order Douglas-Kroll energy
                        ! Given by equation (A20) in J. Chem. Phys. 130, 124103 (2009).
                        !-----------------------------------------------------------------------------
                        APVPA  = MATMUL(AP,MATMUL(PVPp2,AP)) 
                        APVtPA = MATMUL(AP,MATMUL(PVPt,AP)) 
                        AVA    = MATMUL(AP,MATMUL(Vp2,AP))
                        AVtA   = MATMUL(AP,MATMUL(Vt,AP))

                        E3 = 0.50d0*( MATMUL(APVtPA,MATMUL(AVtA,AVA)) - MATMUL(APVtPA,MATMUL(IPP,MATMUL(APVtPA,AVA))) - MATMUL(AVtA,MATMUL(PP,MATMUL(AVtA,AVA))) &
                        &  + MATMUL(AVtA,MATMUL(APVtPA,AVA)) + MATMUL(APVtPA,MATMUL(AVtA,APVPA)) - MATMUL(APVtPA,MATMUL(IPP,MATMUL(APVtPA,APVPA))) &
                        &  - MATMUL(AVtA,MATMUL(PP,MATMUL(AVtA,APVPA))) + MATMUL(AVtA,MATMUL(APVtPA,APVPA)) & 
                        &  - 2*MATMUL(APVtPA,MATMUL(IPP,MATMUL(APVPA,AVtA))) + 2*MATMUL(APVtPA,MATMUL(IPP,MATMUL(APVPA,MATMUL(IPP,APVtPA)))) & 
                        &  - 2*MATMUL(APVtPA,MATMUL(AVA,MATMUL(PP,AVtA))) + 2*MATMUL(APVtPA,MATMUL(AVA,APVtPA)) &
                        &  + 2*MATMUL(AVtA,MATMUL(APVPA,AVtA)) - 2*MATMUL(AVtA,MATMUL(APVPA,MATMUL(IPP,APVtPA))) &
                        &  + 2*MATMUL(AVtA,MATMUL(PP,MATMUL(AVA,MATMUL(PP,AVtA)))) &
                        &  - 2*MATMUL(AVtA,MATMUL(PP,MATMUL(AVA,APVtPA))) + MATMUL(AVA,MATMUL(APVtPA,AVtA)) + MATMUL(APVPA,MATMUL(APVtPA,AVtA)) &
                        &  - MATMUL(AVA,MATMUL(APVtPA,MATMUL(IPP,APVtPA))) - MATMUL(APVPA,MATMUL(AVtA,MATMUL(PP,AVtA))) &
                        &  - MATMUL(AVA,MATMUL(AVtA,MATMUL(PP,AVtA))) - MATMUL(APVPA,MATMUL(APVtPA,MATMUL(IPP,APVtPA))) &
                        &  + MATMUL(AVA,MATMUL(AVtA,APVtPA)) + MATMUL(APVPA,MATMUL(AVtA,APVtPA)) )
                ENDIF
                HDKH = HDKH + E3
        ENDIF

        IF (DKHORDER .GE. 4 ) THEN
                Et = MATMUL(AP,MATMUL(PVPt,MATMUL(AP,MATMUL(AP,MATMUL(Vt,AP))))) & 
                &  + MATMUL(AP,MATMUL(Vt,MATMUL(AP,MATMUL(AP,MATMUL(PVPt,AP))))) &
                &  - MATMUL(AP,MATMUL(PVPt,MATMUL(AP,MATMUL(IPP,MATMUL(AP,MATMUL(PVPt,AP)))))) &
                &  - MATMUL(AP,MATMUL(Vt,MATMUL(AP,MATMUL(PP,MATMUL(AP,MATMUL(Vt,AP)))))) 
         
                APVPA  = MATMUL(AP,MATMUL(PVPp2,AP)) 
                APVtPA = MATMUL(AP,MATMUL(PVPt,AP)) 
                AVA    = MATMUL(AP,MATMUL(Vp2,AP))
                AVtA   = MATMUL(AP,MATMUL(Vt,AP))

                !-------------------------------------------------
                ! The term (a) on pages 23-24 in my notes for DKH:
                !-------------------------------------------------
                E4 = 0.250d0*( MATMUL(Et,E2) + MATMUL(E2,Et) )
                
                TEMP = MATMUL(APVPA,AVtA) + MATMUL(AVA,APVtPA) - MATMUL(APVPA,MATMUL(IPP,APVtPA)) - MATMUL(AVA,MATMUL(PP,AVtA)) 
                E4 = E4 - 0.250d0*MATMUL(Et,TEMP)
                TEMP = -MATMUL(APVtPA,AVA) - MATMUL(AVtA,APVPA) + MATMUL(APVtPA,MATMUL(IPP,APVPA)) + MATMUL(AVtA,MATMUL(PP,AVA))
                E4 = E4 + 0.250d0*MATMUL(TEMP,Et)

                !---------------------------------------------
                ! The term (c) on pages 24-25 in my DKH notes:
                !---------------------------------------------
                DE(:,:) = 0.0d0
                TERM(:,:) = 0.0d0

                RE1R = MATMUL(AP,MATMUL(PVPp2,AP)) + MATMUL(AP,MATMUL(PP,MATMUL(Vp2,MATMUL(PP,AP))))
                
                !----------------------
                ! Equation (9) page 24:
                !----------------------
                DO I=1,NB
                        DO J=1,NB
                                DE(J,J) = 1.0d0/( EP(J,J) + EP(I,I) )
                        ENDDO
                        BE3 = MATMUL(APVtPA,MATMUL(IPP,MATMUL(RE1R,MATMUL(DE,AVtA))))
                        CE3 = MATMUL(AVtA,MATMUL(RE1R,MATMUL(DE,MATMUL(IPP,APVtPA))))
                        DE3 = MATMUL(APVtPA,MATMUL(IPP,MATMUL(RE1R,MATMUL(DE,MATMUL(IPP,APVtPA)))))
                        EE3 = MATMUL(AVtA,MATMUL(RE1R,MATMUL(DE,AVtA)))
                
                        TEMP =  0.50d0*(BE3 + CE3 - DE3 - EE3)
                        TERM(I,:) = TEMP(I,:)
                        TEMP(:,:) = 0.0d0
                ENDDO
        
                E4 = E4 + MATMUL(TERM,E1)

                TERM(:,:) = 0.0d0
                !----------------------
                ! Equation (11) page 24:
                !----------------------
                DO I=1,NB
                        DO J=1,NB
                                DE(J,J) = 1.0d0/( EP(J,J) + EP(I,I) )
                        ENDDO
                        BE3 = MATMUL(APVtPA,MATMUL(DE,AVtA))
                        CE3 = MATMUL(AVtA,MATMUL(DE,APVtPA))
                        DE3 = MATMUL(APVtPA,MATMUL(DE,MATMUL(IPP,APVtPA)))
                        EE3 = MATMUL(AVtA,MATMUL(DE,MATMUL(PP,AVtA)))
                
                        TEMP =  0.50d0*MATMUL(E1, (BE3 + CE3 - DE3 - EE3)  )
                        TERM(I,:) = TEMP(I,:)
                        TEMP(:,:) = 0.0d0
                ENDDO

                E4 = E4 - MATMUL(TERM,E1)
                
                !---------------------------------------------
                ! The term (d) on pages 25-25 in my DKH notes:
                !---------------------------------------------
                
                TERM(:,:) = 0.0d0
                !----------------------
                ! Equation (12) page 25:
                !----------------------
                DO I=1,NB
                        DO J=1,NB
                                DE(J,J) = 1.0d0/( EP(J,J) + EP(I,I) )
                        ENDDO
                        BE3 = MATMUL(APVtPA,MATMUL(DE,AVtA))
                        CE3 = MATMUL(AVtA,MATMUL(DE,APVtPA))
                        DE3 = MATMUL(APVtPA,MATMUL(DE,MATMUL(IPP,APVtPA)))
                        EE3 = MATMUL(AVtA,MATMUL(DE,MATMUL(PP,AVtA)))
                
                        TEMP =  0.50d0*MATMUL((BE3 + CE3 - DE3 - EE3), E1 )
                        TERM(:,I) = TEMP(:,I)
                        TEMP(:,:) = 0.0d0
                ENDDO

                E4 = E4 - MATMUL(E1,TERM)
                
                TERM(:,:) = 0.0d0
                !-----------------------
                ! Equation (13) page 25:
                !-----------------------
                DO I=1,NB
                        DO J=1,NB
                                DE(J,J) = 1.0d0/( EP(J,J) + EP(I,I) )
                        ENDDO
                        BE3 = MATMUL(APVtPA,MATMUL(IPP,MATMUL(DE,MATMUL(RE1R,AVtA))))
                        CE3 = MATMUL(AVtA,MATMUL(DE,MATMUL(RE1R,MATMUL(IPP,APVtPA))))
                        DE3 = MATMUL(APVtPA,MATMUL(IPP,MATMUL(DE,MATMUL(RE1R,MATMUL(IPP,APVtPA)))))
                        EE3 = MATMUL(AVtA,MATMUL(DE,MATMUL(RE1R,AVtA)))
                
                        TEMP =  0.50d0*(BE3 + CE3 - DE3 - EE3)
                        TERM(:,I) = TEMP(:,I)
                        TEMP(:,:) = 0.0d0
                ENDDO
        
                E4 = E4 + MATMUL(E1, TERM) 
                
                !---------------------------------------------
                ! The term (e) on pages 26-27 in my DKH notes:
                !---------------------------------------------
                
                TERM(:,:) = 0.0d0
                TEMP(:,:) = 0.0d0
                !-----------------------
                ! Equation (15) page 27:
                !-----------------------
                DO I=1,NB
                        DO J=1,NB
                                DE(J,J) = 1.0d0/( EP(J,J) + EP(I,I) )
                        ENDDO
                        TEMP = MATMUL(APVtPA,MATMUL(IPP,MATMUL(APVPA,MATMUL(DE,MATMUL(IPP,MATMUL(APVPA,AVtA)))))) &
                        &    + MATMUL(APVtPA,MATMUL(AVA,MATMUL(DE,MATMUL(PP,MATMUL(AVA,MATMUL(PP,AVtA)))))) &
                        &    + MATMUL(APVtPA,MATMUL(IPP,MATMUL(APVPA,MATMUL(DE,MATMUL(AVA,MATMUL(PP,AVtA)))))) &
                        &    + MATMUL(APVtPA,MATMUL(AVA,MATMUL(DE,MATMUL(APVPA,AVtA)))) &
                        &    + MATMUL(AVtA,MATMUL(APVPA,MATMUL(DE,MATMUL(IPP,MATMUL(APVPA,MATMUL(IPP,APVtPA)))))) &
                        &    + MATMUL(AVtA,MATMUL(PP,MATMUL(AVA,MATMUL(DE,MATMUL(PP,MATMUL(AVA,APVtPA)))))) &
                        &    + MATMUL(AVtA,MATMUL(APVPA,MATMUL(DE,MATMUL(AVA,APVtPA)))) &
                        &    + MATMUL(AVtA,MATMUL(PP,MATMUL(AVA,MATMUL(DE,MATMUL(APVPA,MATMUL(IPP,APVtPA)))))) &
                        &    - MATMUL(APVtPA,MATMUL(IPP,MATMUL(APVPA,MATMUL(DE,MATMUL(IPP,MATMUL(APVPA,MATMUL(IPP,APVtPA))))))) &
                        &    - MATMUL(APVtPA,MATMUL(AVA,MATMUL(DE,MATMUL(PP,MATMUL(AVA,APVtPA))))) &
                        &    - MATMUL(APVtPA,MATMUL(IPP,MATMUL(APVPA,MATMUL(DE,MATMUL(AVA,APVtPA))))) &
                        &    - MATMUL(APVtPA,MATMUL(AVA,MATMUL(DE,MATMUL(APVPA,MATMUL(IPP,APVtPA))))) &
                        &    - MATMUL(AVtA,MATMUL(APVPA,MATMUL(DE,MATMUL(IPP,MATMUL(APVPA,AVtA))))) & 
                        &    - MATMUL(AVtA,MATMUL(PP,MATMUL(AVA,MATMUL(DE,MATMUL(PP,MATMUL(AVA,MATMUL(PP,AVtA))))))) &
                        &    - MATMUL(AVtA,MATMUL(APVPA,MATMUL(DE,MATMUL(AVA,MATMUL(PP,AVtA))))) &
                        &    - MATMUL(AVtA,MATMUL(PP,MATMUL(AVA,MATMUL(DE,MATMUL(APVPA,AVtA)))))
                        TERM(I,:) = TEMP(I,:)
                        TEMP(:,:) = 0.0d0
                ENDDO
        
                E4 = E4 - 0.50d0*TERM

                TERM(:,:) = 0.0d0
                !-----------------------
                ! Equation (16) page 27:
                !-----------------------
                DO I=1,NB
                        DO J=1,NB
                                DE(J,J) = 1.0d0/( EP(J,J) + EP(I,I) )
                        ENDDO
                        BE3 = MATMUL(APVtPA,MATMUL(IPP,MATMUL(DE,MATMUL(RE1R,AVtA))))
                        CE3 = MATMUL(AVtA,MATMUL(DE,MATMUL(RE1R,MATMUL(IPP,APVtPA))))
                        DE3 = MATMUL(APVtPA,MATMUL(IPP,MATMUL(DE,MATMUL(RE1R,MATMUL(IPP,APVtPA)))))
                        EE3 = MATMUL(AVtA,MATMUL(DE,MATMUL(RE1R,AVtA)))
                
                        TEMP =  0.50d0*MATMUL(E1, (BE3 + CE3 - DE3 - EE3))
                        TERM(I,:) = TEMP(I,:)
                        TEMP(:,:) = 0.0d0
                ENDDO
        
                E4 = E4 + TERM
                
                !---------------------------------------------
                ! The term (f) on pages 26-28 in my DKH notes:
                !---------------------------------------------
               
                TEMP(:,:) = 0.0d0 
                !-----------------------
                ! Equation (17) page 27:
                !-----------------------
                DO I=1,NB
                        DO J=1,NB
                                DE(J,J) = 1.0d0/( EP(J,J) + EP(I,I) )
                        ENDDO
                        BE3 = MATMUL(APVtPA,MATMUL(IPP,MATMUL(DE,MATMUL(RE1R,AVtA))))
                        CE3 = MATMUL(AVtA,MATMUL(DE,MATMUL(RE1R,MATMUL(IPP,APVtPA))))
                        DE3 = MATMUL(APVtPA,MATMUL(IPP,MATMUL(DE,MATMUL(RE1R,MATMUL(IPP,APVtPA)))))
                        EE3 = MATMUL(AVtA,MATMUL(DE,MATMUL(RE1R,AVtA)))
                
                        TEMP =  0.50d0*MATMUL((BE3 + CE3 - DE3 - EE3),E1)
                        TERM(:,I) = TEMP(:,I)
                        TEMP(:,:) = 0.0d0
                ENDDO
        
                E4 = E4 + TERM
                
                TERM(:,:) = 0.0d0
                TEMP(:,:) = 0.0d0
                !-----------------------
                ! Equation (15) page 27:
                !-----------------------
                DO I=1,NB
                        DO J=1,NB
                                DE(J,J) = 1.0d0/( EP(J,J) + EP(I,I) )
                        ENDDO
                        TEMP = MATMUL(APVtPA,MATMUL(IPP,MATMUL(APVPA,MATMUL(DE,MATMUL(IPP,MATMUL(APVPA,AVtA)))))) &
                        &    + MATMUL(APVtPA,MATMUL(AVA,MATMUL(DE,MATMUL(PP,MATMUL(AVA,MATMUL(PP,AVtA)))))) &
                        &    + MATMUL(APVtPA,MATMUL(IPP,MATMUL(APVPA,MATMUL(DE,MATMUL(AVA,MATMUL(PP,AVtA)))))) &
                        &    + MATMUL(APVtPA,MATMUL(AVA,MATMUL(DE,MATMUL(APVPA,AVtA)))) &
                        &    + MATMUL(AVtA,MATMUL(APVPA,MATMUL(DE,MATMUL(IPP,MATMUL(APVPA,MATMUL(IPP,APVtPA)))))) &
                        &    + MATMUL(AVtA,MATMUL(PP,MATMUL(AVA,MATMUL(DE,MATMUL(PP,MATMUL(AVA,APVtPA)))))) &
                        &    + MATMUL(AVtA,MATMUL(APVPA,MATMUL(DE,MATMUL(AVA,APVtPA)))) &
                        &    + MATMUL(AVtA,MATMUL(PP,MATMUL(AVA,MATMUL(DE,MATMUL(APVPA,MATMUL(IPP,APVtPA)))))) &
                        &    - MATMUL(APVtPA,MATMUL(IPP,MATMUL(APVPA,MATMUL(DE,MATMUL(IPP,MATMUL(APVPA,MATMUL(IPP,APVtPA))))))) &
                        &    - MATMUL(APVtPA,MATMUL(AVA,MATMUL(DE,MATMUL(PP,MATMUL(AVA,APVtPA))))) &
                        &    - MATMUL(APVtPA,MATMUL(IPP,MATMUL(APVPA,MATMUL(DE,MATMUL(AVA,APVtPA))))) &
                        &    - MATMUL(APVtPA,MATMUL(AVA,MATMUL(DE,MATMUL(APVPA,MATMUL(IPP,APVtPA))))) &
                        &    - MATMUL(AVtA,MATMUL(APVPA,MATMUL(DE,MATMUL(IPP,MATMUL(APVPA,AVtA))))) & 
                        &    - MATMUL(AVtA,MATMUL(PP,MATMUL(AVA,MATMUL(DE,MATMUL(PP,MATMUL(AVA,MATMUL(PP,AVtA))))))) &
                        &    - MATMUL(AVtA,MATMUL(APVPA,MATMUL(DE,MATMUL(AVA,MATMUL(PP,AVtA))))) &
                        &    - MATMUL(AVtA,MATMUL(PP,MATMUL(AVA,MATMUL(DE,MATMUL(APVPA,AVtA)))))
                        TERM(:,I) = TEMP(:,I)
                        TEMP(:,:) = 0.0d0
                ENDDO
        
                E4 = E4 - 0.50d0*TERM

                HDKH = HDKH + E4
        ENDIF
        
        IF (DKHORDER .GE. 5 ) THEN
                E5 = 0.0d0
                HDKH = HDKH + E5
        ENDIF

        !-----------------------------------------------------------------
        ! Transforming the HDKH hamiltonian from the p**2 -representation
        ! back to the strandard representation of uquantchem
        !-----------------------------------------------------------------

         HDKH= MATMUL(OMEGA,MATMUL(HDKH,TRANSPOSE(OMEGA)))

        !-------------------------------------------------------
        ! Transforming back to the non-orthogonal representation
        !-------------------------------------------------------

         HDKH = MATMUL(SH,MATMUL(HDKH,SH))

      ELSE
        HDKH = T + V
      ENDIF
END SUBROUTINE DKH

