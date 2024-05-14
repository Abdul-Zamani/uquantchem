SUBROUTINE EP3r(MULTIPLICITY,Cup,Cdown,Ints,NB,Ne,EHFeigenup,EHFeigendown,E0,nuce,SPINCONSERVE)
      ! 
      !


      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NB,Ne,MULTIPLICITY
      DOUBLE PRECISION, INTENT(IN) :: Cup(NB,NB),Cdown(NB,NB),Ints(NB,NB,NB,NB)
      DOUBLE PRECISION, INTENT(IN) :: EHFeigenup(NB),EHFeigendown(NB), E0,nuce
      LOGICAL, INTENT(IN) :: SPINCONSERVE
      INTEGER :: I,J,N,M,O,P,N2,K,L,K1,II,JJ,DIFFER,Q,Nel,NBl
      INTEGER, ALLOCATABLE :: EXC(:,:),V1(:)
      DOUBLE PRECISION, ALLOCATABLE :: Calpha(:),Ca(:),Cbeta(:),Cb(:),EP2TEMP(:)
      DOUBLE PRECISION, ALLOCATABLE :: moIntsAA(:,:,:,:),&
                                       moIntsAB(:,:,:,:),&
                                       moIntsBA(:,:,:,:),&
                                       moIntsBB(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: tempInts1(:,:,:,:), & 
                                       tempInts2(:,:,:,:), &
                                       tempInts3(:,:,:,:)
      DOUBLE PRECISION :: TEMP1,DEIJKL
      LOGICAL :: EXCITE, conver,osd2

      DOUBLE PRECISION :: Sz,twoSP1,Nalpha,Nbeta,SEOld1,SEOld2,&
                          SEOld1AA,SEOld1AB,SEOld2AA,SEOld2AB,numer,denom,&
                          R2,R2AA,R2BB,R2AB,R2BA,R3,R3CS,R3CD,&
                          EMP2,EMP2AA,EMP2AB,EMP2BA,EMP2BB,EPole,EPoleOld,E,PS
      INTEGER :: a,b,c,d,pole,mu,nu,lam,sig,neup,nedown,iter
     
          
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
          
      allocate(tempInts1(nb,nb,nb,nb))
      allocate(tempInts2(nb,nb,nb,nb))
      allocate(tempInts3(nb,nb,nb,nb))

      allocate(MOIntsAA(nb,nb,nb,nb))
      allocate(MOIntsAB(nb,nb,nb,nb))
      allocate(MOIntsBA(nb,nb,nb,nb))
      allocate(MOIntsBB(nb,nb,nb,nb))

      !Quarter transformations: O(N^5)
        print*,' '
        print*,' '
        print*,'==========================================================='
        print*,'           Enter Quarter transformations: O(N^5)           '
        print*,'==========================================================='
        print*,' '
        print*,' '
        print*,' '

!AA
        print*,'   AA   '        
      tempInts1=0.0d0 
      tempInts2=0.0d0
      tempInts3=0.0d0
      moIntsAA = 0.0d0
      do i=1, nb
        do mu=1, nb
          tempInts1(i,:,:,:) = tempInts1(i,:,:,:) + &
          (Cup(mu,i)*Ints(mu,:,:,:))
        enddo 
        do j=1, nb
          do nu=1, nb
            tempInts2(i,j,:,:) = tempInts2(i,j,:,:) + &
            (Cup(nu,j)*tempInts1(i,nu,:,:))
           enddo
          do k=1, nb
            do lam=1, nb
              tempInts3(i,j,k,:) = tempInts3(i,j,k,:) + &
              (Cup(lam,k)*tempInts2(i,j,lam,:))
            enddo
            do l=1, nb
              do sig=1, nb
                MOIntsAA(i,j,k,l) = MOIntsAA(i,j,k,l) + &
                (Cup(sig,l)*tempInts3(i,j,k,sig))
              enddo      
            enddo!end i            
          enddo!end j     
        enddo!end k         
      enddo!end l   
 !AB
        print*,'   AB   '        
      tempInts1=0.0d0 
      tempInts2=0.0d0
      tempInts3=0.0d0
      moIntsAB = 0.0d0
      do i=1, nb
        do mu=1, nb
          tempInts1(i,:,:,:) = tempInts1(i,:,:,:) + &
          (Cup(mu,i)*Ints(mu,:,:,:))
        enddo 
        do j=1, nb
          do nu=1, nb
            tempInts2(i,j,:,:) = tempInts2(i,j,:,:) + &
            (Cup(nu,j)*tempInts1(i,nu,:,:))
           enddo
          do k=1, nb
            do lam=1, nb
              tempInts3(i,j,k,:) = tempInts3(i,j,k,:) + &
              (Cdown(lam,k)*tempInts2(i,j,lam,:))
            enddo
            do l=1, nb
              do sig=1, nb
                MOIntsAB(i,j,k,l) = MOIntsAB(i,j,k,l) + &
                (Cdown(sig,l)*tempInts3(i,j,k,sig))
              enddo      
            enddo!end i            
          enddo!end j     
        enddo!end k         
      enddo!end l   

!AZ no BA or BB for closed shell
 
      print*,'   Tran.Done.   '        
    
      deallocate(tempInts1,tempInts2,tempInts3)


        print*,' '
        print*,' '
        print*,'==========================================================='
        print*,'            Results from the EP2 calculation:              '
        print*,'==========================================================='
        print*,' '              
        print*,' '

        do pole=1,neUp !!!!! begin pole search

        print*,' '
        print*,'orb',pole
        print*,' '

        numer= 0.0d0
        denom= 0.0d0
        R2   = 0.0d0
        R2AA = 0.0d0
        R2BB = 0.0d0
        R2AB = 0.0d0
        R2BA = 0.0d0
        R3   = 0.0d0 
        R3CS = 0.0d0
        R3CD = 0.0d0 

        if(pole.le.nEup) then !ip relax  
!IPs 
!AZ 3rd order AA
        !term 1 AA
        do i=1,Neup
          do j=1,NeUp
            do a=NeUp+1,NB
              do b=NeUp+1,NB
                do c=NeUp+1,NB
                  numer = (MOIntsAA(pole,pole,b,i)-MOIntsAA(pole,i,b,pole))*&
                          (MOIntsAA(i,a,j,c)-MOIntsAA(i,c,j,a))*&
                          (MOIntsAA(a,b,c,j)-MOIntsAA(a,j,c,b)) 
                  denom = (EHFeigenup(i)-EHFeigenup(b))*&
                          (EHFeigenup(i)+EHFeigenup(j)-&
                          EHFeigenup(a)-EHFeigenup(c)) 
                  R2AA  =  R2AA + (numer/denom)*(0.5d0) !+1/2 factor
                  R3CS  =  R3CS + (numer/denom)*(0.5d0)
                enddo 
              enddo 
            enddo 
          enddo
        enddo
        !term 2 AA
        do i=1,Neup
          do j=1,NeUp
            do k=1,NeUp
              do a=NeUp+1,NB
                do b=NeUp+1,NB
                  numer = (MOIntsAA(pole,pole,a,j)-MOIntsAA(pole,j,a,pole))*&
                          (MOIntsAA(i,a,k,b)-MOIntsAA(i,b,k,a))*&
                          (MOIntsAA(j,i,b,k)-MOIntsAA(j,k,b,i)) 
                  denom = (EHFeigenup(j)-EHFeigenup(a))*&
                          (EHFeigenup(i)+EHFeigenup(k)-&
                          EHFeigenup(a)-EHFeigenup(b)) 
                  R2AA  =  R2AA + (numer/denom)*(-0.5d0) !-1/2 factor
                  R3CS  =  R3CS + (numer/denom)*(-0.5d0)
                enddo 
              enddo 
            enddo 
          enddo
        enddo
        !term 3 AA
        do i=1,Neup
          do j=1,NeUp
            do a=NeUp+1,NB
              do b=NeUp+1,NB
                do c=NeUp+1,NB
                  numer = (MOIntsAA(pole,pole,a,b)-MOIntsAA(pole,b,a,pole))*&
                          (MOIntsAA(i,a,j,c)-MOIntsAA(i,c,j,a))*&
                          (MOIntsAA(b,i,c,j)-MOIntsAA(b,j,c,i)) 
                  denom = (EHFeigenup(i)+EHFeigenup(j)-&
                          EHFeigenup(a)-EHFeigenup(c))*&
                          (EHFeigenup(i)+EHFeigenup(j)-&
                          EHFeigenup(b)-EHFeigenup(c)) 
                  R2AA  =  R2AA + (numer/denom)*(0.5d0) !+1/2 factor
                  R3CD =   R3CD + (numer/denom)*(0.5d0) 
                enddo 
              enddo 
            enddo 
          enddo
        enddo     
        !term 4 AA
        do i=1,Neup
          do j=1,NeUp
            do k=1,NeUp
              do a=NeUp+1,NB
                do b=NeUp+1,NB
                  numer = (MOIntsAA(pole,pole,j,i)-MOIntsAA(pole,i,j,pole))*&
                          (MOIntsAA(i,a,k,b)-MOIntsAA(i,b,k,a))*&
                          (MOIntsAA(a,j,b,k)-MOIntsAA(a,k,b,j)) 
                  denom = (EHFeigenup(j)+EHFeigenup(k)-&
                          EHFeigenup(a)-EHFeigenup(b))*&
                          (EHFeigenup(i)+EHFeigenup(k)-&
                          EHFeigenup(a)-EHFeigenup(b)) 
                  R2AA  =  R2AA + (numer/denom)*(-0.5d0) !-1/2 factor
                  R3CD =   R3CD + (numer/denom)*(-0.5d0)
                enddo 
              enddo 
            enddo 
          enddo
        enddo     
        !term 5 AA
        do i=1,Neup
          do j=1,NeUp
            do a=NeUp+1,NB
              do b=NeUp+1,NB
                do c=NeUp+1,NB
                  numer = (MOIntsAA(pole,pole,i,a)-MOIntsAA(pole,a,i,pole))*&
                          (MOIntsAA(b,i,c,j)-MOIntsAA(b,j,c,i))*&
                          (MOIntsAA(a,b,j,c)-MOIntsAA(a,c,j,b)) 
                  denom = (EHFeigenup(i)-EHFeigenup(a))*&
                          (EHFeigenup(i)+EHFeigenup(j)-&
                          EHFeigenup(b)-EHFeigenup(c)) 
                  R2AA  =  R2AA + (numer/denom)*(0.5d0) !+1/2 factor
                  R3CS  =  R3CS + (numer/denom)*(0.5d0)
                enddo 
              enddo 
            enddo 
          enddo
        enddo     
        !term 6 AA
        do i=1,Neup
          do j=1,NeUp
            do k=1,NeUp
              do a=NeUp+1,NB
                do b=NeUp+1,NB
                  numer = (MOIntsAA(pole,pole,i,a)-MOIntsAA(pole,a,i,pole))*&
                          (MOIntsAA(a,j,b,k)-MOIntsAA(a,k,b,j))*&
                          (MOIntsAA(j,i,k,b)-MOIntsAA(j,b,k,i)) 
                  denom = (EHFeigenup(i)-EHFeigenup(a))*&
                          (EHFeigenup(j)+EHFeigenup(k)-&
                          EHFeigenup(a)-EHFeigenup(b)) 
                  R2AA  =  R2AA + (numer/denom)*(-0.5d0) !-1/2 factor
                  R3CS  =  R3CS + (numer/denom)*(-0.5d0)
                enddo 
              enddo 
            enddo 
          enddo
        enddo

!AZ save 3rd order relax into r3
        R3 = R2AA
!AZ 3rd order AA

!2nd order AA
!AA aa||aa  pp||ai 
        do i=1,Neup
          do a=NeUp+1,NB
            denom = EHFeigenup(a) - EHFeigenup(i)
!            numer = MOIntsAA(a,i,pole,pole) - MOIntsAA(a,pole,pole,i)
            numer = MOIntsAA(pole,pole,a,i) - MOIntsAA(pole,i,a,pole)


            numer = numer*numer 
            R2AA  =  R2AA + numer/denom 
          enddo
        enddo



       R2BB = R2AA


!AZ 3rd order AB !1/12 start here
        !term 1 AB; just make them closed shell, alpha=beta 
        do i=1,Neup
          do j=1,NeUp
            do a=NeUp+1,NB
              do b=NeUp+1,NB
                do c=NeUp+1,NB
                  numer = (MOIntsAB(pole,pole,b,i))*&
                          (MOIntsAB(i,a,j,c))*&
                          (MOIntsAB(a,b,c,j)) 
                  denom = (EHFeigenup(i)-EHFeigenup(b))*&
                          (EHFeigenup(i)+EHFeigenup(j)-&
                          EHFeigenup(a)-EHFeigenup(c)) 
                  R2AB  =  R2AB + (numer/denom)*(0.5d0) !+1/2 factor
                  R3CS  =  R3CS + (numer/denom)*(0.5d0)
                enddo 
              enddo 
            enddo 
          enddo
        enddo
        !term 2 AB
        do i=1,Neup
          do j=1,NeUp
            do k=1,NeUp
              do a=NeUp+1,NB
                do b=NeUp+1,NB
                  numer = (MOIntsAB(pole,pole,a,j))*&
                          (MOIntsAB(i,a,k,b))*&
                          (MOIntsAB(j,i,b,k)) 
                  denom = (EHFeigenup(j)-EHFeigenup(a))*&
                          (EHFeigenup(i)+EHFeigenup(k)-&
                          EHFeigenup(a)-EHFeigenup(b)) 
                  R2AB  =  R2AB + (numer/denom)*(-0.5d0) !-1/2 factor
                  R3CS  =  R3CS + (numer/denom)*(-0.5d0)
                enddo 
              enddo 
            enddo 
          enddo
        enddo
        !term 3 AB
        do i=1,Neup
          do j=1,NeUp
            do a=NeUp+1,NB
              do b=NeUp+1,NB
                do c=NeUp+1,NB
                  numer = (MOIntsAB(pole,pole,a,b))*&
                          (MOIntsAB(i,a,j,c))*&
                          (MOIntsAB(b,i,c,j)) 
                  denom = (EHFeigenup(i)+EHFeigenup(j)-&
                          EHFeigenup(a)-EHFeigenup(c))*&
                          (EHFeigenup(i)+EHFeigenup(j)-&
                          EHFeigenup(b)-EHFeigenup(c)) 
                  R2AB  =  R2AB + (numer/denom)*(0.5d0) !+1/2 factor
                  R3CD  =  R3CD + (numer/denom)*(0.5d0)
                enddo 
              enddo 
            enddo 
          enddo
        enddo     
        !term 4 AB
        do i=1,Neup
          do j=1,NeUp
            do k=1,NeUp
              do a=NeUp+1,NB
                do b=NeUp+1,NB
                  numer = (MOIntsAA(pole,pole,j,i))*&
                          (MOIntsAA(i,a,k,b))*&
                          (MOIntsAA(a,j,b,k)) 
                  denom = (EHFeigenup(j)+EHFeigenup(k)-&
                          EHFeigenup(a)-EHFeigenup(b))*&
                          (EHFeigenup(i)+EHFeigenup(k)-&
                          EHFeigenup(a)-EHFeigenup(b)) 
                  R2AB  =  R2AB + (numer/denom)*(-0.5d0) !-1/2 factor
                  R3CD  =  R3CD + (numer/denom)*(-0.5d0)
                enddo 
              enddo 
            enddo 
          enddo
        enddo     
        !term 5 AB
        do i=1,Neup
          do j=1,NeUp
            do a=NeUp+1,NB
              do b=NeUp+1,NB
                do c=NeUp+1,NB
                  numer = (MOIntsAB(pole,pole,i,a))*&
                          (MOIntsAB(b,i,c,j))*&
                          (MOIntsAB(a,b,j,c)) 
                  denom = (EHFeigenup(i)-EHFeigenup(a))*&
                          (EHFeigenup(i)+EHFeigenup(j)-&
                          EHFeigenup(b)-EHFeigenup(c)) 
                  R2AB  =  R2AB + (numer/denom)*(0.5d0) !+1/2 factor
                  R3CS  =  R3CS + (numer/denom)*(0.5d0)
                enddo 
              enddo 
            enddo 
          enddo
        enddo     
        !term 6 AB
        do i=1,Neup
          do j=1,NeUp
            do k=1,NeUp
              do a=NeUp+1,NB
                do b=NeUp+1,NB
                  numer = (MOIntsAB(pole,pole,i,a))*&
                          (MOIntsAB(a,j,b,k))*&
                          (MOIntsAB(j,i,k,b)) 
                  denom = (EHFeigenup(i)-EHFeigenup(a))*&
                          (EHFeigenup(j)+EHFeigenup(k)-&
                          EHFeigenup(a)-EHFeigenup(b)) 
                  R2AB  =  R2AB + (numer/denom)*(-0.5d0) !-1/2 factor
                  R3CS  =  R3CS + (numer/denom)*(-0.5d0)
                enddo 
              enddo 
            enddo 
          enddo
        enddo
 
!AZ 3rd order AA

!AZ save r3 correction

        R3=R3+R2AB

!2nd order AB 
!AB aa||bb pp||ai  
        do i=1,NeDown
          do a=NeDown+1,NB 
            denom = EHFeigendown(a) - EHFeigendown(i)
!            numer = MOIntsAB(a,i,pole,pole)
            numer = MOIntsAB(pole,pole,a,i)
            numer = numer*numer 
            R2AB  =  R2AB + numer/denom 
           enddo
         enddo 


        elseif(pole.gt.neUp) then !ea relax
!EAs
!AA aa||aa  pp||ai 
        do i=1,Neup
          do a=NeUp+1,NB
            denom = EHFeigenup(i) - EHFeigenup(a)
!            numer = MOIntsAA(a,i,pole,pole) - MOIntsAA(a,pole,pole,i)
            numer = MOIntsAA(pole,pole,i,a) - MOIntsAA(pole,a,i,pole)


            numer = numer*numer 
            R2AA  =  R2AA + numer/denom 
          enddo
        enddo

       R2BB = R2AA

!AB aa||bb pp||ai  
        do i=1,NeDown
          do a=NeDown+1,NB 
            denom = EHFeigendown(i) - EHFeigendown(a)
!            numer = MOIntsAB(a,i,pole,pole)
            numer = MOIntsAB(pole,pole,i,a)
            numer = numer*numer 
            R2AB  =  R2AB + numer/denom 
           enddo
         enddo 


        endif!end relax 

       R2BA = R2AB

       R2 = R2AA + R2AB 

       print*,'Koopmans =',EHFeigenup(pole)
       print*,'R3 only (Ha)',R3
       print*,'R3 CS:  (Ha)',R3CS
       print*,'R3 CD:  (Ha)',R3CD     
       print*,'R2+R3 (Ha) =',R2
       print*,'R2+R3 AA (Ha) =',R2AA
       print*,'R2+R3 AB (Ha) =',R2AB
!       print*,'eps+R2 (Ha) =',((0.0d0-EHFeigenup(pole))+R2)
!       print*,'eps+R2 (eV) =',((0.0d0-EHFeigenup(pole))+R2)*27.2114
       print*,'eps+R3 (Ha) =',((0.0d0-EHFeigenup(pole))-R2)
       print*,'eps+R3 (eV) =',((0.0d0-EHFeigenup(pole))-R2)*27.2114
       enddo !!!!! end pole search

       print*,' '
       print*,' '
       print*,'==========================================================='
       print*,'            Results from the MP2 calculation:              '
       print*,'==========================================================='
       print*,' '

!
       emp2=0.0d0
       emp2aa=0.0d0
       emp2ab=0.0d0
       emp2ba=0.0d0
       emp2bb=0.0d0 
        do i=1,neup
          do j=1,neup
            do a=neup+1,nB
              do b=neup+1,nB
              EMP2AA = EMP2AA +  &
                (moIntsAA(i,a,j,b)-moIntsAA(i,b,j,a))**2 / &  
                ((ehfeigenup(i) + ehfeigenup(j)) - &
                ehfeigenup(a) - ehfeigenup(b))
!OS K int is 0 no double counting  
              EMP2AB = EMP2AB +  &
                (moIntsAB(i,a,j,b))**2 / &  
                ((ehfeigenup(i) + ehfeigendown(j)) - &
                ehfeigenup(a) - ehfeigendown(b))

              EMP2BA = EMP2BA +  &
                (moIntsBA(i,a,j,b))**2 / &  
                ((ehfeigendown(i) + ehfeigenup(j)) - &
                ehfeigendown(a) - ehfeigenup(b))

!              EMP2BB = EMP2BB +  &
!                (moIntsBB(i,a,j,b)-moIntsBB(i,b,j,a))**2 / &  
!                ((ehfeigendown(i) + ehfeigendown(j)) - &
!                ehfeigendown(a) - ehfeigendown(b))
               
              enddo 
            enddo 
          enddo 
        enddo 

        EMP2BB = EMP2AA 

        print*,' '
        EMP2AA =  EMP2AA/4.0d0
        EMP2AB =  EMP2AB
        EMP2BA =  EMP2BA
        EMP2BB =  EMP2BB/4.0d0
        print*,'E2AA (Ha) =',EMP2AA
        print*,'E2AB (Ha) =',EMP2AB
        print*,'E2BB (Ha) =',EMP2BB
        
        EMP2=EMP2AA+EMP2AB+EMP2BB

        print*,'E2 (Ha) =',EMP2

        print*,' '
 
        WRITE(*,'(A27,F30.20,A3)')'      Correlation Energy =',EMP2,' au'
        WRITE(*,'(A27,F30.20,A3)')'     Hartree-Fock Energy =',E0,' au'
        IF ( nuce .NE. 0) WRITE(*,'(A27,F30.20,A3)')' Nuclear repulsion Energy =',nuce,' au'
        IF ( nuce .NE. 0) WRITE(*,'(A27,F30.20,A3)')'            Total Energy =',E0+EMP2+nuce,' au'
        IF ( nuce .EQ. 0) WRITE(*,'(A27,F30.20,A3)')'            Total Energy =',E0+EMP2,' au'
        print*,' '
        print*,' '

        deallocate(MOIntsAA,MOIntsAB,MOIntsBA,MOIntsBB)
 
        !-----------------------------
        ! HERE THE OUTPUT IS GENERATED
        !-----------------------------
!        WRITE(*,'(A27,F30.20,A3)')'      Correlation Energy =',EP2,' au'
!        WRITE(*,'(A27,F30.20,A3)')'     Hartree-Fock Energy =',E0,' au'
!        print*,' '
!        print*,' '
       
        END SUBROUTINE EP3r
