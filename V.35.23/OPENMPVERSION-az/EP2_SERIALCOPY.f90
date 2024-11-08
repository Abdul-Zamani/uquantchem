SUBROUTINE EP2(MULTIPLICITY,Cup,Cdown,Ints,NB,Ne,EHFeigenup,EHFeigendown,E0,nuce,SPINCONSERVE)
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
                                       tempInts3(:,:,:,:),&
                                       tijab(:,:,:,:),&!11/6/24
                                       pocc(:,:),&!11/6/24
                                       pvirt(:,:),&!11/6/24
                                       tempVec(:),&!11/7/24
                                       tempMat(:,:)!11/7/24
      DOUBLE PRECISION :: TEMP1,DEIJKL
      LOGICAL :: EXCITE, conver,osd2

      DOUBLE PRECISION :: Sz,twoSP1,Nalpha,Nbeta,SEOld1,SEOld2,&
                          SEOld1AA,SEOld1AB,SEOld2AA,SEOld2AB,&
                          EMP2,EMP2AA,EMP2AB,EMP2BA,EMP2BB,EPole,EPoleOld,E,PS
!AZ 10/3/24
     DOUBLE PRECISION :: trPOcc,trPVirt
!AZ 10/3/24
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

!AZ 11/6/24 
      allocate(tijab(nb,nb,nb,nb))
      allocate(pOcc(neup,neup))
      allocate(pVirt(nb-neup,nb-neup))
      allocate(tempVec(nb))
      allocate(tempMat(nb,nb))
!AZ 11/6/24 

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
            (Cdown(nu,j)*tempInts1(i,nu,:,:))
           enddo
          do k=1, nb
            do lam=1, nb
              tempInts3(i,j,k,:) = tempInts3(i,j,k,:) + &
              (Cup(lam,k)*tempInts2(i,j,lam,:))
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
 !BA
        print*,'   BA   '        
      tempInts1=0.0d0 
      tempInts2=0.0d0
      tempInts3=0.0d0
      moIntsBA = 0.0d0
      do i=1, nb
        do mu=1, nb
          tempInts1(i,:,:,:) = tempInts1(i,:,:,:) + &
          (Cdown(mu,i)*Ints(mu,:,:,:))
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
                MOIntsBA(i,j,k,l) = MOIntsBA(i,j,k,l) + &
                (Cup(sig,l)*tempInts3(i,j,k,sig))
              enddo      
            enddo!end i            
          enddo!end j     
        enddo!end k         
      enddo!end l   
  !BB
        print*,'   BB   '        
      tempInts1=0.0d0 
      tempInts2=0.0d0
      tempInts3=0.0d0
      moIntsBB = 0.0d0
      do i=1, nb
        do mu=1, nb
          tempInts1(i,:,:,:) = tempInts1(i,:,:,:) + &
          (Cdown(mu,i)*Ints(mu,:,:,:))
        enddo 
        do j=1, nb
          do nu=1, nb
            tempInts2(i,j,:,:) = tempInts2(i,j,:,:) + &
            (Cdown(nu,j)*tempInts1(i,nu,:,:))
           enddo
          do k=1, nb
            do lam=1, nb
              tempInts3(i,j,k,:) = tempInts3(i,j,k,:) + &
              (Cdown(lam,k)*tempInts2(i,j,lam,:))
            enddo
            do l=1, nb
              do sig=1, nb
                MOIntsBB(i,j,k,l) = MOIntsBB(i,j,k,l) + &
                (Cdown(sig,l)*tempInts3(i,j,k,sig))
              enddo      
            enddo!end i            
          enddo!end j     
        enddo!end k         
      enddo!end l   
 
      print*,'   Tran.Done.   '        
    
      deallocate(tempInts1,tempInts2,tempInts3)


        print*,' '
        print*,' '
        print*,'==========================================================='
        print*,'            Results from the EP2 calculation:              '
        print*,'==========================================================='
        print*,' '              
        print*,' '

        

        osd2 = .false. !flag for OS

        print*,'OSD2 is',osd2 

        do pole=1,NB !!!!! begin pole search

        iter=0 !max iter 15 for now

        print*,' '
        print*,'orb',pole
        print*,' '
        
        E=0.0d0
        SEOld1AA = 0.0d0
        SEOld1AB = 0.0d0
        SEOld2AA = 0.0d0
        SEOld2AB = 0.0d0


        SEold1 = 0.0d0
        SEold2 = 0.0d0
        PS=0.0d0
        E = EHFeigenup(pole)*0.92
        conver = .false. 

        EPoleOld=E
        do while(conver.eqv..false.)
       !! print*,'conver',conver 
!SE1  
        if(osd2.eqv..false.) then !!os
!AA
        do i=1,Neup
          do a=NeUp+1,NB 
            do b=NeUp+1,NB
!antisymm?
               SEold1AA = SEold1AA +( ((moIntsAA(pole,a,i,b)-moIntsAA(pole,b,i,a))**2.0d0) / &
               (EPoleOld + EHFeigenup(i) -EHFeigenup(a)-EHFeigenup(b)) )
            enddo
          enddo 
        enddo

       endif !!os
!AB        
        do i=1,Nedown
          do a=NeUp+1,NB 
            do b=Nedown+1,NB
               SEold1AB = SEold1AB + ( ((moIntsAB(pole,a,i,b))**2.0d0) / &


               (EPoleOld + EHFeigendown(i) -EHFeigenup(a)-EHFeigendown(b)))
            enddo
          enddo 
        enddo

!SE2
        if(osd2.eqv..false.) then !!os
!AA
        do a=Neup+1,Nb
          do i=1,NeUp 
            do j=1,NeUp
               SEold2AA = SEold2AA + ( ((moIntsAA(pole,i,a,j)-moIntsAA(pole,j,a,i))**2.0d0) / &
               (EPoleOld + EHFeigenup(a) -EHFeigenup(i)-EHFeigenup(j)) )
            enddo
          enddo 
        enddo
     
        endif !!os 
!AB        
        do a=Nedown+1,NB
          do i=1,NeUp 
            do j=1,Nedown
               SEold2AB = SEold2AB + ( ((moIntsAB(pole,i,a,j))**2.0d0) / &


               (EPoleOld + EHFeigendown(a) -EHFeigenup(i)-EHFeigendown(j)))
            enddo
          enddo 
        enddo

!new NR instead
       SEOld1 = SEOld1AA/2.0d0 +  SEOld1AB
       SEOld2 = SEOld2AA/2.0d0 + SEOld2AB
       EPole = EHFeigenup(pole) + SEOld1+SEOld2
!!       print*,'sigma(2)',SEOld1+SEOld2
       print*,'SEOld1',SEOld1
       print*,'SEOld2',SEOld2
         

!derivatives
        SEOld1AA = 0.0d0
        SEOld1AB = 0.0d0
        SEOld2AA = 0.0d0
        SEOld2AB = 0.0d0

       SEOld1=0.0d0
       SEOld2=0.0d0
!SE1  
        if(osd2.eqv..false.) then !!os
!AA
        do i=1,Neup
          do a=NeUp+1,NB 
            do b=NeUp+1,NB
!antisymm?
               SEold1AA = SEold1AA +( ((moIntsAA(pole,a,i,b)-moIntsAA(pole,b,i,a))**2.0d0) / &
               (EPoleOld + EHFeigenup(i) -EHFeigenup(a)-EHFeigenup(b))**2 )
            enddo
          enddo 
        enddo

        endif !!os 

!AB        
        do i=1,Nedown
          do a=NeUp+1,NB 
            do b=Nedown+1,NB
               SEold1AB = SEold1AB + ( ((moIntsAB(pole,a,i,b))**2.0d0) / &


               (EPoleOld + EHFeigendown(i) -EHFeigenup(a)-EHFeigendown(b))**2)
            enddo
          enddo 
        enddo

!SE2
        if(osd2.eqv..false.) then !!os
!AA
        do a=Neup+1,Nb
          do i=1,NeUp 
            do j=1,NeUp
               SEold2AA = SEold2AA + ( ((moIntsAA(pole,i,a,j)-moIntsAA(pole,j,a,i))**2.0d0) / &
               (EPoleOld+ EHFeigenup(a) -EHFeigenup(i)-EHFeigenup(j))**2 )
            enddo
          enddo 
        enddo

        endif !!os 

!AB        
        do a=Nedown+1,NB
          do i=1,NeUp 
            do j=1,Nedown
               SEold2AB = SEold2AB + ( ((moIntsAB(pole,i,a,j))**2.0d0) / &


               (EPoleOld + EHFeigendown(a) -EHFeigenup(i)-EHFeigendown(j))**2)
            enddo
          enddo 
        enddo

       SEold1 = -1.0d0*(SEold1AA/2.0d0  +  SEOld1AB)
       SEold2 = -1.0d0*(SEold2AA/2.0d0  +  SEOld2AB)

       E = (EpoleOld - ((EpoleOld-Epole)/(1-(SEold1+SEold2))))
       PS = 1/(1-(SEold1+SEold2))
       !print*,'E after NRstep',E      

       iter=iter+1

       if(abs(EPole-EPoleOld).lt.0.0001.or.iter.eq.15) then
       print*,'Koopmans =',EHFeigenup(pole)
       print*,'D2 (Ha) =',E
       print*,'D2 (eV) =',E*27.2114
       print*,'PS =',PS
       conver=.true.
       if(iter.eq.15) then 
         print*,'pole not converged after',iter,'iter'
       endif 
       endif 

       EPoleOld=E
       SEOld1=0.0d0
       SEOld2=0.0d0
        SEOld1AA = 0.0d0
        SEOld1AB = 0.0d0
        SEOld2AA = 0.0d0
        SEOld2AB = 0.0d0


       enddo!while

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

              EMP2BB = EMP2BB +  &
                (moIntsBB(i,a,j,b)-moIntsBB(i,b,j,a))**2 / &  
                ((ehfeigendown(i) + ehfeigendown(j)) - &
                ehfeigendown(a) - ehfeigendown(b))
 
              !AZ 11/6/24
              tijab(i,a,j,b) = (moIntsAA(i,a,j,b)-moIntsAA(i,b,j,a)) / &
              ((ehfeigenup(i) + ehfeigenup(j)) - &
              ehfeigenup(a) - ehfeigenup(b))
              !AZ 11/6/24
              enddo 
            enddo 
          enddo 
        enddo 
!AZ 10/3/24
       print*,'Enter MP2 Density (closed shell)'
!       a=0
!       b=0
!       pVirt=0 
!       do i=1,neup
!         do j=1,neup
!           do c=1,nb-neup
!           pVirt = tijab(i,:,j,c)*tijab(i,:,j,c) 
!           enddo 
!        enddo 
!       enddo 
!       i=0
!       j=0
!       pOcc=0 
!       do a=1,nb-neup
!         do b=1,nb-neup
!           do k=1,neup
!           pOcc = tijab(:,a,k,b)*tijab(:,a,k,b) 
!           enddo 
!        enddo 
!       enddo

        POcc= 0.0d0
        do i = 1, neup
          do j = 1, neup
            do a = 1, nb-neup
              do b = 1, nb-neup
                do k = 1, neup 
                POcc(i, j) = -0.5d0*POcc(i, j) + tijab(i,a,k,b) * tijab(j,a,k,b)
                end do
             end do
           end do
         end do
       end do
 
        PVirt = 0.0d0
        do a = 1, nb-neup
          do b = 1, nb-neup
            do i = 1, neup
              do j = 1, neup
                do c = 1, nb-neup 
                PVirt(a, b) = 0.5d0*PVirt(a, b) + tijab(i,a,j,c) * tijab(i,b,j,c)
                end do
             end do
           end do
         end do
       end do

!       write(*,*)"Pocc"
!       write(*,*)
!       call print_mat(POcc,neup) 
!       write(*,*)"PVirt"
!       write(*,*)
!       call print_mat(PVirt,nb-neup) 
 
       trPOcc = 0
       do i=1,neup
         trPOcc = trPOcc + pOcc(i,i)  
       enddo
       trPVirt = 0
       do i=1,nb-neup
         trPVirt = trPVirt + pVirt(i,i)  
       enddo


!       write(*,*) 'trPocc meowwwwwwwww', trPOcc 
!       write(*,*) 'trPVirt meowwwwwwwww', trPVirt
 

        tempVec = 0.0
        tempMat = 0.0
        call diagh(POcc,neup,tempVec,tempMat)
        print*,'NOONs of Pocc MP2'
        do i=1,neup
          write(*, '(F12.9)')tempVec(i)
        enddo  
        print*,'sum of NO occs'
        write(*, '(F12.9)')sum(tempVec)
        tempVec = 0.0
        tempMat = 0.0 
        call diagh(PVirt,nb-neup,tempVec,tempMat)
        print*,'NOONs of PVirt MP2'
        do i=1,nb-neup 
          write(*, '(F11.9)')tempVec(i)
        enddo  
        print*,'sum of NO virts'
        write(*, '(F12.9)')sum(tempVec)
!AZ 10/3/24

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
        deallocate(tijab,tempVec,tempMat) 
        !-----------------------------
        ! HERE THE OUTPUT IS GENERATED
        !-----------------------------
!        WRITE(*,'(A27,F30.20,A3)')'      Correlation Energy =',EP2,' au'
!        WRITE(*,'(A27,F30.20,A3)')'     Hartree-Fock Energy =',E0,' au'
!        print*,' '
!        print*,' '
       
        END SUBROUTINE EP2

        subroutine print_mat(mat,n) 
        implicit none
        integer, intent(in) :: n            ! Dimension of the matrix (n x n)
        double precision, intent(in) :: mat(n, n)       ! Input matrix of size n x n
        integer :: i, j                     ! Loop counters
  
    ! Print column headers
        do i = 1, n 
          write(*, '(10F8.3)') (mat(i, j), j = 1, n)
        end do
      
        end subroutine print_mat
