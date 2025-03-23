SUBROUTINE gcEP2(MULTIPLICITY,Cup,Cdown,Ints,NB,Ne,EHFeigenup,EHFeigendown,E0,nuce,SPINCONSERVE)
      ! 
      !

!
      USE omp_lib
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
                                       moIntsBB(:,:,:,:),&
                                       moInts(:,:,:,:),&
                                       tei(:,:,:,:),eps(:),&
                                       tijab(:,:,:,:),&!11/6/24
                                       pocc(:,:),&!11/6/24
                                       poccA(:,:),&!11/6/24
                                       poccB(:,:),&!11/6/24
                                       pvirt(:,:),&!11/6/24
                                       pvirtA(:,:),&!11/6/24
                                       pvirtB(:,:),&!11/6/24
                                       tempVec(:),&!11/7/24
                                       tempMat(:,:),&!11/7/24
                                       onoon(:),&!3/19/25 
                                       vnoon(:),&!3/19/25
                                       Nvec(:),&!3/23/25
                                       gcON(:)!3/23/25

      DOUBLE PRECISION, ALLOCATABLE :: tempInts1(:,:,:,:), & 
                                       tempInts2(:,:,:,:), &
                                       tempInts3(:,:,:,:)
      DOUBLE PRECISION :: TEMP1,DEIJKL
      LOGICAL :: EXCITE, conver,osd2

      DOUBLE PRECISION :: Sz,twoSP1,Nalpha,Nbeta,SEOld1,SEOld2,&
                          SEOld1AA,SEOld1AB,SEOld2AA,SEOld2AB,&
                          EMP2,EMP2AA,EMP2AB,EMP2BA,EMP2BB,EPole,EPoleOld,E,PS,&
                          X1,X2,Np,EPoleOldgc,EPolegc,Egc,PSgc,SEOldgc,SEgc 
!AZ 3/23/25
      DOUBLE PRECISION :: numer,denom,R2,R2AA,R2BB,R2AB,R2BA
!AZ 3/23/25 
!AZ 10/3/24
     DOUBLE PRECISION :: trPOcc,trPVirt                                              
     DOUBLE PRECISION :: wInts, wIntsAA,wIntsAB,wIntsBB
!AZ 10/3/24
      INTEGER :: a,b,c,d,r,s,t,pole,mu,nu,lam,sig,neup,nedown,iter
     
          
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

      allocate(MOInts(nb,nb,nb,nb))
      allocate(tei(nb*2,nb*2,nb*2,nb*2))

!AZ 11/6/24 
      allocate(tijab(nb*2,nb*2,nb*2,nb*2))
      allocate(pOcc(neup*2,neup*2))
      allocate(pOccA(neup,neup))
      allocate(pOccB(neup,neup))
      allocate(pVirt(nb*2-neup*2,nb*2-neup*2))
      allocate(pVirtA(nb-neup,nb-neup))
      allocate(pVirtB(nb-neup,nb-neup))
      allocate(tempVec(nb*2))
      allocate(tempMat(nb*2,nb*2))
!AZ 11/6/24 

!AZ 3/19/25
      allocate(onoon(neup))
      allocate(vnoon(nb-neup))
!AZ 3/19/25
!AZ 3/23/25
      allocate(Nvec(nb*2)) 
      allocate(gcON(nb))
!AZ 3/23/25

      !Quarter transformations: O(N^5)
        print*,' '
        print*,' '
        print*,'==========================================================='
        print*,'           Enter Quarter transformations: O(N^5)           '
        print*,'==========================================================='
        print*,' '
        print*,' '
        print*,' '

!AZ
      !$OMP PARALLEL SHARED(nb,MOInts,Cup,tempInts1,tempInts2,tempInts3) &
      !$OMP & PRIVATE(i,j,k,l,mu,nu,lam,sig)
!
      Write(*,*) 'Hello'
      Write(*,*) omp_get_num_threads()
!

!AZ

!MO ints
!$OMP DO
do i=1, nb
  do mu=1, nb
    tempInts1(i,:,:,:) = tempInts1(i,:,:,:) + &
    (Cup(mu,i)*Ints(mu,:,:,:))
  enddo
enddo
!$OMP END DO

! Add a barrier here to synchronize threads before continuing
!$OMP BARRIER

!$OMP DO
do i=1, nb
  do j=1, nb
    do nu=1, nb
      tempInts2(i,j,:,:) = tempInts2(i,j,:,:) + &
      (Cup(nu,j)*tempInts1(i,nu,:,:))
    enddo
  enddo
enddo
!$OMP END DO

! Add another barrier here
!$OMP BARRIER

!$OMP DO
do i=1, nb
  do j=1, nb
    do k=1, nb
      do lam=1, nb
        tempInts3(i,j,k,:) = tempInts3(i,j,k,:) + &
        (Cup(lam,k)*tempInts2(i,j,lam,:))
      enddo
      do l=1, nb
        do sig=1, nb
          MOInts(i,j,k,l) = MOInts(i,j,k,l) + &
          (Cup(sig,l)*tempInts3(i,j,k,sig))
        enddo
      enddo!end i            
    enddo!end j     
  enddo!end k         
enddo!end l   

!$OMP END DO
!$OMP END PARALLEL

!spin-basis double bar integrals <12||12>

!11/17 this seems to work

!$OMP PARALLEL SHARED(nb,MOInts,tei) &
!$OMP & PRIVATE(p,q,r,s)

!$OMP DO
do s=1, nB*2
  do r=1, nB*2
    do q=1, nB*2
      do p=1, nB*2
        !(pq|rs)
        tei(p,q,r,s) = & !int-->floor division
        ( MOInts(int((p+1)/2),int((r+1)/2),&
          int((q+1)/2),int((s+1)/2))*&
          (kronecker(mod((p),2),mod((r),2)))*& !p r
          (kronecker(mod((q),2),mod((s),2))) ) -& ! q s 
        !(pq|sr)
        ( MOInts(int((p+1)/2),int((s+1)/2),&
          int((q+1)/2),int((r+1)/2))*&
          (kronecker(mod((p),2),mod((s),2)))*& ! p s
          (kronecker(mod((q),2),mod((r),2))) )! q r 

        !Test Print
        !print*,'pqrs',p,q,r,s
        !print*,' (kronecker(mod((p),2),mod((r),2)))',&
        !         (kronecker(mod((p),2),mod((r),2)))
        !print*,' (kronecker(mod((q),2),mod((s),2)))',&
        !         (kronecker(mod((q),2),mod((s),2)))

      enddo
    enddo
  enddo
enddo
!$OMP END DO

! Add a barrier to synchronize threads before ending the parallel region
!$OMP BARRIER

! End the parallel region
!$OMP END PARALLEL

     print*,'   Tran.Done.   '        

!tile spin for eigenvalues
     allocate(eps(nb*2))
     do i=1,(nB*2)
       eps(i) = EHFeigenup(int((i+1)/2))
       print*,'eps',i,eps(i)
     enddo 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
       emp2=0.0d0
       emp2aa=0.0d0
       emp2ab=0.0d0
       emp2ba=0.0d0
       emp2bb=0.0d0 
        do i=1,neup*2
          do j=1,neup*2
            do a=(neup*2)+1,nB*2
              do b=(neup*2)+1,nB*2
              EMP2AA = EMP2AA +  &
                (tei(i,j,a,b))**2 / &  
                ((eps(i) + eps(j)) - &
                eps(a) - eps(b))
              !AZ 11/6/24
                tijab(i,j,a,b) =  &
                (tei(i,j,a,b)) / &  
                ((eps(i) + eps(j)) - &
                eps(a) - eps(b))
              !AZ 11/6/24
              enddo 
            enddo 
          enddo 
        enddo 

        print*,' '
        EMP2AA =  EMP2AA/4.0d0
!        EMP2AB =  EMP2AB
!        EMP2BA =  EMP2BA
!        EMP2BB =  EMP2BB/4.0d0
        print*,'E2AA (Ha) =',EMP2AA
!        print*,'E2AB (Ha) =',EMP2AB
!        print*,'E2BB (Ha) =',EMP2BB
        
        EMP2=EMP2AA!+EMP2AB+EMP2BB
!AZ 11/7
        !amplitudes are tiled by spin
        POcc= 0.0d0
        do i = 1, neup*2
          do j = 1, neup*2
            do a = (neup*2)+1,nB*2
              do b = (neup*2)+1,nB*2
                do k = 1, neup*2
!                POcc(i, j) = POcc(i, j) -0.5d0 * tijab(i,k,a,b) * tijab(k,j,b,a)
                POcc(i, j) = POcc(i, j) -0.5d0 * tijab(k,i,b,a) * tijab(j,k,a,b)
                end do
             end do
           end do
         end do
       end do
 
        PVirt = 0.0d0
        do a = (neup*2)+1,nB*2
          do b = (neup*2)+1,nB*2
            do i = 1, neup*2
              do j = 1, neup*2
                do c = (neup*2)+1,nB*2
                 !pq are virt
                 p = mod(a - 1, neup*2 - nb*2) + 1 
                 q = mod(b - 1, neup*2 - nb*2) + 1
!                PVirt(a, b) = 0.5d0*PVirt(a, b) + tijab(i,j,a,c) * tijab(i,j,b,c)
                 PVirt(p,q) = PVirt(p,q) + &
                                    0.5d0 * tijab(j,i,c,a) * tijab(i,j,b,c)
                end do
             end do
           end do
         end do
       end do

!!      write(*,*)"Pocc"
      write(*,*)
!!      call print_mat(POcc,neup*2) 
!!      write(*,*)"PVirt"
!!      write(*,*)
!!      call print_mat(PVirt,nb*2-neup*2) 
 
!!       trPOcc = 0
       print*,'alpha diag of POcc'
       do i=1,neup*2, 2
!!         trPOcc = trPOcc + pOcc(i,i)
         print*,pOcc(i,i)   
       enddo
!!       trPVirt = 0
!!       do i=1,nb*2-neup*2
!!         trPVirt = trPVirt + pVirt(i,i)  
!!       enddo


       write(*,*) 'trPocc meowwwwwwwww', trPOcc 
!       write(*,*) 'trPVirt meowwwwwwwww', trPVirt

       print*,'take alpha tiles of POcc'
       do i=1,(neup*2),2
         do j=1,(neup*2),2
               pOccA((i+1)/2,(j+1)/2) = pOcc(i,j)
         enddo 
       enddo 
       print*,'alpha diag of POccA'
       do i=1,neup
!!         trPOcc = trPOcc + pOcc(i,i)
         print*,pOccA(i,i)
       enddo
       do i=2,neup*2, 2
         do j=2,neup*2, 2
           do k=1,neup
             do l=1,neup
               pOccB(k,l) = pOcc(i,j)
             enddo 
           enddo 
         enddo 
       enddo 

       print*,'take alpha tiles of PVirt'
       do i=1,(nb*2)-(neup*2),2
         do j=1,(nb*2)-(neup*2),2
               pVirtA((i+1)/2,(j+1)/2) = PVirt(i,j)
         enddo 
       enddo 

       print*,'alpha diag of PVirtA'
       do i=1,nB-neup
         print*,pVirtA(i,i)
       enddo
!       do i=(neup*2)+2,nb*2,2
!         do j=(neup*2)+2,nb*2,2
!           do k=1,nb-neup
!             do l=1,nb-neup
!               pVirtB(k,l) = pVirt(i,j)
!             enddo
!           enddo
!         enddo
!       enddo


        tempVec = 0.0
        tempMat = 0.0
        call diagh(POccA,neup,tempVec,tempMat)
        print*,'NOONs of Pocc MP2'
        onoon(1:neup) = tempVec(1:neup) + 1.0d0 ! add HF occ
        do i=1,neup
          write(*, '(F12.9)')onoon(i) 
        enddo  
        print*,'sum of NO occs'
        write(*, '(F12.9)')sum(onoon) 
        tempVec = 0.0
        tempMat = 0.0 
        call diagh(PVirtA,nb-neup,tempVec,tempMat)
        vnoon(1:nb-neup) = tempVec(1:nb-neup) + 0.0d0 ! add HF virt
        print*,'NOONs of PVirt MP2'
        do i=1,nb-neup
          write(*, '(F11.9)')vnoon(i)
        enddo  
        print*,'sum of NO virts'
        write(*, '(F12.9)')sum(vnoon) 
     
        write(*,*)
        print*,'sum of alpha NOONs'
        write(*, '(F12.9)')sum(vnoon)+sum(onoon) 

!AZ 11/7
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!AZ 3/23/25 tom orb
     !sort occs
     tempVec = 0.0d0 
     tempVec(1:neup) = sort_descending(onoon) 
     gcON(1:neup) = 1.0d0!tempVec(1:neup) !1.0d0
     tempVec = 0.0d0 
     tempVec(1:nb-neup) = sort_descending(vnoon) 
     gcON(neup+1:nB) = 0.0d0!tempVec(1:nb-neup) !0.0d0
     
!tile spin for eigenvalues
     do i=1,(nB*2)
       Nvec(i) = gcON(int((i+1)/2))
       print*,'Nvec before',i,Nvec(i)
     enddo

      Nvec(8) = 0.5d0  !AZ alpha (odd?) is 1/2 occ
      
!     do i=1,(neup*2)
!       eps(i)  = eps(i)*gcON(int((i+1)/2)) !try shifting eigs?
!       print*,'scaled eps(i)',i,eps(i) 
!     enddo 
!AZ 3/23/25

     deallocate(tempInts1,tempInts2,tempInts3)

        !AZ if tiled by spin, 2*poleIndex
        do pole=1,(neup*2)+1,2 !!!!! begin pole search


        iter=0 !max iter 15 for now

        print*,' '
        print*,'orb',int(pole/2)+1 !pole
        print*,' '
        
        E=0.0d0
        Egc=0.0d0 !AZ 3/23 for gc d2 
        SEOld1AA = 0.0d0
        SEOld1AB = 0.0d0
        SEOld2AA = 0.0d0
        SEOld2AB = 0.0d0

        SEOldgc =0.0d0 !AZ 3/23/25
        SEgc = 0.0d0 !AZ 3/23/25 

        SEold1 = 0.0d0
        SEold2 = 0.0d0
        PS=0.0d0
        E = eps(pole)*0.92 !EHFeigenup(pole)*0.92
        Egc = eps(pole)*0.92 !AZ 3/23 for gc d2
        conver = .false. 

        EPoleOld=E
        EPoleOldgc=Egc !AZ 3/23 for gc d2
        do while(conver.eqv..false.)
       !! print*,'conver',conver 



!SE1 2ph   
!AA
        do i=1,Neup*2
          do a=(NeUp*2)+1,NB*2
            do b=(NeUp*2)+1,NB*2
!antisymm?
               SEold1AA = SEold1AA +( ((tei(pole,i,a,b))**2.0d0) / &
               (EPoleOld + eps(i) -eps(a)-eps(b)) )
            enddo
          enddo 
        enddo


!SE2 2hp 
!AA
        do a=(Neup*2)+1,Nb*2
          do i=1,NeUp*2 
            do j=1,NeUp*2
               SEold2AA = SEold2AA + ( ((tei(pole,a,i,j))**2.0d0) / &
               (EPoleOld + eps(a) -eps(i)-eps(j)) )
            enddo
          enddo 
        enddo

!AZ 3/23/25     
!gcD2
!http://dx.doi.org/10.1016/bs.arcc.2017.06.002
!Nqst = nq(1-ns-nt)+ns*nt
!(1-nq)*ns*nt + nq*(1-ns)*(1-nt)
      do q=1,Nb*2   
        do s=1,Nb*2 
          do t=1,Nb*2 
          SEOldgc = SEOldgc + &
!                    (Nvec(q)*(1-Nvec(s)-Nvec(t)))+(Nvec(s)*Nvec(t)) * &
                    (( Nvec(q)*(1-Nvec(s))*(1-Nvec(t)) )+&
                    ( (1-Nvec(q))*(Nvec(s)*(Nvec(t))) ) )* &  
!                    ( ((tei(pole,q,s,t))**2.0d0) / &
                    ( ((tei(pole,q,s,t))*(tei(s,t,pole,q))) / &
                    (EPoleOldgc + eps(q) -eps(s)-eps(t)) )
!          print*,"Nqst", ( Nvec(q)*(1-Nvec(s))*(1-Nvec(t)) )+&
!                         ( (1-Nvec(q))*(Nvec(s)*(Nvec(t))) )
          enddo 
        enddo 
      enddo 
!AZ 3/23/25     

!new NR instead
       SEOld1 = SEOld1AA/2.0d0 
       SEOld2 = SEOld2AA/2.0d0
       EPole = eps(pole) + SEOld1+SEOld2
!!       print*,'sigma(2)',SEOld1+SEOld2
       print*,'SEOld1',SEOld1
       print*,'SEOld2',SEOld2

!AZ 3/23/25     
       SEgc = SEOldgc/2.0d0
       EPolegc = eps(pole) + SEgc
!AZ 3/23/25     

!derivatives
        SEOld1AA = 0.0d0
        SEOld1AB = 0.0d0
        SEOld2AA = 0.0d0
        SEOld2AB = 0.0d0

!AZ 3/23/25    
       SEOldgc = 0.0d0 !zero-out, re-use var 
       SEgc = 0.0d0 !zero-out, re-use var  
!AZ 3/23/25    

       SEOld1=0.0d0
       SEOld2=0.0d0
!SE1 2ph   
!AA
        do i=1,Neup*2
          do a=(NeUp*2)+1,NB*2
            do b=(NeUp*2)+1,NB*2
!antisymm?
               SEold1AA = SEold1AA +( ((tei(pole,i,a,b))**2.0d0) / &
               (EPoleOld + eps(i) -eps(a)-eps(b))**2 )
            enddo
          enddo 
        enddo

!SE2 2hp 
!AA
        do a=(Neup*2)+1,Nb*2
          do i=1,NeUp*2
            do j=1,NeUp*2
               SEold2AA = SEold2AA + ( ((tei(pole,a,i,j))**2.0d0) / &
               (EPoleOld+ eps(a) -eps(i)-eps(j))**2 )
            enddo
          enddo 
        enddo

       SEold1 = -1.0d0*(SEold1AA/2.0d0)
       SEold2 = -1.0d0*(SEold2AA/2.0d0)

!AZ 3/23/25     
!gcD2 deriv
!Nqst = nq(1-ns-nt)+ns*nt
      do q=1,Nb*2
        do s=1,Nb*2
          do t=1,Nb*2
          SEOldgc = SEOldgc + &
                    !(Nvec(q)*(1-Nvec(s)-Nvec(t)))+(Nvec(s)*Nvec(t)) * &
                    (( Nvec(q)*(1-Nvec(s))*(1-Nvec(t)) )+&
                    ( (1-Nvec(q))*(Nvec(s)*(Nvec(t))) ) )* &
                    ( ((tei(pole,q,s,t))**2.0d0) / &
                    (EPoleOldgc + eps(q) -eps(s)-eps(t))**2 )
          enddo
        enddo
      enddo
     SEgc = -1.0d0*(SEoldgc/2.0d0)
!AZ 3/23/25 
    
!ORX AA
        numer= 0.0d0
        denom= 0.0d0
        R2   = 0.0d0
        R2AA = 0.0d0
        R2BB = 0.0d0
        R2AB = 0.0d0
        R2BA = 0.0d0
        

!<pa||pi>/e(a)-e(i) IP 
!<pi||pa>/e(i)-e(a) EA
        if(pole.le.nEup*2) then !ip relax  

!AA aa||aa  pp||ai 
        do i=1,Neup*2
          do a=(NeUp*2)+1,NB*2
            denom = eps(a) - eps(i)
            numer = (tei(pole,a,pole,i))
            numer = numer*numer !*(1-Np)  
            R2AA  =  R2AA + numer/denom 
          enddo
        enddo

       R2BB = R2AA

!AB aa||bb pp||ai  
        do i=1,NeDown*2
          do a=(NeDown*2)+1,NB *2
            denom = eps(a) - eps(i)
            numer = (tei(pole,a,pole,i)) 
            numer = numer*numer !*(1-Np) 
!            R2AB  =  R2AB + numer/denom 
           enddo
         enddo 

        elseif(pole.gt.neUp) then !ea relax

!AA aa||aa  pp||ai 
        do i=1,Neup*2
          do a=(NeUp*2)+1,NB*2
            denom = eps(i) - eps(a)
            numer = (tei(pole,i,pole,a)) 
            numer = numer*numer !*(Np) 
            R2AA  =  R2AA + numer/denom 
          enddo
        enddo

       R2BB = R2AA

!AB aa||bb pp||ai  
        do i=1,NeDown*2
          do a=(NeDown*2)+1,NB*2
            denom = eps(i) - eps(a)
            numer = (tei(pole,i,pole,a))     
            numer = numer*numer !*(Np) 
!            R2AB  =  R2AB + numer/denom 
           enddo
         enddo 


        endif!end relax 

       R2BA = R2AB

       R2 = R2AA + R2AB 

!


       E = (EpoleOld - ((EpoleOld-Epole)/(1-(SEold1+SEold2))))
       Egc = (EpoleOldgc - ((EpoleOldgc-Epolegc)/(1-(SEgc)))) !AZ gc d2 
       PSgc = 1/(1-(SEgc)) !AZ gc d2
       PS = 1/(1-(SEold1+SEold2))
       !print*,'E after NRstep',E      

       iter=iter+1

       if(abs(EPole-EPoleOld).lt.0.0001.or.iter.eq.15) then
       print*,'Koopmans =',eps(pole)
       print*,'D2 (Ha) =',E
       print*,'D2 (eV) =',E*27.2114
       print*,'R2 (Ha) =',R2
       print*,'R2 AA (Ha) =',R2AA
       print*,'R2 AB (Ha) =',R2AB
       print*,'eps+R2 (Ha) =',(0.0d0-(eps(pole)+R2))
       print*,'eps+R2 (eV) =',(0.0d0-(eps(pole)+R2))*27.2114
!       print*,'eps+R2 (Ha) =',((0.0d0-eps(pole))-R2)
!       print*,'eps+R2 (eV) =',((0.0d0-eps(pole))-R2)*27.2114
       print*,'PS =',PS
       print*,'gcD2 (Ha) =',Egc
       print*,'gcD2 (eV) =',Egc*27.2114
       print*,'gsPS =',PSgc
       conver=.true.
       if(iter.eq.15) then 
         print*,'pole not converged after',iter,'iter'
       endif 
       endif 

!AZ 3/23/25
       if(abs(EPolegc-EPoleOldgc).lt.0.0001.or.iter.eq.15) then
       print*,'gcD2 (Ha) =',Egc
       print*,'gcD2 (eV) =',Egc*27.2114
       print*,'gsPS =',PSgc
       conver=.true.
       if(iter.eq.15) then 
         print*,'pole not converged after',iter,'iter'
       endif 
       endif 
       EPoleOldgc = Egc 
       SEgc = 0.0d0 
       SEOldgc = 0.0d0 
!AZ 3/23/25 


       EPoleOld=E
       SEOld1=0.0d0
       SEOld2=0.0d0
        SEOld1AA = 0.0d0
        SEOld1AB = 0.0d0
        SEOld2AA = 0.0d0
        SEOld2AB = 0.0d0


       enddo!while

       enddo !!!!! end pole search

!AZ 10/3/24

       print*,' '
       print*,' '
       print*,'==========================================================='
       print*,'            Results from the MP2 calculation:              '
       print*,'==========================================================='
       print*,' '


        print*,'E2 (Ha) =',EMP2

        print*,' '
 
        WRITE(*,'(A27,F30.20,A3)')'      Correlation Energy =',EMP2,' au'
        WRITE(*,'(A27,F30.20,A3)')'     Hartree-Fock Energy =',E0,' au'
        IF ( nuce .NE. 0) WRITE(*,'(A27,F30.20,A3)')' Nuclear repulsion Energy =',nuce,' au'
        IF ( nuce .NE. 0) WRITE(*,'(A27,F30.20,A3)')'            Total Energy =',E0+EMP2+nuce,' au'
        IF ( nuce .EQ. 0) WRITE(*,'(A27,F30.20,A3)')'            Total Energy =',E0+EMP2,' au'
        print*,' '
        print*,' '

              
!AZ 10/3/24

        deallocate(MOInts,tei)
        deallocate(tijab,tempVec,tempMat)   
contains 

       function kronecker(i,j) result(dij)

      implicit none
      integer,intent(in) :: i,j
      real :: dij

      if(i.eq.j) then
        dij = 1.0
      else
        dij = 0.0
      endif

      end function kronecker

    ! Function to sort a vector in descending order
    !GPT+Stack Overflow
    function sort_descending(vec) result(sorted_vec)
        double precision, dimension(:), intent(in) :: vec
        double precision, dimension(:), allocatable :: sorted_vec
        integer :: i, j
        double precision :: temp
        integer :: N

        N = size(vec)  ! Get the size of the vector
        allocate(sorted_vec(N))  ! Allocate the result array

        ! Copy the original vector to the sorted vector
        sorted_vec = vec

        ! Bubble sort to sort in descending order
        do i = 1, N-1
            do j = 1, N-i
                if (sorted_vec(j) < sorted_vec(j+1)) then
                    ! Swap elements if they are in the wrong order
                    temp = sorted_vec(j)
                    sorted_vec(j) = sorted_vec(j+1)
                    sorted_vec(j+1) = temp
                endif
            enddo
        enddo
    end function sort_descending 
        END SUBROUTINE gcEP2
    
