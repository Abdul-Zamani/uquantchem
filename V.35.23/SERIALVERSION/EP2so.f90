SUBROUTINE EP2so(MULTIPLICITY,Cup,Cdown,Ints,NB,Ne,EHFeigenup,EHFeigendown,E0,nuce,SPINCONSERVE)
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
                                       moIntsBB(:,:,:,:),&
                                       moInts(:,:,:,:),&
                                       tei(:,:,:,:),eps(:)
      DOUBLE PRECISION, ALLOCATABLE :: tempInts1(:,:,:,:), & 
                                       tempInts2(:,:,:,:), &
                                       tempInts3(:,:,:,:)
      DOUBLE PRECISION :: TEMP1,DEIJKL
      LOGICAL :: EXCITE, conver,osd2

      DOUBLE PRECISION :: Sz,twoSP1,Nalpha,Nbeta,SEOld1,SEOld2,&
                          SEOld1AA,SEOld1AB,SEOld2AA,SEOld2AB,&
                          EMP2,EMP2AA,EMP2AB,EMP2BA,EMP2BB,EPole,EPoleOld,E,PS,&
                          X1,X2
      INTEGER :: a,b,c,d,r,s,pole,mu,nu,lam,sig,neup,nedown,iter
     
          
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

      !Quarter transformations: O(N^5)
        print*,' '
        print*,' '
        print*,'==========================================================='
        print*,'           Enter Quarter transformations: O(N^5)           '
        print*,'==========================================================='
        print*,' '
        print*,' '
        print*,' '

!MO ints
      print*,'   AA   '        
      tempInts1=0.0d0 
      tempInts2=0.0d0
      tempInts3=0.0d0
      moInts = 0.0d0
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
                MOInts(i,j,k,l) = MOInts(i,j,k,l) + &
                (Cup(sig,l)*tempInts3(i,j,k,sig))
              enddo      
            enddo!end i            
          enddo!end j     
        enddo!end k         
      enddo!end l   

      !spin integration 
      do i=1,(nb*2) 
        do j=1,(nb*2) 
          do k=1,(nb*2) 
            do l=1,(nb*2) 
              if( (mod(l,2).eq.mod(j,2)).and.(mod(i,2).eq.mod(k,2)) ) then
!              X1 = MOInts(int((i+1)/2),int((k+1)/2),int((j+1)/2),int((l+1)/2))
              else
!              X1 = 0.0d0
              endif
              if( (mod(i,2).eq.mod(l,2)).and.(mod(j,2).eq.mod(k,2)) ) then
!              X2 = MOInts(int((i+1)/2),int((l+1)/2),int((j+1)/2),int((k+1)/2))
              else
!              X2 = 0.0d0
              endif
!              tei(i,k,j,l) = X1-X2 !ijkl - ijlk
            enddo
          enddo 
        enddo
     enddo 

!spin-basis double bar integrals <12||12>

!11/17 this seems to work

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

     print*,'   Tran.Done.   '        

!tile spin for eigenvalues
     allocate(eps(nb*2))
     do i=1,(nB*2)
       eps(i) = EHFeigenup(int((i+1)/2))
       print*,'eps',i,eps(i)
     enddo 
  
     deallocate(tempInts1,tempInts2,tempInts3)
        !AZ if tiled by spin, 2*poleIndex
        do pole=1,neup*2 !!!!! begin pole search

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
        E = eps(pole)*0.92 !EHFeigenup(pole)*0.92
        conver = .false. 

        EPoleOld=E
        do while(conver.eqv..false.)
       !! print*,'conver',conver 
!SE1  
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


!SE2
!AA
        do a=(Neup*2)+1,Nb*2
          do i=1,NeUp*2 
            do j=1,NeUp*2
               SEold2AA = SEold2AA + ( ((tei(pole,a,i,j))**2.0d0) / &
               (EPoleOld + eps(a) -eps(i)-eps(j)) )
            enddo
          enddo 
        enddo
     

!new NR instead
       SEOld1 = SEOld1AA/2.0d0 
       SEOld2 = SEOld2AA/2.0d0
       EPole = eps(pole) + SEOld1+SEOld2
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

!SE2
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

       E = (EpoleOld - ((EpoleOld-Epole)/(1-(SEold1+SEold2))))
       PS = 1/(1-(SEold1+SEold2))
       !print*,'E after NRstep',E      

       iter=iter+1

       if(abs(EPole-EPoleOld).lt.0.0001.or.iter.eq.15) then
       print*,'Koopmans =',eps(pole)
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



        deallocate(MOInts,tei)
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

 
        END SUBROUTINE EP2so
    
