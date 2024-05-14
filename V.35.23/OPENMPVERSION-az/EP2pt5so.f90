SUBROUTINE EP2pt5so(MULTIPLICITY,Cup,Cdown,Ints,NB,Ne,EHFeigenup,EHFeigendown,E0,nuce,SPINCONSERVE)
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
      DOUBLE PRECISION :: D2, D3,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,&
                          B1,B2,B3,B4,B5,B6,Aterms,Bterms,&
                          secondOrder,thirdOrder,&
                          secondOrderDeriv,thirdOrderDeriv,deriv
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
!!      print*,'   AA   '        
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
!!       print*,'eps',i,eps(i)
     enddo 
  
     deallocate(tempInts1,tempInts2,tempInts3)

        !AZ if tiled by spin, 2*poleIndex
        do pole=1,(neup*2)+1,2 !!!!! begin pole search

        iter=0 !max iter 15 for now

        print*,' '
        print*,'orb',int(pole/2)+1 !pole
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
       !!print*,'SEOld1',SEOld1
       !!print*,'SEOld2',SEOld2

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
       D2 = E
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

!
       print*,'ENTER 3rd Order'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! BEGIN 3rd Order

        iter=0 !reset iter for 3rd order 

        PS=0.0d0 !reset PS for 3rd order 

        !reset 2nd and 3rd order collective terms
        secondOrder=0.0d0
        thirdOrder=0.0d0
        secondOrderDeriv=0.0d0
        thirdOrderDeriv=0.0d0

        E = ((eps(pole)*0.92)+D2)/2.0d0 !average pole: shifted HF plus D2
        conver = .false. 

        EPoleOld=E
        do while(conver.eqv..false.)
       !! print*,'conver',conver 

!Second Order
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
!save 2nd order terms
       secondOrder =  SEOld1+SEOld2
       EPole = eps(pole) + SEOld1+SEOld2
!!       print*,'sigma(2)',SEOld1+SEOld2
!!       print*,'SEOld1',SEOld1
!!       print*,'SEOld2',SEOld2

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

       !save 2nd order derivatives
       secondOrderDeriv = SEold1+SEold2

       !B terms: energy independent
       B1=0.0d0 
       B2=0.0d0
       B3=0.0d0
       B4=0.0d0
       B5=0.0d0
       B6=0.0d0 
       !Closed-shell, diagonal approx p=q
!B1     
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do a=(NeUp*2)+1,Nb*2
              do b=(NeUp*2)+1,Nb*2
                do c=(NeUp*2)+1,Nb*2
                   B1 = B1 + ( & 
                        (tei(pole,b,pole,i)*tei(i,j,a,c)*tei(a,c,b,j))/& 
                        ((eps(i)-eps(b))*(eps(i)+ eps(j) -eps(a)-eps(c)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        B1=(0.5d0)*B1
!B2
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do k=1,NeUp*2
              do a=(NeUp*2)+1,Nb*2
                do b=(NeUp*2)+1,Nb*2
                   B2 = B2 + ( & 
                        (tei(pole,a,pole,j)*tei(i,k,a,b)*tei(j,b,i,k))/& 
                        ((eps(j)-eps(a))*(eps(i)+ eps(k) -eps(a)-eps(b)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        B2=(-0.5d0)*B2
!B3
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do a=(NeUp*2)+1,Nb*2
              do b=(NeUp*2)+1,Nb*2
                do c=(NeUp*2)+1,Nb*2
                   B3 = B3 + ( & 
                        (tei(pole,a,pole,b)*tei(i,j,a,c)*tei(b,c,i,j))/& 
                        ((eps(i)+eps(j)-eps(a)-eps(c))*&
                        (eps(i)+eps(j)-eps(b)-eps(c)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        B3=(0.5d0)*B3
!B4
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do k=1,NeUp*2
              do a=(NeUp*2)+1,Nb*2
                do b=(NeUp*2)+1,Nb*2
                   B4 = B4 + ( & 
                        (tei(pole,j,pole,i)*tei(i,k,a,b)*tei(a,b,j,k))/& 
                        ((eps(j)+eps(k)-eps(a)-eps(b))*&
                        (eps(i)+eps(k)-eps(a)-eps(b)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        B4=(-0.5d0)*B4

!B5
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do a=(NeUp*2)+1,Nb*2
              do b=(NeUp*2)+1,Nb*2
                do c=(NeUp*2)+1,Nb*2
                   B5 = B5 + ( & 
                        (tei(pole,i,pole,a)*tei(b,c,i,j)*tei(a,j,b,c))/& 
                        ((eps(i)-eps(a))*(eps(i)+ eps(j) -eps(b)-eps(c)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        B5=(0.5d0)*B5

!B6 
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do k=1,NeUp*2
              do a=(NeUp*2)+1,Nb*2
                do b=(NeUp*2)+1,Nb*2
                   B6 = B6 + ( & 
                        (tei(pole,i,pole,a)*tei(a,b,j,k)*tei(j,k,i,b))/& 
                        ((eps(i)-eps(a))*(eps(j)+ eps(k) -eps(a)-eps(b)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        B6=(-0.5d0)*B6

       !no B term derivatives

       !Total of B terms
       Bterms=B1+B2+B3+B4+B5+B6 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!

       !A terms: energy dependent
       A1=0.0d0
       A2=0.0d0
       A3=0.0d0
       A4=0.0d0
       A5=0.0d0
       A6=0.0d0
       A7=0.0d0
       A8=0.0d0
       A9=0.0d0
       A10=0.0d0
       A11=0.0d0
       A12=0.0d0
!A1 :: has <VV||VV> term
        do i=1,NeUp*2 !remember spin has tiled indices
          do a=(NeUp*2)+1,Nb*2
            do b=(NeUp*2)+1,Nb*2
              do c=(NeUp*2)+1,Nb*2
                do d=(NeUp*2)+1,Nb*2
                   A1 = A1 + ( & 
                        (tei(pole,i,a,c)*tei(a,c,b,d)*tei(b,d,pole,i))/& 
                        ((EPoleOld+eps(i)-eps(a)-eps(c))*&
                        (EPoleOld+ eps(i) -eps(b)-eps(d)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        A1=(0.25d0)*A1

!A2
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do a=(NeUp*2)+1,Nb*2
              do b=(NeUp*2)+1,Nb*2
                do c=(NeUp*2)+1,Nb*2
                   A2 = A2 + ( & 
                        (tei(pole,i,a,c)*tei(a,j,b,i)*tei(b,c,pole,j))/& 
                        ((EPoleOld+eps(i)-eps(a)-eps(c))*&
                        (EPoleOld+ eps(j) -eps(b)-eps(c)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        A2=(-1.0d0)*A2
!A3
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do a=(NeUp*2)+1,Nb*2
              do b=(NeUp*2)+1,Nb*2
                do c=(NeUp*2)+1,Nb*2
                   A3 = A3 + ( & 
                        (tei(pole,b,i,c)*tei(i,j,a,b)*tei(a,c,pole,j))/& 
                        ((EPoleOld+eps(j)-eps(a)-eps(c))*&
                        (eps(i)+ eps(j) -eps(a)-eps(b)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        A3=(-1.0d0)*A3
!A4
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do k=1,NeUp*2
              do a=(NeUp*2)+1,Nb*2
                do b=(NeUp*2)+1,Nb*2
                   A4 = A4 + ( & 
                        (tei(pole,k,i,j)*tei(i,j,a,b)*tei(a,b,pole,k))/& 
                        ((EPoleOld+eps(k)-eps(a)-eps(b))*&
                        (eps(i)+ eps(j) -eps(a)-eps(b)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        A4=(0.25d0)*A4
!A5
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do a=(NeUp*2)+1,Nb*2
              do b=(NeUp*2)+1,Nb*2
                do c=(NeUp*2)+1,Nb*2
                   A5 = A5 + ( & 
                        (tei(pole,j,a,b)*tei(a,c,i,j)*tei(i,b,pole,c))/& 
                        ((EPoleOld+eps(j)-eps(a)-eps(b))*&
                        (eps(i)+ eps(j) -eps(a)-eps(c)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        A5=(-1.0d0)*A5
!A6
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do k=1,NeUp*2
              do a=(NeUp*2)+1,Nb*2
                do b=(NeUp*2)+1,Nb*2
                   A6 = A6 + ( & 
                        (tei(pole,j,a,b)*tei(a,b,i,k)*tei(i,k,pole,j))/& 
                        ((EPoleOld+eps(j)-eps(a)-eps(b))*&
                        (eps(i)+ eps(k) -eps(a)-eps(b)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        A6=(0.25d0)*A6
!A7
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do a=(NeUp*2)+1,Nb*2
              do b=(NeUp*2)+1,Nb*2
                do c=(NeUp*2)+1,Nb*2
                   A7 = A7 + ( & 
                        (tei(pole,c,i,j)*tei(i,j,a,b)*tei(a,b,pole,c))/& 
                        ((eps(i)+eps(j)-EPoleOld-eps(c))*&
                        (eps(i)+ eps(j) -eps(a)-eps(b)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        A7=(-0.25d0)*A7
!A8
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do k=1,NeUp*2
              do a=(NeUp*2)+1,Nb*2
                do b=(NeUp*2)+1,Nb*2
                   A8 = A8 + ( & 
                        (tei(pole,b,i,k)*tei(i,j,a,b)*tei(a,k,pole,j))/& 
                        ((eps(i)+eps(k)-EPoleOld-eps(b))*&
                        (eps(i)+ eps(j) -eps(a)-eps(b)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        A8=(1.0d0)*A8

!A9
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do a=(NeUp*2)+1,Nb*2
              do b=(NeUp*2)+1,Nb*2
                do c=(NeUp*2)+1,Nb*2
                   A9 = A9 + ( & 
                        (tei(pole,b,a,c)*tei(a,c,i,j)*tei(i,j,pole,b))/& 
                        ((eps(i)+eps(j)-EPoleOld-eps(b))*&
                        (eps(i)+ eps(j) -eps(a)-eps(c)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        A9=(-0.25d0)*A9

!A10
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do k=1,NeUp*2
              do a=(NeUp*2)+1,Nb*2
                do b=(NeUp*2)+1,Nb*2
                   A10 = A10 + ( & 
                        (tei(pole,k,a,j)*tei(a,b,i,k)*tei(i,j,pole,b))/& 
                        ((eps(i)+eps(j)-EPoleOld-eps(b))*&
                        (eps(i)+ eps(k) -eps(a)-eps(b)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        A10=(1.0d0)*A10

!A11
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do k=1,NeUp*2
              do a=(NeUp*2)+1,Nb*2
                do b=(NeUp*2)+1,Nb*2
                   A11 = A11 + ( & 
                        (tei(pole,b,i,k)*tei(i,a,j,b)*tei(j,k,pole,a))/& 
                        ((eps(j)+eps(k)-EPoleOld-eps(a))*&
                        (eps(i)+ eps(k) -EPoleOld-eps(b)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        A11=(1.0d0)*A11
!A12 :: has <OO||OO> term
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do k=1,NeUp*2
              do l=1,NeUp*2
                do a=(NeUp*2)+1,Nb*2
                   A12 = A12 + ( & 
                        (tei(pole,a,i,l)*tei(i,l,j,k)*tei(j,k,pole,a))/& 
                        ((eps(j)+eps(k)-EPoleOld-eps(a))*&
                        (eps(i)+ eps(l) -EPoleOld-eps(a)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        A12=(-0.25d0)*A12


       !Total of A terms 
       Aterms=A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+A11+A12 
      
       !Total of 3rd order terms
       thirdOrder=Aterms+Bterms
       !New pole: epsHF + Sigma(2) + Sigma(3)
       EPole = eps(pole) + secondOrder + (0.5d0*thirdOrder)

       !A term derivatives
       !like 2nd order, use the 3rd order A term vars for deriv vars
       !square denoms with E
       A1=0.0d0
       A2=0.0d0
       A3=0.0d0
       A4=0.0d0
       A5=0.0d0
       A6=0.0d0
       A7=0.0d0
       A8=0.0d0
       A9=0.0d0
       A10=0.0d0
       A11=0.0d0
       A12=0.0d0
!AZ fix deriv terms u'v+uv'
!derivA1 :: derivs with 2 EPoleOld have factor (-1)*(-1)=1
        do i=1,NeUp*2 !remember spin has tiled indices
          do a=(NeUp*2)+1,Nb*2
            do b=(NeUp*2)+1,Nb*2
              do c=(NeUp*2)+1,Nb*2
                do d=(NeUp*2)+1,Nb*2
                   A1 = A1 + (-0.25d0)*( & 
                        (tei(pole,i,a,c)*tei(a,c,b,d)*tei(b,d,pole,i))/& 
                        (((EPoleOld+eps(i)-eps(a)-eps(c)))*&
                        ((EPoleOld+ eps(i) -eps(b)-eps(d))**2))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !A1=(-1.0d0)*A1 
        !A1=(-0.25d0)*A1 !0.25
        do i=1,NeUp*2 !remember spin has tiled indices
          do a=(NeUp*2)+1,Nb*2
            do b=(NeUp*2)+1,Nb*2
              do c=(NeUp*2)+1,Nb*2
                do d=(NeUp*2)+1,Nb*2
                   A1 = A1 +  (-0.25d0)*( & 
                        (tei(pole,i,a,c)*tei(a,c,b,d)*tei(b,d,pole,i))/& 
                        (((EPoleOld+eps(i)-eps(a)-eps(c))**2)*&
                        ((EPoleOld+ eps(i) -eps(b)-eps(d))))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
       ! A1=(-0.25d0)*A1 !0.25
        !deriv factor 
        !A1=(-1.0d0)*A1 !change sign? AZ 2/22
!derivA2 (-1)*(-1)
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do a=(NeUp*2)+1,Nb*2
              do b=(NeUp*2)+1,Nb*2
                do c=(NeUp*2)+1,Nb*2
                   A2 = A2 +  (1.0d0)*( & 
                        (tei(pole,i,a,c)*tei(a,j,b,i)*tei(b,c,pole,j))/& 
                        (((EPoleOld+eps(i)-eps(a)-eps(c)))*&
                        ((EPoleOld+ eps(j) -eps(b)-eps(c))**2))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !A2=(-1.0d0)*A2
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do a=(NeUp*2)+1,Nb*2
              do b=(NeUp*2)+1,Nb*2
                do c=(NeUp*2)+1,Nb*2
                   A2 = A2 + (1.0d0)*( & 
                        (tei(pole,i,a,c)*tei(a,j,b,i)*tei(b,c,pole,j))/& 
                        (((EPoleOld+eps(i)-eps(a)-eps(c))**2)*&
                        ((EPoleOld+ eps(j) -eps(b)-eps(c))))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        !A2=(-1.0d0)*A2 !-1 
        !deriv factor
        !A2=(-1.0d0)*A2 !change sign? AZ 2/22
!derivA3 :: one E in deriv, (-1) factor 
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do a=(NeUp*2)+1,Nb*2
              do b=(NeUp*2)+1,Nb*2
                do c=(NeUp*2)+1,Nb*2
                   A3 = A3 + (1.0d0)*( & 
                        (tei(pole,b,i,c)*tei(i,j,a,b)*tei(a,c,pole,j))/& 
                        (((EPoleOld+eps(j)-eps(a)-eps(c))**2)*&
                        (eps(i)+ eps(j) -eps(a)-eps(b)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        !A3=(-1.0d0)*A3 !-1
        !deriv factor
        !A3=(1.0d0)*A3 !change sign? AZ 2/22
!derivA4 :: (-1)
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do k=1,NeUp*2
              do a=(NeUp*2)+1,Nb*2
                do b=(NeUp*2)+1,Nb*2
                   A4 = A4 + (-0.25d0)*( & 
                        (tei(pole,k,i,j)*tei(i,j,a,b)*tei(a,b,pole,k))/& 
                        (((EPoleOld+eps(k)-eps(a)-eps(b))**2)*&
                        (eps(i)+ eps(j) -eps(a)-eps(b)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        !A4=(0.25d0)*A4 !0.25
        !deriv factor
        !A4=(-1.0d0)*A4 
!derivA5 :: (-1)
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do a=(NeUp*2)+1,Nb*2
              do b=(NeUp*2)+1,Nb*2
                do c=(NeUp*2)+1,Nb*2
                   A5 = A5 + (1.0d0)*( & 
                        (tei(pole,j,a,b)*tei(a,c,i,j)*tei(i,b,pole,c))/& 
                        (((EPoleOld+eps(j)-eps(a)-eps(b))**2)*&
                        (eps(i)+ eps(j) -eps(a)-eps(c)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        !A5=(-1.0d0)*A5 !-1
        !deriv factor
        !A5=(-1.0d0)*A5
!derivA6 :: (-1)
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do k=1,NeUp*2
              do a=(NeUp*2)+1,Nb*2
                do b=(NeUp*2)+1,Nb*2
                   A6 = A6 + (-0.25d0)*( & 
                        (tei(pole,j,a,b)*tei(a,b,i,k)*tei(i,k,pole,j))/& 
                        (((EPoleOld+eps(j)-eps(a)-eps(b))**2)*&
                        (eps(i)+ eps(k) -eps(a)-eps(b)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        !A6=(0.25d0)*A6 !0.25
        !deriv factor 
        !A6=(-1.0d0)*A6
!derivA7 :: (+1) since -EPole 
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do a=(NeUp*2)+1,Nb*2
              do b=(NeUp*2)+1,Nb*2
                do c=(NeUp*2)+1,Nb*2
                   A7 = A7 + (-0.25d0)*( & 
                        (tei(pole,c,i,j)*tei(i,j,a,b)*tei(a,b,pole,c))/& 
                        (((eps(i)+eps(j)-EPoleOld-eps(c))**2)*&
                        (eps(i)+ eps(j) -eps(a)-eps(b)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        !A7=(-0.25d0)*A7 !-0.25
        !deriv factor
        !A7=(1.0d0)*A7 
!derivA8 :: (+1) since -EPole
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do k=1,NeUp*2
              do a=(NeUp*2)+1,Nb*2
                do b=(NeUp*2)+1,Nb*2
                   A8 = A8 + (1.0d0)*( & 
                        (tei(pole,b,i,k)*tei(i,j,a,b)*tei(a,k,pole,j))/& 
                        (((eps(i)+eps(k)-EPoleOld-eps(b))**2)*&
                        (eps(i)+ eps(j) -eps(a)-eps(b)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        !A8=(1.0d0)*A8  !1
        !deriv factor 
        !A8=(1.0d0)*A8
!derivA9 :: (+1) since -Epole
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do a=(NeUp*2)+1,Nb*2
              do b=(NeUp*2)+1,Nb*2
                do c=(NeUp*2)+1,Nb*2
                   A9 = A9 + (-0.25d0)*( & 
                        (tei(pole,b,a,c)*tei(a,c,i,j)*tei(i,j,pole,b))/& 
                        (((eps(i)+eps(j)-EPoleOld-eps(b))**2)*&
                        (eps(i)+ eps(j) -eps(a)-eps(c)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        !A9=(-0.25d0)*A9 !-0.25
        !deriv factor
        !A9=(1.0d0)*A9 
!derivA10 :: (+1) since -EPole
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do k=1,NeUp*2
              do a=(NeUp*2)+1,Nb*2
                do b=(NeUp*2)+1,Nb*2
                   A10 = A10 + (1.0d0)*( & 
                        (tei(pole,k,a,j)*tei(a,b,i,k)*tei(i,j,pole,b))/& 
                        (((eps(i)+eps(j)-EPoleOld-eps(b))**2)*&
                        (eps(i)+ eps(k) -eps(a)-eps(b)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        !A10=(1.0d0)*A10 !1
        !deriv factor
        !A10=(1.0d0)*A10
!derivA11 :: (+1)*(+1)=(+1) since -Epole and -Epole
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do k=1,NeUp*2
              do a=(NeUp*2)+1,Nb*2
                do b=(NeUp*2)+1,Nb*2
                   A11 = A11 + (1.0d0)*( & 
                        (tei(pole,b,i,k)*tei(i,a,j,b)*tei(j,k,pole,a))/& 
                        (((eps(j)+eps(k)-EPoleOld-eps(a))**2)*&
                        (eps(i)+ eps(k) -EPoleOld-eps(b)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !A11=(1.0d0)*A11 
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do k=1,NeUp*2
              do a=(NeUp*2)+1,Nb*2
                do b=(NeUp*2)+1,Nb*2
                   A11 = A11 + (1.0d0)*( & 
                        (tei(pole,b,i,k)*tei(i,a,j,b)*tei(j,k,pole,a))/& 
                        ((eps(j)+eps(k)-EPoleOld-eps(a))*&
                        ((eps(i)+ eps(k) -EPoleOld-eps(b))**2))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        !A11=(1.0d0)*A11 !1
        !deriv factor
        !A11=(1.0d0)*A11
!derivA12 ::  (+1)*(+1)=(+1) since -Epole and -Epole
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do k=1,NeUp*2
              do l=1,NeUp*2
                do a=(NeUp*2)+1,Nb*2
                   A12 = A12 + (-0.25d0)*( & 
                        (tei(pole,a,i,l)*tei(i,l,j,k)*tei(j,k,pole,a))/& 
                        (((eps(j)+eps(k)-EPoleOld-eps(a))**2)*&
                        (eps(i)+ eps(l) -EPoleOld-eps(a)))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !A12=(1.0d0)*A12
        do i=1,NeUp*2 !remember spin has tiled indices
          do j=1,NeUp*2
            do k=1,NeUp*2
              do l=1,NeUp*2
                do a=(NeUp*2)+1,Nb*2
                   A12 = A12 + (-0.25d0)*( & 
                        (tei(pole,a,i,l)*tei(i,l,j,k)*tei(j,k,pole,a))/& 
                        ((eps(j)+eps(k)-EPoleOld-eps(a))*&
                        ((eps(i)+ eps(l) -EPoleOld-eps(a))**2))&
                             )
                enddo 
              enddo 
            enddo
          enddo
        enddo
        !factor 
        !A12=(-0.25d0)*A12 !-0.25
        !deriv factor
        !A12=(1.0d0)*A12 !change sign?
!AZ 2/21 
!now sum derivs, check the terms with 2 pole variables, check ()'s
!add 2nd order derivs as well, put into Newton and PS formula 

       !save 3rd order derivs
       thirdOrderDeriv=A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+A11+A12
       !Compute new pole NR step 
       deriv = secondOrderDeriv+ 0.5d0*(thirdOrderDeriv)
       E = (EpoleOld - ((EpoleOld-Epole)/(1-(deriv))))
       PS = 1/(1-(deriv))
       !print*,'E after NRstep',E      

       iter=iter+1

       if(abs(EPole-EPoleOld).lt.0.00001.or.iter.eq.15) then
       D3 = E
       print*,'Koopmans =',eps(pole)
       print*,'D2.5 (Ha) =',D3
       print*,'D2.5 (eV) =',D3*27.2114
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


!!!!!!!!! END 3rd Order


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

 
        END SUBROUTINE EP2pt5so
    
