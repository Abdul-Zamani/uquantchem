SUBROUTINE mom(COEFF,COEFFold,SAO,NB,NOcc,spin)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NB,Nocc,spin
      DOUBLE PRECISION, INTENT(INOUT) :: COEFF(NB,NB) !cnew !occ and virt
      DOUBLE PRECISION, INTENT(IN) :: COEFFold(NB,NB) !cold !occ and virt
      DOUBLE PRECISION, INTENT(IN) :: SAO(NB,NB)
      DOUBLE PRECISION :: s(NB),O(NB,NB),nVirt,tempVec(NB),tempmat(NB,NB),temp
      DOUBLE PRECISION :: tempmat2(NB,NB)
      INTEGER :: i,j,k,mu,nu,tempidx
      INTEGER :: sorted_indices(NB)



    tempmat2 = COEFFold
    if(spin.eq.1) tempmat2(:,1)= tempmat2(:,1)*0.5d0 

    do i=1, NB
!!!      write(*,'(A,1x,F10.4,1x,A,1x,F10.4)')&
!!!      'Cold',tempmat2(i,1),'C', COEFF(i,1)
    enddo 

    s=0.0d0 
    j=0


     !do projection 
     !Cold^T * S * C 
     O = matmul(transpose(tempmat2),matmul(SAO,COEFF))
     do j=1,nB 
       s(j) = 0.0d0 
       do i=1,nocc
       !s(j) = s(j) + (O(i,j))
       s(j) = s(j) + (O(i,j))**2.0d0  
       enddo
     enddo 
     s=s**0.5d0 


     tempVec=0.0d0

      do j=1,NB!NOcc
        write(*,'(A,I5,8X,F20.6)')'sj',j,s(j)  
      enddo

      nVirt = NOcc - sum(s(1:NOcc)) 
!      write(*,*)
      write(*,'(A,F10.4)'),'nVirt',nVirt
!      write(*,*)

      sorted_indices = [(i, i=1, NB)]
      tempVec = s

  ! Sort the array s_temp in descending order using a simple bubble sort (for
  ! illustration)
  do i = 1, NB-1
    do j = i+1, NB
      if (tempVec(i) < tempVec(j)) then !gt ascending lt descending
        ! Swap s_temp(i) and s_temp(j)
        temp = tempVec(i)
        tempVec(i) = tempVec(j)
        tempVec(j) = temp

        ! Swap the corresponding indices in index
        tempidx = sorted_indices(i)
        sorted_indices(i) = sorted_indices(j)
        sorted_indices(j) = tempidx
      end if
    end do
  end do

  ! Output the sorted vector and the indices of the original positions
  do j=1,NB!NOcc
    write(*,'(A,I5,8X,F20.6)')'sorted sj',j,tempVec(j)
  enddo
  do j=1,NB!NOcc
    write(*,*)'sorted idx',j,sorted_indices(j) 
  enddo 
     
  tempMat = 0.0d0
  tempVec = 0.0d0 

  !tempMat = COEFF   ! make it equal, coeff intent in 

  ! Rearrange the columns of C according to the index array
  do i=1, NB
   tempMat(:, i) = COEFF(:,sorted_indices(i))
   if(sorted_indices(i).eq.1.and.spin.eq.1) then
   tempidx = i 
   print*,'tom is idx and sort idx',i,sorted_indices(i)
   tempMat(:,i) = COEFF(:,sorted_indices(i))*0.5d0 
   endif 
  end do
    do i=1, NB
!      write(*,'(A,1x,F10.4)')&
!      'tempC', tempMat(i,tempidx)
    enddo
!!  COEFF=tempMat

  !find tom 
  print*,'tom sorted index #',tempidx
    do i=1, NB
!!!      write(*,'(A,1x,F10.4)')&
!!!      'C(:,tom)', tempmat(i,tempidx)
    enddo


       COEFF = tempMat !spit out reordered C  

!AZ 3/28/35
!This subroutine performs the maximum overlap method
!The occupied orbitals are reordered according to a metric
!Let s = sum_j^nocc Oij*Oij 
!pj = sum_i_nocc Oij = sum_nu[sum_mu(sum_i^nOcc Cold_i,mu)*S_mu,nu]*Cnew_nu,j
!Or as a matmul over occupied orbitals:
!O = (Cold)^T S Cnew 
END SUBROUTINE mom
      subroutine Print_Matrix_Full_Real(IOut,AMat,M,N)

!This subroutine prints a real matrix that is fully dimension - i.e.,
!not stored in packed form. AMat is the matrix, which is dimensioned
!(M,N).

      implicit none
      integer,intent(in)::IOut,M,N
      double precision,dimension(M,N),intent(in)::AMat
      integer,parameter::NColumns=5
      integer::i,j,IFirst,ILast

 1000 Format(1x,A)
 2000 Format(5x,5(7x,I7))
 2010 Format(1x,I7,5F14.6)

      do IFirst = 1,N,NColumns
        ILast = Min(IFirst+NColumns-1,N)
        write(IOut,2000) (i,i=IFirst,ILast)
        do i = 1,M
          write(IOut,2010) i,(AMat(i,j),j=IFirst,ILast)
        enddo
      enddo

      return
      end subroutine Print_Matrix_Full_Real

 
