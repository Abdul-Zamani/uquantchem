SUBROUTINE mom(COEFF,COEFFold,SAO,NB,NOcc,P,spin)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NB,Nocc,spin
      DOUBLE PRECISION, INTENT(INOUT) :: COEFF(NB,NB) !cnew
      DOUBLE PRECISION, INTENT(IN) :: COEFFold(NB,NB) !cold
      DOUBLE PRECISION, INTENT(IN) :: SAO(NB,NB)
      DOUBLE PRECISION, INTENT(OUT) :: P(NB,NB)
      DOUBLE PRECISION :: s(NB),O(NB,NB),nVirt,tempVec(NB),tempmat(NB,NB),temp
      INTEGER :: i,j,mu,nu,tempidx
      INTEGER :: sorted_indices(NB)




    ! Open the file for reading (replace 24 with your actual unit number)
!     open(24, FILE='alphaMO.dat', ACTION='READ')

    ! Read the matrix back into Cup
!    DO I = 1, NB
!        DO J = 1, NB
            ! Read each value from the file into Cup(J, I)
 !           read(24, *) COEFFOld(I, J) !COEFFOld(J, I)
 !       END DO
 !   END DO

    ! Close the file
!    close(24)

    !COEFFold(:,1) =  COEFFold(:,1)*0.5d0 
    !COEFF(:,1) =  COEFFOld(:,1)*0.5d0 

    do i=1, NB
      write(*,'(A,1x,F10.4,1x,A,1x,F10.4)')&
      'Cold',COEFFold(i,1),'C', COEFF(i,1)
    enddo 

    ! Optionally, print the matrix to verify
    !print *, "Regenerated Matrix C (Cup):"
    !DO I = 1, NB
    !    print *, COEFFold(:,1) !, I)  ! Print each column of the matrix
    !END DO

    s=0.0d0 
    j=0

        ! Outer loop over j (index for columns of the result)
        do j = 1, NB
            ! Loop over i (index for rows of the result)
            do i = 1, NOcc
                ! Temporary variable to accumulate the sum
                temp = 0.0
                ! Loop over mu (rows of COEFFold and SAO)
                do mu = 1, NB
                    ! Loop over nu (columns of SAO and COEFF)
                    do nu = 1, NB
                        ! Accumulate the product in temp_sum
                        temp = temp+ &
                        ((COEFFold(mu,i)*SAO(mu,nu)*COEFF(nu,j))**2.0d0) 
                    enddo
                enddo
                ! Store the result in s(j)
                s(j) = s(j) + temp
            enddo
        enddo       
!^recent
        s=s**0.5d0 

      tempVec=0.0d0
!      tempmat=0.0d0 
!      tempVec = tempmat(:,8)
!      tempmat(:,8) = tempmat(:,7)
!      tempmat(:,7) = tempVec
!      coeff = tempmat 
!      O = matmul((matmul(transpose(COEFFold),SAO)),COEFF)
!      O = matmul(transpose(COEFFold),matmul(SAO,COEFF))

!      s = 0.0d0 
      do j = 1, nB!NOcc
        do i = 1, NB
          !s(j) = (s(j) + ((O(i, j))**2.0d0))**(0.5d0)
          !s(j) = (s(j) + ((O(i, j))))
          !s(j) = (s(j) + ((O(i, j))**2.0d0))

        end do
      end do

      do j=1,NB!NOcc
        write(*,'(A,I5,8X,F20.6)')'sj',j,s(j)  
      enddo

! print*,'c(:,1) after',Coeff(:, 1)
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
      if (tempVec(i) < tempVec(j)) then
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
  do i = 1, NB
!goes here   tempMat(:, i) = COEFF(:,sorted_indices(i))
   if(sorted_indices(i).eq.1.and.spin.eq.1) then
   print*,'tom is idx and sort dix',i,sorted_indices(i)
   tempidx = i 
   endif 
   tempMat(:,i) = COEFF(:,sorted_indices(i))
  end do
    do i=1, NB
      write(*,'(A,1x,F10.4)')&
      'tempC', tempMat(i,tempidx)
    enddo
  COEFF=tempMat

  !find tom 
  print*,'tom sorted index #',tempidx
  coeff(:,tempidx) =  coeff(:,tempidx)/2.0d0 
    do i=1, NB
      write(*,'(A,1x,F10.4)')&
      'newCtom', coeff(i,tempidx)
    enddo
  !COEFF(:,1) =COEFF(:,1)*0.5 !this works, why not inside ^ 

      P = MATMUL(COEFF,TRANSPOSE(COEFF))

!AZ 3/28/35
!This subroutine performs the maximum overlap method
!The occupied orbitals are reordered according to a metric
!Let s = sum_j^nocc Oij*Oij 
!pj = sum_i_nocc Oij = sum_nu[sum_mu(sum_i^nOcc Cold_i,mu)*S_mu,nu]*Cnew_nu,j
!Or as a matmul over occupied orbitals:
!O = (Cold)^T S Cnew 
END SUBROUTINE mom
          
