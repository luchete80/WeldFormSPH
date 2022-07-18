PROGRAM WORKSECTIONS
      INTEGER N, I, NTHREADS, TID, OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
      PARAMETER (N=10000000)
      REAL A(N), B(N), C(N), D(N)

!     Some initializations
      DO I = 1, N
        A(I) = I * 1.5
        B(I) = I + 22.35
        C(N) = 0.0
        D(N) = 0.0 
      ENDDO

!$OMP PARALLEL SHARED(A,B,C,D,NTHREADS), PRIVATE(I,TID)
      TID = OMP_GET_THREAD_NUM()
      IF (TID .EQ. 0) THEN
        NTHREADS = OMP_GET_NUM_THREADS()
        PRINT *, 'Number of threads =', NTHREADS
      END IF
      PRINT *, 'Thread',TID,' starting...'

!$OMP SECTIONS

!$OMP SECTION
      PRINT *, 'Thread',TID,' doing section 1'
      DO I = 1, N
         C(I) = A(I) + B(I)
         if (i.lt.10) then
            WRITE(*,100) TID,I,C(I)
            end if
 100     FORMAT(' Thread',I2,': C(',I2,')=',F8.2)
      ENDDO

!$OMP SECTION
      PRINT *, 'Thread',TID,' doing section 2'
      DO I = 1, N
         if (i.lt.10) then
         D(I) = A(I) * B(I)
         WRITE(*,200) TID,I,D(I)
 200     FORMAT(' Thread',I2,': D(',I2,')=',F8.2)
         endif
      ENDDO


!$OMP END SECTIONS NOWAIT
      PRINT *, 'Thread',TID,' done.'

!$OMP END PARALLEL
END PROGRAM WORKSECTIONS