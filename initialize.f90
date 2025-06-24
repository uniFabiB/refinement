SUBROUTINE initialize
   use, intrinsic :: iso_c_binding     ! Newly Added July 14, 2017
   USE global_variables
   USE mpi
   IMPLICIT NONE
   INCLUDE "fftw3-mpi.f03"             ! Needed, as there is a fftw_mpi_local_size_3d command below

   INTEGER :: i,j,k
   real(pr) :: mode
         

   C_local_alloc = fftw_mpi_local_size_3d(C_n(3), C_n(2), C_n(1)/2+1, MPI_COMM_WORLD, C_local_N, C_local_k_offset)   ! Newly Added July 14, 2017
   local_N = int( C_local_N )
   local_k_offset = int( C_local_k_offset )
   total_local_size = n(1)*n(2)*local_N


   if (rank == 0) THEN
      allocate(global_u(1:n(1), 1:n(2), 1:n(3)))
   end if


   PI = 4.0_pr*ATAN2(1.0_pr,1.0_pr)
   n_dim = 3*n(1)*n(2)*local_N
   dV = 1.0_pr/PRODUCT(REAL(n,pr))

   !  Kcut = 4.0_pr*PI*REAL(n(1),pr)/5.0_pr

   Kcut = 2.0_pr*PI*10_pr*REAL(n(1),pr)/21.0_pr

   Kmax = Kcut/2.0_pr

   !--Set up wavenumbers
   DO i = 0, n(1)-1
      IF (i<=n(1)/2) THEN
         K1(i+1) = 2.0_pr*PI*REAL(i,pr)
      ELSE
         K1(i+1) = 2.0_pr*PI*REAL(i-n(1),pr)
      END IF
   END DO

   DO i = 0,n(2)-1
      IF (i <= n(2)/2) THEN
      K2(i+1) = 2.0_pr*PI*REAL(i,pr)
      ELSE
      K2(i+1) = 2.0_pr*PI*REAL(i-n(2),pr)
   END IF
   END DO

   DO i = 0, n(3)-1
      IF (i<=n(3)/2) THEN
         K3(i+1) = 2.0_pr*PI*REAL(i,pr)
      ELSE
         K3(i+1) = 2.0_pr*PI*REAL(i-n(3),pr)
      END IF
   END DO

   kkmax = ceiling(sqrt(real(n(1),pr)**2/4_pr+real(n(2),pr)**2/4_pr+real(n(3),pr)**2/4_pr)) 


   IF (n(1)<64) THEN
      parallel_data = .FALSE.
   ELSE
      parallel_data = .TRUE.
   END IF
 
END SUBROUTINE
