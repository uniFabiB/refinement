!======================================
! MODULE CONTAINING INTERFACES FOR 
! FFTW
!
! (*) ffftwfunction
!
! CONTAINS SUBROUTINE
!
! (*) fftfwd(u, fu): perform r2c fft
! (*) fftbwd(fu, u): perform c2r fft
! (*) fftfwd_m(u, fu, m)
! (*) fftbwd_m(fu, u, m)
! (*) init_fft()
! (*) fft_deallocate()
!======================================



MODULE fftwfunction   ! Newly added on March 20, 2017
   
  use, intrinsic :: iso_c_binding
!      USE global_variables
!      IMPLICIT NONE 
!      INCLUDE "mpif.h"
!      include 'fftw3-mpi.f03'

  type(C_PTR) :: fwdplan , bwdplan, tmpdata_r, tmpdata_cx
  real(c_double), pointer :: tmppointer_r(:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: tmppointer_cx(:,:,:)

  ! note that for mpi-fftw, the size of tmpointer_r is 2*(Nx/2+1) * Ny * Nz
  ! suppose that u is a real input of size Nx * Ny * Nz
  ! tmpointer_r(i,j,k) = tmpointer[(k*M+j)*(N+2) + i] = u(i,j,k)

CONTAINS

!=======================================
! Forward Fourier transform
!=======================================
  SUBROUTINE fftfwd(u,fu)

    USE global_variables
    IMPLICIT NONE 
    INCLUDE "mpif.h"
    include 'fftw3-mpi.f03'


    real(pr), DIMENSION (1:n(1),1:n(2),1:local_N), INTENT(IN) :: u
    COMPLEX(pr), DIMENSION (1:n(1)/2+1,1:n(2),1:local_N), INTENT(OUT) :: fu
    INTEGER :: i, j, k

    do k = 1, local_N
       do j = 1, n(2)
          do i = 1, n(1)+2
             if (i <= n(1)) then
                tmppointer_r(i,j,k) = u(i,j,k)
             else
                tmppointer_r(i,j,k) = real(0.0_pr)
             end if
          end do
       end do
    end do

      
    CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
    
    call fftw_mpi_execute_dft_r2c(fwdplan, tmppointer_r, tmppointer_cx)
    
    fu(:,:,:) = tmppointer_cx(:,:,:)
    fu = fu/(PRODUCT(REAL(n,pr)))

    CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

  END SUBROUTINE fftfwd



!=======================================
! Inverse Fourier transform
!=======================================
   SUBROUTINE fftbwd(fu,u)
 !     use, intrinsic :: iso_c_binding
      USE global_variables 
      IMPLICIT NONE 
      INCLUDE "mpif.h"
      include 'fftw3-mpi.f03'

      COMPLEX(pr), DIMENSION (1:n(1)/2+1,1:n(2),1:local_N), INTENT(IN) :: fu
      real(pr), DIMENSION (1:n(1),1:n(2),1:local_N), INTENT(OUT) :: u
      integer :: i, j, k

      tmppointer_cx(:,:,:) = fu(:,:,:)

      call fftw_mpi_execute_dft_c2r(bwdplan, tmppointer_cx, tmppointer_r)
      do k = 1, local_N
         do j = 1, n(2)
            do i = 1, n(1)
               u(i,j,k) = tmppointer_r(i,j,k)
            end do
         end do
      end do
          



      CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

   END SUBROUTINE fftbwd

!==========================================

  SUBROUTINE init_fft
!==========================================             
!      use, intrinsic :: iso_c_binding
      USE global_variables
      IMPLICIT NONE 
      INCLUDE "mpif.h"
      include 'fftw3-mpi.f03'

      tmpdata_r = fftw_alloc_real(2*C_local_alloc)
      tmpdata_cx = fftw_alloc_complex(C_local_alloc)
   

      call c_f_pointer(tmpdata_r, tmppointer_r, [C_n(1)+2, C_n(2),C_local_N])
      call c_f_pointer(tmpdata_cx, tmppointer_cx, [C_n(1)/2+1, C_n(2),C_local_N])

      fwdplan = fftw_mpi_plan_dft_r2c_3d(C_n(3), C_n(2), C_n(1), tmppointer_r, tmppointer_cx, MPI_COMM_WORLD, FFTW_MEASURE)
      bwdplan = fftw_mpi_plan_dft_c2r_3d(C_n(3), C_n(2), C_n(1), tmppointer_cx, tmppointer_r, MPI_COMM_WORLD, FFTW_MEASURE)



    END SUBROUTINE init_fft

    

  

   !==================================================
   ! Forward Fourier transform for multiple dimensions
   !==================================================
  SUBROUTINE fftfwd_m(u, fu, m)
     USE global_variables 
      IMPLICIT NONE 
      INCLUDE "mpif.h"
      include 'fftw3-mpi.f03'
     INTEGER, INTENT(IN) :: m
     real(pr), DIMENSION (1:n(1),1:n(2),1:local_N, 1:m), INTENT(IN) :: u
     COMPLEX(pr), DIMENSION (1:n(1)/2+1,1:n(2),1:local_N, 1:m), INTENT(OUT) :: fu
     INTEGER :: nn
     DO nn = 1,m
        call fftfwd(u(:,:,:,nn), fu(:,:,:,nn))
     END DO
     
     
   END SUBROUTINE fftfwd_m
   !==================================================
   ! Forward Fourier transform for multiple dimensions
   !==================================================
   SUBROUTINE fftbwd_m(fu, u, m)
      USE global_variables 
      IMPLICIT NONE 
      INCLUDE "mpif.h"
      include 'fftw3-mpi.f03'
     INTEGER, INTENT(IN) :: m
     COMPLEX(pr), DIMENSION (1:n(1)/2+1,1:n(2),1:local_N, 1:m), INTENT(IN) :: fu
     real(pr), DIMENSION (1:n(1),1:n(2),1:local_N, 1:m), INTENT(OUT) :: u
     INTEGER :: nn
     DO nn = 1,m
        call fftbwd(fu(:,:,:,nn), u(:,:,:,nn))
     END DO

   END SUBROUTINE fftbwd_m

!=============================================
 SUBROUTINE fft_deallocate
   use, intrinsic :: iso_c_binding
   IMPLICIT NONE 
   INCLUDE "mpif.h"
   include 'fftw3-mpi.f03'
   call fftw_destroy_plan(fwdplan)
   call fftw_destroy_plan(bwdplan)
   call fftw_free(tmpdata_r)
   call fftw_free(tmpdata_cx)
   call fftw_mpi_cleanup()
 END SUBROUTINE  fft_deallocate

END MODULE
 
