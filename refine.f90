!-----------------------------------------------------!
! Program used to solve the Cauchy problem            !
! in the  Navier-Stokes system                        !
!                                                     !
!                                                     !
! Parallel version, Complex-Complex FFT               !   
!                                                     !
! July, 2014.                                         !
!                                                     !
! Author: Diego Ayala                                 !
! Department of Mathematics and Statistics            !
! McMaster University                                 !
!-----------------------------------------------------!

PROGRAM main

   USE global_variables
   USE data_ops
   use fftwfunction
   USE function_ops
   use mpi
   IMPLICIT NONE
   INCLUDE "fftw3-mpi.f03"             ! Needed, as there is a fftw_mpi_local_size_3d command below

   COMPLEX(pr), DIMENSION (:,:,:,:), ALLOCATABLE :: fuPre
   integer :: res, res_pre, res_refined, nn, step, deleteStat

   



   !=============================================
   ! MPI
   !=============================================
   CALL MPI_INIT(Statinfo)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,Statinfo)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD,np,Statinfo)
   call fftw_mpi_init()
   if(rank==0) print*, "start"
   






   filename_initial = "u_result_0613_q4_B026_iter0800"          ! without folderpath or .nc
   res_pre = 256







   res_refined = res_pre * 2
   do step = 1,2
      if(step==1) then
         res = res_pre
         C_n = res
         n = res
      
         ALLOCATE ( K1(1:n(1)) )
         ALLOCATE ( K2(1:n(2)) )
         ALLOCATE ( K3(1:n(3)) )

         CALL initialize
         
         ALLOCATE( UvecPre(1:n(1),1:n(2),1:local_N,1:3) )

         CALL init_fft
      
      
         if(rank==0) print*, "n(:)", n(:)
      
         CALL read_field_R3toR3_ncdf(UvecPre, "./"//filename_initial//".nc", "Ux", "Uy", "Uz")

         call mpi_barrier(mpi_comm_world, statinfo)

         if(read_nx /= res_pre) then
            if(rank==0) print*, "ERROR res_pre", res_pre, " probably not matching the files resolution ", read_nx, "STOPPING"
            errorIndicator = .true.
            exit
         end if
         
      
         
         call mpi_barrier(mpi_comm_world, statinfo)
         if(rank==0) print*, "read done"
         
         !CALL save_field_R3toR3_ncdf(UvecPre(:,:,:,1), UvecPre(:,:,:,2), UvecPre(:,:,:,3), "Ux", "Uy", "Uz", "./"//filename_initial//"_ie.nc", "netCDF")
      
         !if(rank==0) print*, "save done"
         
         allocate( fuPre(1:n(1)/2+1,1:n(2),1:local_N, 1:3) )
         DO nn = 1,3
            call fftfwd(uvecPre(:,:,:,nn), fuPre(:,:,:,nn))
         END DO
         
         call mpi_barrier(mpi_comm_world, statinfo)
         if(rank==0) print*, "fourier done"
         
         call save_velocity_cx(fuPre, "./"//filename_initial//"_cx.dat")
         
         deallocate( UvecPre )
         if (rank==0) deallocate( global_u )
         deallocate( K1 )
         deallocate( K2 )
         deallocate( K3 )
         
         call mpi_barrier(mpi_comm_world, statinfo)
         
         if (rank == 0) then
            if(errorIndicator) then
               print *, "THERE WAS AN ERROR SOMEWHERE"
            else
               print *, "step 1 exited normally"
            end if
         end if

      elseif(step == 2) then
         
         res = res_refined
         C_n = res
         n = res

         ALLOCATE ( K1(1:n(1)) )
         ALLOCATE ( K2(1:n(2)) )
         ALLOCATE ( K3(1:n(3)) )

         CALL initialize

         ALLOCATE( UvecRefined(1:n(1),1:n(2),1:local_N,1:3) )

         CALL init_fft

         if(rank==0) print*, "n(:)", n(:)


         call initial_condition_refine((/res_pre, res_pre, res_pre/), "./"//filename_initial//"_cx.dat")

         call mpi_barrier(mpi_comm_world, statinfo)
         if(rank==0) print*, "ic_refine done"

         if(rank==0) print*, "UvecRefined(1,2,1,1)", UvecRefined(1,2,1,1), "UvecRefined(3,2,1,3)", UvecRefined(3,2,1,3), "UvecRefined(2,2,2,2)", UvecRefined(2,2,2,2)
         
         CALL save_field_R3toR3_ncdf(UvecRefined(:,:,:,1), UvecRefined(:,:,:,2), UvecRefined(:,:,:,3), "Ux", "Uy", "Uz", "./"//filename_initial//"_refined"//fb_int_to_str(res_pre)//"->"//fb_int_to_str(res_refined)//".nc", "netCDF")
         if(rank==0) print*, "save_ref done"


         if(rank==0 .and. .not.errorIndicator) then
            open(unit=1234, iostat=deleteStat, file="./"//filename_initial//"_cx.dat", status='old')
            if (deleteStat == 0) then
               close(1234, status='delete')
               print*, "deleting temporary file done"
            else
               errorIndicator = .true.
            end if
         end if

         deallocate( UvecRefined )
         if (rank==0) deallocate( global_u )
         deallocate( K1 )
         deallocate( K2 )
         deallocate( K3 )
         
         if (rank == 0) then
            if(errorIndicator) then
               print *, "THERE WAS AN ERROR SOMEWHERE"
            else
               print *, "step 2 exited normally"
            end if
         end if

      end if

   end do
   
   CALL MPI_FINALIZE (Statinfo)
   if(rank==0) print*, "end"

END PROGRAM main



