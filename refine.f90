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

   integer :: startIndex, endIndex
   character(len=:), allocatable :: tempString1, tempString2

   



   !=============================================
   ! MPI
   !=============================================
   CALL MPI_INIT(Statinfo)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,Statinfo)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD,np,Statinfo)
   call fftw_mpi_init()
   if(rank==0) print*, "start"
   



   !!! instrucitons https://github.com/uniFabiB/refinement#readme
   ! basically
   !  salloc --time=0-0:30 --mem-per-cpu=1800M --ntasks=16 --account=rrg-bprotas
   !  change filename_initial
   !  make clean; make; srun




   filename_initial = "u_result_q4_n512_B018_iterend.nc"          ! without folderpath
   								  ! .nc does not matter







   ! extract resolution value: for filename_initial = "u_result_q5_n256_B014_iterend" >> res_pre = 256
   startIndex = index(filename_initial,'_n')   ! first occurance
   tempString1 = filename_initial(startIndex+2:)
   endIndex = index(tempString1,"_")
   tempString2 = tempString1(:endIndex-1)
   read(tempString2,*) res_pre
   if(rank==0) print*, "assigned res_pre=", res_pre

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
      

         if(rank==0) print*, "reading ..."
         if(filename_initial(len(filename_initial)-2:len(filename_initial)) == ".nc") then
            filename_initial = filename_initial(1:len(filename_initial)-3)
         end if
         CALL read_field_R3toR3_ncdf(UvecPre, "./"//filename_initial//".nc", "Ux", "Uy", "Uz")

         call mpi_barrier(mpi_comm_world, statinfo)

         if(read_nx /= res_pre) then
            if(rank==0) print*, "ERROR res_pre", res_pre, " probably not matching the files resolution ", read_nx, "STOPPING"
            errorIndicator = .true.
            exit
         end if
         
      
         
         call mpi_barrier(mpi_comm_world, statinfo)
         if(rank==0) print*, achar(9), "done"
         
         !CALL save_field_R3toR3_ncdf(UvecPre(:,:,:,1), UvecPre(:,:,:,2), UvecPre(:,:,:,3), "Ux", "Uy", "Uz", "./"//filename_initial//"_ie.nc", "netCDF")
      
         !if(rank==0) print*, "save done"
         
         if(rank==0) print*, "fourier ..."
         allocate( fuPre(1:n(1)/2+1,1:n(2),1:local_N, 1:3) )
         DO nn = 1,3
            call fftfwd(uvecPre(:,:,:,nn), fuPre(:,:,:,nn))
         END DO
         
         call mpi_barrier(mpi_comm_world, statinfo)
         if(rank==0) print*, achar(9), "done"
         
         if(rank==0) print*, "saving dat file ..."
         call save_velocity_cx(fuPre, "./"//filename_initial//"_cx.dat")
         call mpi_barrier(mpi_comm_world, statinfo)
         if(rank==0) print*, achar(9), "done"
         
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


         if(rank==0) print*, "starting ic_refine ..."
         call initial_condition_refine((/res_pre, res_pre, res_pre/), "./"//filename_initial//"_cx.dat")
         call mpi_barrier(mpi_comm_world, statinfo)
         if(rank==0) print*, achar(9), "ic_refine done"

         !if(rank==0) print*, "UvecRefined(1,2,1,1)", UvecRefined(1,2,1,1), "UvecRefined(3,2,1,3)", UvecRefined(3,2,1,3), "UvecRefined(2,2,2,2)", UvecRefined(2,2,2,2)
         
         if(rank==0) print*, "saving nc file"
         CALL save_field_R3toR3_ncdf(UvecRefined(:,:,:,1), UvecRefined(:,:,:,2), UvecRefined(:,:,:,3), "Ux", "Uy", "Uz", "./"//filename_initial//"_refined"//fb_int_to_str(res_pre)//"->"//fb_int_to_str(res_refined)//".nc", "netCDF")
         call mpi_barrier(mpi_comm_world, statinfo)
         if(rank==0) print*, achar(9), "done"

         if(rank==0) print*, "deleting dat file"
         if(rank==0 .and. .not.errorIndicator) then
            open(unit=1234, iostat=deleteStat, file="./"//filename_initial//"_cx.dat", status='old')
            if (deleteStat == 0) then
               close(1234, status='delete')
               print*, achar(9), "done"
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



