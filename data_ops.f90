MODULE data_ops

  IMPLICIT NONE

  CONTAINS

        !=================================
        !    SAVE VELOCITY FROM NS SYSTEM
        !=================================
        SUBROUTINE save_velocity_cx(u_cx, filename)
          USE global_variables  
          use mpi
          IMPLICIT NONE
          
          complex(pr), DIMENSION(1:n(1)/2+1,1:n(2),1:local_N,1:3), INTENT(IN) :: u_cx
          CHARACTER(*), intent(in) :: filename
          complex(pr), DIMENSION(:,:,:), allocatable :: aux_cx

          integer :: local_count, j1, j2, j3

          local_count = (n(1)/2+1)*n(2)*local_N
  
          
          if(rank == 0) then
             allocate(aux_cx(1:n(1)/2+1, 1:n(2), 1:n(3)))
          else
            allocate(aux_cx(1:1,1:1,1:1)) ! otherwise the compiler complains that it is not allocated but called on the ranks where nothing happens
          end if
          call mpi_barrier(mpi_comm_world, statinfo)
          call mpi_gather(u_cx(:,:,:,1), local_count, MPI_DOUBLE_COMPLEX, aux_cx, local_count, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, statinfo)
          
          if (rank == 0) then
             open(11, file = filename, status = 'replace')
             do j3 = 1,n(3)
                do j2 = 1,n(2)
                   do j1 = 1, n(1)/2+1
                      write(11, "(2 G20.12)") real(aux_cx(j1,j2,j3)), aimag(aux_cx(j1,j2,j3))
                   end do
                end do
             end do
             close(11)
             aux_cx = cmplx(0.0_pr)
          end if
          call mpi_barrier(mpi_comm_world, statinfo)
          call mpi_gather(u_cx(:,:,:,2), local_count, mpi_double_complex, aux_cx, local_count, mpi_double_complex, 0, mpi_comm_world, statinfo)
          
          if (rank == 0) then
             open(11, file = filename, status = 'old', position = 'append')
             do j3 = 1,n(3)
                do j2 = 1,n(2)
                   do j1 = 1, n(1)/2+1
                      write(11, "(2 G20.12)") real(aux_cx(j1,j2,j3)), aimag(aux_cx(j1,j2,j3))
                   end do

                end do
             end do
             close(11)
             aux_cx = cmplx(0.0_pr)
          end if
          call mpi_barrier(mpi_comm_world, statinfo)
          call mpi_gather(u_cx(:,:,:,3), local_count, mpi_double_complex, aux_cx, local_count, mpi_double_complex, 0, mpi_comm_world, statinfo)
          
          if (rank == 0) then
             open(11, file = filename, status = 'old', position = 'append')
             do j3 = 1,n(3)
                do j2 = 1,n(2)
                   do j1 = 1, n(1)/2+1
                      write(11, "(2 G20.12)") real(aux_cx(j1,j2,j3)), aimag(aux_cx(j1,j2,j3))
                   end do
                end do
             end do
             close(11)
             !call system('gzip -f '//filename)
             deallocate(aux_cx)
          end if
          call mpi_barrier(mpi_comm_world, statinfo)



          
          
        END SUBROUTINE save_velocity_cx


         !==========================================
        ! SAVE FIELD IN R3
        !==========================================
        SUBROUTINE save_field_R3toR3_ncdf(f1, f2, f3, f1_name, f2_name, f3_name, file_name, myformat)
          USE global_variables
          USE netcdf
          use mpi
          IMPLICIT NONE

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N), INTENT(IN) :: f1, f2, f3
          CHARACTER(len=*) :: file_name
          CHARACTER(len=*) :: f1_name, f2_name, f3_name
          CHARACTER(len=*) :: myformat
         
          REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: u1, u2, u3
          REAL(pr), DIMENSION(:,:,:,:), ALLOCATABLE :: myfield

          REAL(pr) :: local_maxf1, local_maxf2, local_maxf3, maxf1, maxf2, maxf3 
          
          INTEGER :: ncout, ncid, varids(3), dimids(3), f_id
          INTEGER :: x_dimid, y_dimid, z_dimid, ux_id, uy_id, uz_id, uvec_id, maxUx_id, maxUy_id, maxUz_id       

          CHARACTER(200) :: parallel_file
          CHARACTER(2) :: RANKtxt
          INTEGER :: fname_len, ii
          INTEGER, DIMENSION(1:3) :: starts, counts
          

          IF (parallel_data) THEN
             
             IF (rank==0) THEN
                ncout = nf90_create(file_name, NF90_CLOBBER, ncid=ncid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_def_dim(ncid, "x", n(1), x_dimid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_def_dim(ncid, "y", n(2), y_dimid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_def_dim(ncid, "z", NF90_UNLIMITED, z_dimid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                dimids = (/ x_dimid, y_dimid, z_dimid /)

                ncout = nf90_def_var(ncid, TRIM(f1_name), NF90_DOUBLE, dimids, ux_id)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_def_var(ncid, TRIM(f2_name), NF90_DOUBLE, dimids, uy_id)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_def_var(ncid, TRIM(f3_name), NF90_DOUBLE, dimids, uz_id)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                     
                ncout = nf90_enddef(ncid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_close(ncid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
             END IF
             CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
             starts = (/ 1, 1, local_k_offset+1 /)

            
             
             counts = (/ n(1), n(2), local_N /)            
 
             !!--------------------------
             !! START netCDF ROUTINES
             !!--------------------------
             DO ii=0,np-1
                IF (rank==ii) THEN 
                   ncout = nf90_open(file_name, NF90_WRITE, ncid)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

                   ncout = nf90_inq_varid(ncid, TRIM(f1_name), f_id)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   ncout = nf90_put_var(ncid, f_id, f1, start = starts, count = counts)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   
                   ncout = nf90_inq_varid(ncid, TRIM(f2_name), f_id)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   ncout = nf90_put_var(ncid, f_id, f2, start = starts, count = counts)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   
                   ncout = nf90_inq_varid(ncid, TRIM(f3_name), f_id)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   ncout = nf90_put_var(ncid, f_id, f3, start = starts, count = counts)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

                   ncout = nf90_close(ncid)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                END IF
                CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo) 
             END DO
             
             CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
 
          ELSE 
             IF (rank == 0) THEN 
                ALLOCATE( u1(1:n(1),1:n(2),1:n(3)) )
                ALLOCATE( u2(1:n(1),1:n(2),1:n(3)) )
                ALLOCATE( u3(1:n(1),1:n(2),1:n(3)) )
             END IF
             CALL MPI_GATHER(f1, total_local_size, MPI_DOUBLE_PRECISION, u1, total_local_size, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, Statinfo)
             CALL MPI_GATHER(f2, total_local_size, MPI_DOUBLE_PRECISION, u2, total_local_size, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, Statinfo)
             CALL MPI_GATHER(f3, total_local_size, MPI_DOUBLE_PRECISION, u3, total_local_size, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, Statinfo)

             IF (rank==0) THEN
                ncout = nf90_create(file_name, NF90_CLOBBER, ncid)
                ncout = nf90_def_dim(ncid, "x", n(1), x_dimid)
                ncout = nf90_def_dim(ncid, "y", n(2), y_dimid)
                ncout = nf90_def_dim(ncid, "z", n(3), z_dimid)
                dimids =  (/ x_dimid, y_dimid, z_dimid /)
                
                ncout = nf90_def_var(ncid, TRIM(f1_name), NF90_DOUBLE, dimids, ux_id)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_def_var(ncid, TRIM(f2_name), NF90_DOUBLE, dimids, uy_id)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_def_var(ncid, TRIM(f3_name), NF90_DOUBLE, dimids, uz_id)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_enddef(ncid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
          
                ncout = nf90_put_var(ncid, ux_id, u1)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_put_var(ncid, uy_id, u2)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_put_var(ncid, uz_id, u3)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

                ncout = nf90_close(ncid)
             
                DEALLOCATE( u1 ) 
                DEALLOCATE( u2 )
                DEALLOCATE( u3 )

             END IF  
             CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

          END IF 
 
        END SUBROUTINE save_field_R3toR3_ncdf


        !============================================================
        ! READ VORTICITY IN netCDF FORMAT
        !============================================================
        SUBROUTINE read_field_R3toR3_ncdf(myfield, filename, Fx_txt, Fy_txt, Fz_txt)
          USE global_variables
          USE netcdf
          USE mpi
          IMPLICIT NONE

          REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(OUT) :: myfield
          CHARACTER(len=*), INTENT(IN) :: filename, Fx_txt, Fy_txt, Fz_txt
         
          REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: local_f, global_f
 
          INTEGER :: ncout, ncid, fid, dimids(3)
          INTEGER :: x_dimid, y_dimid, z_dimid
          INTEGER :: fname_len, ii, nx_ncdf, ny_ncdf, nz_ncdf

          INTEGER, DIMENSION(1:3) :: starts, counts

          
          
          IF (parallel_data) THEN
             ALLOCATE( local_f(1:n(1),1:n(2),1:local_N) )
             starts = (/ 1, 1, local_k_offset+1 /)
             
             counts = (/ n(1), n(2), local_N /)            
 
             !--------------------------
             ! START netCDF ROUTINES
             !--------------------------
             DO ii=0,np-1
                IF (rank == ii) THEN 
                   ncout = nf90_open(filename, NF90_NOWRITE, ncid)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

                   ncout = nf90_inq_dimid(ncid, "x", x_dimid)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   ncout = nf90_inq_dimid(ncid, "y", y_dimid)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   ncout = nf90_inq_dimid(ncid, "z", z_dimid)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

                   ncout = nf90_inquire_dimension(ncid, x_dimid, len = nx_ncdf)
                   read_nx = nx_ncdf
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   ncout = nf90_inquire_dimension(ncid, y_dimid, len = ny_ncdf)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   ncout = nf90_inquire_dimension(ncid, z_dimid, len = nz_ncdf)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

                   ncout = nf90_inq_varid(ncid, Fx_txt, fid)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   ncout = nf90_get_var(ncid, fid, local_f, start = starts, count = counts)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   myfield(:,:,:,1) = local_f
 
                   ncout = nf90_inq_varid(ncid, Fy_txt, fid)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   ncout = nf90_get_var(ncid, fid, local_f, start = starts, count = counts)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   myfield(:,:,:,2) = local_f
 
                   ncout = nf90_inq_varid(ncid, Fz_txt, fid)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   ncout = nf90_get_var(ncid, fid, local_f, start = starts, count = counts)
                   IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                   myfield(:,:,:,3) = local_f
 
                   ncout = nf90_close(ncid)
 
                   DEALLOCATE(local_f)
                END IF
                CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
             END DO

             CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

          ELSE
             IF (rank == 0) THEN
                ALLOCATE( global_f(1:n(1),1:n(2),1:n(3)) )
             END IF
            
             ALLOCATE( local_f(1:n(1),1:n(2),1:local_N) )

             IF (rank == 0) THEN 
                ncout = nf90_open(filename, NF90_NOWRITE, ncid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

                ncout = nf90_inq_dimid(ncid, "x", x_dimid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_inq_dimid(ncid, "y", y_dimid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_inq_dimid(ncid, "z", z_dimid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

                ncout = nf90_inq_varid(ncid, Fx_txt, fid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_get_var(ncid, fid, global_f)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
             END IF
             CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
             CALL MPI_SCATTER(global_f, total_local_size, MPI_DOUBLE_PRECISION, local_f, total_local_size, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, Statinfo)
             myfield(:,:,:,1) = local_f
 
             IF (rank == 0) THEN 
                ncout = nf90_inq_varid(ncid, Fy_txt, fid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_get_var(ncid, fid, global_f)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
             END IF
             CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
             CALL MPI_SCATTER(global_f, total_local_size, MPI_DOUBLE_PRECISION, local_f, total_local_size, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, Statinfo)
             myfield(:,:,:,2) = local_f
 
             IF (rank == 0) THEN 
                ncout = nf90_inq_varid(ncid, Fz_txt, fid)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_get_var(ncid, fid, global_f)
                IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                ncout = nf90_close(ncid)
             END IF
             CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
             CALL MPI_SCATTER(global_f, total_local_size, MPI_DOUBLE_PRECISION, local_f, total_local_size, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, Statinfo)
             myfield(:,:,:,3) = local_f
 
             IF (rank == 0) THEN
                DEALLOCATE( global_f )
             END IF
             DEALLOCATE( local_f )
             CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
          END IF  

        END SUBROUTINE read_field_R3toR3_ncdf



        !============================================================================================
!============================================================================================
        SUBROUTINE read_field_R3toR3_ncdf2(myfield, filename, Fx_txt, Fy_txt, Fz_txt)
         USE global_variables
         USE netcdf
         use mpi
         IMPLICIT NONE

         REAL(pr), DIMENSION(1:n(1),1:n(2),1:local_N,1:3), INTENT(OUT) :: myfield
         CHARACTER(len=*), INTENT(IN) :: filename, Fx_txt, Fy_txt, Fz_txt
        
         !REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: u

         INTEGER :: ncout, ncid, fid, dimids(3)
         INTEGER :: x_dimid, y_dimid, z_dimid
         INTEGER :: fname_len, ii, nx_ncdf, ny_ncdf, nz_ncdf, nn

         INTEGER, DIMENSION(1:3) :: starts, counts

         if  (rank == 0) then
            ncout = nf90_open(filename, NF90_NOWRITE, ncid)
            IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

            ncout = nf90_inq_dimid(ncid, "x", x_dimid)
            IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
            ncout = nf90_inq_dimid(ncid, "y", y_dimid)
            IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
            ncout = nf90_inq_dimid(ncid, "z", z_dimid)
            IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

            ncout = nf90_inquire_dimension(ncid, x_dimid, len = nx_ncdf)
            IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
            ncout = nf90_inquire_dimension(ncid, y_dimid, len = ny_ncdf)
            IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
            ncout = nf90_inquire_dimension(ncid, z_dimid, len = nz_ncdf)
            IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

            global_u = 0.0_pr  
            !allocate(u(1:n(1), 1:n(2), 1:n(3)))
         end if
         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
         do nn = 1,3
            if (rank == 0) then
               do ii = 0, np-1
                  starts = (/ 1, 1, local_N*ii+1 /)
                  counts = (/ n(1), n(2), local_N /)
                  select case (nn)
                  case (1)
                     ncout = nf90_inq_varid(ncid, Fx_txt, fid)
                     IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                     ncout = nf90_get_var(ncid, fid, global_u(:,:,local_N*ii+1:local_N*(ii+1)), start = starts, count = counts)
                     IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)

                  case (2)
                     ncout = nf90_inq_varid(ncid, Fy_txt, fid)
                     IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                     ncout = nf90_get_var(ncid, fid, global_u(:,:,local_N*ii+1:local_N*(ii+1)), start = starts, count = counts)
                     IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                     
                  case (3)
                     ncout = nf90_inq_varid(ncid, Fz_txt, fid)
                     IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                     ncout = nf90_get_var(ncid, fid, global_u(:,:,local_N*ii+1:local_N*(ii+1)), start = starts, count = counts)
                     IF (ncout /= NF90_NOERR) CALL ncdf_error_handle(ncout)
                  end select
               end do
               
               
            end if
                  
            CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
            CALL MPI_SCATTER(global_u, total_local_size, MPI_DOUBLE_PRECISION, myfield(:,:,:,nn), total_local_size, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, Statinfo)
         end do
            

         if (rank == 0) then
            ncout = nf90_close(ncid)
            !deallocate(u)
         end if
                  
         CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)  
                        
                           
                           
             
          END SUBROUTINE read_field_R3toR3_ncdf2


        !===============================================
        ! NETCDF ERROR HANDLE ROUTINE
        !===============================================
        SUBROUTINE ncdf_error_handle(nerror)
          USE global_variables
          USE netcdf
          IMPLICIT NONE

          INTEGER, INTENT(IN) :: nerror
          CHARACTER(80) :: error_string
          CHARACTER(2) :: K0txt, E0txt
 
          error_string = NF90_STRERROR(nerror)
 
          OPEN(10, FILE="./LOGFILES/maxET_info.log", STATUS='OLD', POSITION='APPEND')
          WRITE(10,*) " Error reading netCDF file. "//error_string
          CLOSE(10)
  
        END SUBROUTINE ncdf_error_handle



END MODULE 
