!=======================================================
! MODULE CONTAINS ROUTINES REPRESENTING OPERATIONS 
! APPLIED TO FUNCTIONS.
!
! (*) function_ops_allocate
! (*) function_ops_deallocate
! (*) derivative_fourier
! (*) G_sigma_fourier
! (*) laplacian_fourier
! (*) abs_deriv_fourier
! (*) vel2vort_fourier
! (*) vort2vel_fourier
! (*) L2_product
! (*) calculate_total_energy
! (*) calculate_spectrum
! (*) div_fourier
 

! (*) vel2vort
! (*) Energy
! (*) Enstrophy
! (*) laplacian
! (*) and more...
! (*) calculate_spectral
!=======================================================

MODULE function_ops
  use global_variables
  IMPLICIT NONE

  real(pr), dimension(:,:,:), allocatable :: Epoint ! needed in calculate_total_energy
  real(pr), dimension(:,:), allocatable :: vort_plane
  complex(pr), dimension(:,:,:), allocatable :: temp1_function_cx ! needed in vel2vort_fourier
  complex(pr), dimension(:,:,:), allocatable :: temp2_function_cx ! needed in vel2vort_fourier
  



CONTAINS
  !========================================================= 
  ! SUBROUTINE: function_ops_allocate()
  !
  ! allocate the variables needed for solver
  !=========================================================
  SUBROUTINE function_ops_allocate()
    use global_variables
    IMPLICIT NONE
    if (.not. allocated(Epoint)) allocate(Epoint(1:n(1), 1:n(2), 1:local_N))
    if (.not. allocated(temp1_function_cx)) allocate(temp1_function_cx(1:n(1)/2+1, 1:n(2), 1:local_N))
    if (.not. allocated(temp2_function_cx)) allocate(temp2_function_cx(1:n(1)/2+1, 1:n(2), 1:local_N))
    if (rank == 0) then
       if (.not. allocated(vort_plane)) allocate(vort_plane(1:n(1), 1:n(3)))
    end if
    
    
  END SUBROUTINE function_ops_allocate

  !========================================================= 
  ! SUBROUTINE: solver_deallocate()
  !
  ! deallocate the variables needed for solver
  !=========================================================
  SUBROUTINE function_ops_deallocate()
     use global_variables
    IMPLICIT NONE
    
    if(allocated(Epoint)) deallocate(Epoint)
    if(allocated(temp1_function_cx)) deallocate(temp1_function_cx)
    if(allocated(temp2_function_cx)) deallocate(temp2_function_cx)
    if(rank == 0) then
       if (allocated(vort_plane)) deallocate(vort_plane)
    end if
    
  END SUBROUTINE function_ops_deallocate

  
  !===================================
  !  Initial guess
  !===================================
  SUBROUTINE initial_condition(filename)
    USE global_variables
    USE data_ops
  
    
    IMPLICIT NONE
    INCLUDE "mpif.h"
    CHARACTER(*), intent(in) :: filename
    REAL(pr), DIMENSION(1:3) :: dx
    CHARACTER(2) :: Fx_txt, Fy_txt, Fz_txt
    INTEGER :: nn,i,j,k
    INTEGER :: ii,jj,kk
    if (rank == 0) print *, filename
    Fx_txt = "Ux"
    Fy_txt = "Uy"
    Fz_txt = "Uz"
    CALL read_field_R3toR3_ncdf2(Uvec0, filename, Fx_txt, Fy_txt, Fz_txt)
    
    
  END SUBROUTINE initial_condition



   !===================================
   !  my int to string
   !===================================
   function fb_int_to_str(i) result(str)
      integer, intent(in) :: i
      character(len=100) :: temp
      character(len=:), allocatable :: str
      write(temp, '(I0)') i
      str = trim(temp)
   end function fb_int_to_str

  !===================================
  !  Initial guess
  !===================================
  SUBROUTINE initial_condition_refine(n_pre, path)
    USE global_variables
    USE fftwfunction
    IMPLICIT NONE
    INCLUDE "mpif.h"

    integer, dimension(1:3), intent(in) :: n_pre
    complex(pr), dimension(:,:,:), allocatable :: global_cx
    character(*), intent(in) :: path
    real(pr) :: val1, val2
    integer :: j1, j2, j3, jj2, jj3, nn
    integer :: local_count
    integer :: ios
    integer :: filehandle
    logical :: openedOk

    local_count = (n(1)/2+1)*n(2)*local_N
    
    
    if (rank == 0) print*, achar(9), path

    if (rank == 0) then
       allocate(global_cx(1:n(1)/2+1, 1:n(2), 1:n(3)))
       print *, achar(9), "allocation OK!"
       open(filehandle, file=path, iostat=ios)
       INQUIRE( UNIT=filehandle, OPENED=openedOk ) 
       IF ( openedOk ) print*, achar(9), "file opened ok"
       if(ios /= 0) then
         print*, achar(9), "error opening file"
         errorIndicator = .true.
       end if
    end if
    allocate( temp1_function_cx(1:n(1), 1:n(2), 1:local_N) )
    
    if(rank==0) print*, achar(9), "lines in file", 3*n_pre(3)*n_pre(2)*(n_pre(1)/2+1)
    do nn = 1,3
      if(rank==0) print*, achar(9), "starting reading dat file ", nn, "/3", achar(9), "..."
       if(rank == 0) then
          global_cx = cmplx(0.0_pr)
          do j3 = 1, n_pre(3)
             do j2 = 1, n_pre(2)
                do j1 = 1, n_pre(1)/2+1
                   read(filehandle, "(2 G20.12)", iostat=ios) val1, val2
                   if(ios < 0) then
                     print*, achar(9), "error reading file: 0 /= iostat =", ios, "at j1=",j1
                     errorIndicator = .true.
                     exit
                   end if
                   jj2 = j2
                   jj3 = j3
                   if (j2 > n_pre(2)/2+1) jj2 = j2-1-n_pre(2)+n(2)+1
                   if (j3 > n_pre(3)/2+1) jj3 = j3-1-n_pre(3)+n(3)+1
                   global_cx(j1,jj2,jj3) = cmplx(val1, val2)
                end do
                if(ios < 0) then
                  print*, achar(9), "error reading file: 0 /= iostat =", ios, "at j2=",j2
                  errorIndicator = .true.
                  exit
                end if
             end do
             if(ios < 0) then
               print*, achar(9), "error reading file: 0 /= iostat =", ios, "at j3=",j3
               errorIndicator = .true.
               exit
             end if
          end do
       end if
       CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)
       if(rank==0) print*, achar(9), achar(9), "done"
       call MPI_Scatter(global_cx(:,:,:), local_count, MPI_DOUBLE_COMPLEX, temp1_function_cx, local_count, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, Statinfo)
       call fftbwd(temp1_function_cx, uvecRefined(:,:,:,nn))
    end do

    if (rank == 0) then
       close(filehandle)
       deallocate(global_cx)
       print *, achar(9), "initial_data_refine OK!"
    end if
    
    CALL MPI_BARRIER(MPI_COMM_WORLD, Statinfo)

    END SUBROUTINE initial_condition_refine


  
     

END MODULE function_ops
