MODULE global_variables
  use, intrinsic :: iso_c_binding   ! Newly Added April 17, 2017
  IMPLICIT NONE
!   INCLUDE "mpif.h"
!   include 'fftw3-mpi.f03'
 
  INTEGER, PARAMETER :: pr = KIND (1.0d0)


  integer :: RESOL = 128
  CHARACTER(len=:), allocatable :: filename_initial
  !- Run on graham

  LOGICAL, parameter :: save_binary2nc = .FALSE.


  LOGICAL :: kappaTest
  LOGICAL :: toDealias
  LOGICAL :: timing 
  LOGICAL :: save_diag_NS
  LOGICAL :: save_data_NS
  LOGICAL :: calc_geom_NS
  LOGICAL :: calc_ExactSol
  LOGICAL :: save_diag_Constr
  LOGICAL :: save_data_Constr
  LOGICAL :: save_diag_Optim
  LOGICAL :: save_data_Optim
  LOGICAL :: save_diag_lineMin
  LOGICAL :: save_data_lineMin
  LOGICAL :: parallel_data
  LOGICAL :: save_null_vortex

  logical :: errorIndicator = .false.
  integer :: read_nx = 0


  INTEGER, DIMENSION(3), SAVE :: n
  INTEGER, SAVE :: n_dim  
  INTEGER, SAVE :: DT_index, NU_index, ConsType, iniIndex 
  REAL(pr), SAVE :: E0, K0, PI, visc, dV, Kcut, Kmax
  Integer :: kkmax
  !--NOTE: Kcut = cut frequency used for dealiasing. 
  !-       Kmax = maximal frequency present in solution.
  !-       kkmax = ceiling(sqrt(real(n(1),pr)**2/4_pr+real(n(2),pr)**2/4_pr+real(n(3),pr)**2/4_pr)) 

  REAL(pr), DIMENSION (:), ALLOCATABLE, SAVE :: K1, K2, K3
  REAL(pr), DIMENSION (:), ALLOCATABLE, SAVE :: K1_filter, K2_filter, K3_filter ! added 26/09/21
  REAL(pr), DIMENSION (:), ALLOCATABLE, SAVE :: spectral_k ! added 28/09/21
  REAL(pr), DIMENSION (:,:,:,:), ALLOCATABLE, SAVE :: uvecPre, uvecRefined
  
  !========================================================================== 
  !                            MPI VARIABLES
  !==========================================================================
  INTEGER, SAVE :: rank, Statinfo, np
!  INTEGER, SAVE :: local_nlast, local_last_start, local_nlast_after_trans 
!  INTEGER, SAVE :: local_last_start_after_trans, total_local_size 
!  INTEGER, SAVE :: local_start, local_n

  INTEGER, SAVE :: local_nlastFixres, local_last_startFixres, local_nlast_after_transFixres
  INTEGER, SAVE :: local_last_start_after_transFixres, total_local_sizeFixres, local_startFixres, local_nFixres

  INTEGER, SAVE :: local_nlastHighres, local_last_startHighres, local_nlast_after_transHighres
  INTEGER, SAVE :: local_last_start_after_transHighres, total_local_sizeHighres, local_startHighres, local_nHighres


  !========================================================================== 
  !                            FFTW VARIABLES
  !==========================================================================



  INTEGER(C_INTPTR_T), DIMENSION(3), SAVE :: C_n
  INTEGER(C_INTPTR_T), SAVE :: C_local_alloc, C_local_N, C_local_k_offset       ! Newly Added July 14, 2017
  INTEGER, SAVE :: local_alloc, local_N, local_k_offset, total_local_size       ! Newly Added July 14, 2017


  INTEGER, SAVE :: final_time_iter, reclen
  REAL(pr), DIMENSION (:,:,:,:), ALLOCATABLE, SAVE :: fwd_Field1, fwd_Field2, adj_Uvec, adj_Wvec, Uvec0, adj_Uvec0, adj_Uvec0_direction
  real(pr), dimension(:,:,:), allocatable, save :: global_u




END MODULE
