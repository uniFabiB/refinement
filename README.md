# Instructions
1. Load the needed modules
   
   `module load fftw-mpi netcdf-fortran;`
1. Start an interactive session
   
   `salloc --time=0-0:30 --mem-per-cpu=1800M --ntasks=16 --account=rrg-bprotas`
1. Copy the `nc` file to refine into the working directory
1. in `refine.f90` change
   1. `filename_initial = "u_result_q4_n256_B026_iter0800"` to the copied file name
1. `make; srun`
