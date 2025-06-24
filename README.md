# Instructions
1. Copy the `nc` file to refine into the working directory
2. in `refine.f90` change
   1. `filename_initial = "u_result_0613_q4_B026_iter0800"` to the copied file name
   2. `res_pre = 256` to the unrefined resolution
3. `make; mpiexec -n 32 prog`
