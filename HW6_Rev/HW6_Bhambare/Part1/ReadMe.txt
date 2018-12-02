1. Name of student- Kamalesh Bhambare 

2. Compilation - icpc -openmp main_03.cpp -o pj.out 

3. Number of processors - Provided in pbs_pj.sh script file

4. Execution - qsub pbs_pj.sh

5. code status- operational

6. If the input.dat file doesnot exist code will report error. Code also checks for appropriateness of input for non-negative, zero values and terminates the program if error found. 

7. Program solves for heat equation with heat source for a cube using finite difference method. Dirichlet boundary condition was used at the boundary surface with uniform temperature (zero). Program generates creates grid from the input number of cells and generates a algebraic equation at each grid point. Set of algebraic equation generates a matrix. This matrix is further solved using selected iterative method.

8. Options for selction of iterative scheme
	Point Jacobi = 0
 	Gauss Seidel = 1 
	SOR  = 2

9. When SOR is selcted, Wopt is calculated using standard formula in the program using know mesh size. Its reported in the output file.

10. Inputs required are number of cells in each direction, heat source per cell and conductivity. Number of cells must be greater than one. Conductivity has be to greater than zero.  Heat source has to positive. Number of iterations must be positive. Option for selection of iterative scheme has to be either of 0,1,2.

11. Program report the number of iterations required, maximum error at the end the of iterations and execution time. 

12. Temperature at three different locations is reported in the program. 

13. Banded matrix approach is used in solving. At each PJ iteration TH number of threads are assigned to number of processors provide in the script file. 

