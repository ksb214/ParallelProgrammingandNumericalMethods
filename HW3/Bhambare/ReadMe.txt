1. To compile the program
g++ main.cpp

2. To run
./a.out

3. If the input.dat file doesnot exist code will report error. Code also checks for appropriateness of input for non-negative, zero values and terminates the program if error found. 

4. Program solves for heat equation with heat source for a cube using finite difference method. Dirichlet boundary condition was used at the boundary surface with uniform temperature (zero). Program generates creates grid from the input number of cells and generates a algebraic equation at each grid point. Set of algebraic equation generates a matrix. This matrix is further solved using Gaussian Elmination method in the header file func.h. 

5. Inputs required are number of cells in each direction, heat source per cell and conductivity. Number of cells must be greater than one. Conductivity has be to greater than zero.  Heat source has to positive.

5. Program also report the execution time. 