Name - Kamalesh S. Bhambare

1. To compile the program
g++ main.cpp

2. To run
./a.out

3. If the input.dat file doesnot exist code will report error. Code also checks for appropriateness of input for non-negative, zero values and terminates the program if error found. 

4. Program solves for heat equation with heat source for a cube using finite difference method. Dirichlet boundary condition was used at the boundary surface with uniform temperature (zero). Program generates creates grid from the input number of cells and generates a algebraic equation at each grid point. Set of algebraic equation generates a matrix. This matrix is further solved using selected iterative method.

5. Options for selction of iterative scheme
	Point Jacobi = 0
 	Gauss Seidel = 1 
	SOR  = 2

6. When SOR is selcted, Wopt is calculated using standard formula in the program using know mesh size. Its reported in the output file.

7. Inputs required are number of cells in each direction, heat source per cell and conductivity. Number of cells must be greater than one. Conductivity has be to greater than zero.  Heat source has to positive. Number of iterations must be positive. Option for selection of iterative scheme has to be either of 0,1,2.

8. Program report the number of iterations required, maximum error at the end the of iterations and execution time. 

9. Temperature at three different locations is reported in the program.

10. Scheme Comparison

			Number of iterations	Execution time (s)	T(0.5,0.5,0.5)	T(0.25,0.25,0.25)	T(0.75,0.75,0.75)
Point Jacobi		146 			2.96509			0.054917 	 0.0291564		0.0291564
Gauss-Seidel		79			1.53037			0.0549174	 0.0291565		0.0291566
SOR			26			0.533872		0.0549177	 0.0291566	  	0.0291567
