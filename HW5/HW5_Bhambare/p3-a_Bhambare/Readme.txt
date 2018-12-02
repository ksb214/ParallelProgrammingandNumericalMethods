Name of student- Kamalesh Bhambare 

Compilation - ifort -c ran5.f 
              icpc -openmp main_03.cpp ran5.o

Execution - ./a.out

This program reads the radius of sphere and radi of shell along with other Monte Carlo simulation parameters from the input file. Input file also gives required number of threads to be created. The number of processors are decided by using the pbs script file.  

Code reads all input (along with the number of thread) from <input.dat> file and write the output in <output.dat>. It performs necessary check for the input variables.

Program parallizes the for loop in the Monte Carlo simulation used for volume calculation. Both sphere and shell volume program are parallalized. 

For psuedo random generations thread id is used to seed the random number generator. This ensures that each thread has different seed. The program used for psuedo random number generations is ran5.f. Fortran program libray file is linked to the main c++ program using extern method. Use of program rather than c++ function for random number generators reduces program overheads in calling intrinsic function srand used in earlier part 2.
