Name of student- Kamalesh Bhambare 

Compilation -  icpc -openmp main.cpp
Execution - ./a.out

This program reads the radius of sphere and radi of shell along with other Monte Carlo simulation parameters from the input file. Input file also gives required number of threads to be created. The number of processors are decided by using the pbs script file.  

Code reads all input (along with the number of thread) from <input.dat> file and write the output in <output.dat>. It performs necessary check for the input variables.

Program parallizes the for loop in the Monte Carlo simulation used for volume calculation. Both sphere and shell volume program are parallalized. Program uses srand() function of C++ for pseudo random number genration. 