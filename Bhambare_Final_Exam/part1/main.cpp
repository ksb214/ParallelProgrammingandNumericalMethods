#include<iostream>
#include<mpi.h>
#include<fstream>
using namespace std;


main(int argc, char *argv[])
{
  int np, rank; // Number of processes and rank of process   
  int next, prev;
  int rrank;   // rrank is a variable that stores the recieved rank 

  char file_name[] = "out.?";

  ofstream outfile1;

  MPI_Status status;

  // Initializing MPI
  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &np);  // Determining the number of processes
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Determining the rank of the process

  if( np == 1) {
    cout<<"Only one process is specified \n No communication \n Exiting...\n";
    exit (0);
}

  // Finding next and previous in a cyclic manner
  next = (rank+1) % np ;
  prev = (rank+np -1) % np ;


  if ( (rank % 2) == 1 ) {   // This is for odd numbered ranks
    MPI_Send(&rank, 1, MPI_INT, next, 1, MPI_COMM_WORLD);
    MPI_Recv(&rrank, 1, MPI_INT, prev, 1, MPI_COMM_WORLD, &status);
  }
  else {  
    MPI_Recv(&rrank, 1, MPI_INT, prev, 1, MPI_COMM_WORLD, &status);
    MPI_Send(&rank, 1, MPI_INT, next, 1, MPI_COMM_WORLD);
     }

  //Writing output to file

  file_name[4] = rank + '0';

  outfile1.open(file_name,ios::out);
  outfile1<<" Process " << rrank << " comminicated successfully with process "<<rank << endl;
  outfile1.close();

  MPI_Finalize();

}
