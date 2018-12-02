// Code written by Kamalesh Bhambare to solve heat equation with source in cube with banded matrices and a structure to store temperture 
// and location of grid point
#include<iostream>
#include<fstream>
#include<math.h>
#include<time.h>
#include<mpi.h>

using namespace std;

// Defining the structure for storing variables at each grid point
struct Element
{	
	double T;          
	double x;
	double y;
	double z;
};

void read_in();         // Reads input data 
void check_in(void);    // Checks input for negative and non zero values as well as for input file existance
void echoin(void);      // Echos input
void define_domain();   // Domain decomposition function

void grid(Element *);   // Calculate x and y for grid point
void init(Element *);   // Initialise the grid variables
void SetEqn(Element *); // Setting the coefficients of the matrix
void SolveEqn();
void create_array(void); //Creation and initialization of temperature

inline double err(double *, double *);

void print_node (Element *);  // For initial checking
void plot_node(Element *);
void print_output(void);

double *C, *W, *N, *E, *S, *T, *B, *Q; 
int nx, ny, nz, n, I, N1;
double nxd;
double dx,dy,dz;
double xlength, ylength, zlength;
double cnd, q, Tb;
int info;
time_t date;
double exec_time;
int east, west, top, bot, north, south;
bool flag;
double eps,omega,maxRes;              // Parameters used for iterative methods
int itr,opt,TH, TID,noit; 
double *Tn, *oldT;

//MPI variables
int n_proc, myrank;
int IN_STATUS; // Status for checking input. Set to zero if input is not correct else set to 1. 

int *proc; // Process index
int *k_min, *k_max; // Domain is divided into n_proc processors. Min and max of each domain belonging to process is identified using these two. 
int *left_proc, *right_proc; // process on left and right of process under consideration 

// Index calculator
int INDEX( int i, int j, int k){
return (i + (N1+2)*j + (N1+2)*(N1+2)*k );
}

void define_domains(void);  // Function to define domains, deciding the left and right for the given proces and bounds of the given process

void PJ_mpi(void);              // PJ solver
void free_array(void);      // Subroutine to free various arrays
void solve(void);           // Solver which calls PJ or SOR
void SOR_mpi(void);

ofstream out;
ofstream outfile_seteq;

int main(int argc, char *argv[]){
	clock_t start = clock();
	int p,d;

	MPI_Status status;
	MPI_Init(&argc, &argv);    // Initialising MPI

	MPI_Comm_size(MPI_COMM_WORLD, &n_proc);     
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	IN_STATUS = 1;

	//	Element *node;
	
	if(myrank == 0) {
		read_in();
		check_in();
		if(IN_STATUS != 0) {
			IN_STATUS = 1;
			echoin();           //master is writting the input echo to file             
			
			//	node = new Element[n+1];		
			//	init(node);   // Initialises the node   
			//grid(node);  // Generating grid with stored information of x,y,z of each node	
		}
	}
	else read_in();              // All child processes have the access to the input variables after c

	if(IN_STATUS == 0) exit (0);    //master will terminate the process if inputs are not corret

	MPI_Barrier(MPI_COMM_WORLD);    //Barrier untill all processes read the input

	create_array();
	define_domains();
	//	PJ_mpi();
	solve();



	//Sending the final T matrices to the process 0
	if( myrank!= 0 ) {
	  if( (left_proc[myrank] >= 0 && left_proc[myrank] < n_proc) || (right_proc[myrank] >= 0 && right_proc[myrank] < n_proc) ) {
	    d = k_max[myrank] - k_min[myrank] + 1;
	    MPI_Send(Tn + INDEX(0, 0, k_min[myrank]), (N1+2)*(N1+2)*d, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	    //	    printf("Message sent processor: %d, from: %ld, size: %ld\n",myrank, INDEX(0, 0, k_min[myrank]), (N+2)*(N+2)*d);
	  }
	}

	if( myrank == 0) {
	  for( p =1; p < n_proc; p++ )
	    {
	      if( (left_proc[p] >= 0 && left_proc[p] < n_proc) || (right_proc[p] >= 0 && right_proc[p] < n_proc) ) {
		//printf("p is :%d and n_proc: %d\n",p, n_proc);
		d = k_max[p] - k_min[p] + 1;
		MPI_Recv(Tn + INDEX(0, 0, k_min[p]), (N1+2)*(N1+2)*d, MPI_DOUBLE, p, 1, MPI_COMM_WORLD, &status);
		//printf("Message received Processor: %d, from: %ld, size: %ld\n",p, INDEX(0, 0, k_min[p]), (N+2)*(N+2)*d);
	      }// end of if
	    }// for loop 'p'
	  	  print_output();
	}// if myrank

     

	if(myrank == 0) 
	{
	  //	  for(int i=0; i<(n+1); i++ )
	  //  node[i].T = Tn[i];
	  
	  //  plot_node(node); 
	  //	  print_node(node);
	  //	  delete[] node;	  
	  exec_time = ((double)clock() - start)/CLOCKS_PER_SEC;	
	  out<<"Execution time = "<<exec_time << " seconds \n";
	  out.close(); 
	}

		
	free_array();
	MPI_Finalize();

} // main ends here


void read_in(){
	
	ifstream inp;
  	inp.open("input.dat",ios::in);
	
	inp >> nxd >> q >> cnd >> eps >> itr >> opt ;
  	inp.close();
  	
  	I = (int) nxd;
  	nx = I + 1 ;
	
	N1 = I - 1 ;
	ny = nx;
  	nz = nx; 	
  	n = nx*ny*nz;
  	xlength = 1; 
  	ylength = 1;
  	zlength = 1; 	
  	Tb = 273;	
  	omega = 2/(1+sin(3.14/I));
}

// Program to read and check input
void check_in(){
	
 	if( I <= 1) {
		cout<<"Number of computational cells should be greater than one\nExiting....\n";
		IN_STATUS = 0;
	}	
    
  	if( q <= 0) {
		cout<<" Negative or zero source\nExiting....\n";
		IN_STATUS = 0;
	}

	if( cnd <= 0) {
			cout<<"Thermal condictivity can't be negative or zero\nExiting....\n";
			IN_STATUS = 0;
	}

	if(I%4 != 0) {
		cout<<"Number of cells should be multiple of 4 to place grid at 0.25, 0.5 and 0.75 location\n";
		cout<<"Still running the code for "<< nx << " cells \n";
	}
	
	if( itr <= 0) {
			cout<<"Number of iterations can't be negative or zero\nExiting....\n";
			IN_STATUS = 0;
	}
	
	if( opt >2 || opt< 0) {
			cout<<"Not a valid input for iterative method selection, Use 0 for Point Jacobi, \n";
			cout<<"1 for Gauss-Seidel, 2 for SOR ....\n";
			IN_STATUS = 0;
	}
}

void echoin()
{
	out.open("output.dat",ios::out);
  	time(&date);	
  	out<<"Code Name: "<< "main" <<'\n';
  	out<<"Version Number: "<< "1.0" <<'\n';
  	out<<"Author Name: "<< "Kamalesh Bhambare " <<'\n';
  	out<<"Date and time of execution: " << ctime(&date) <<'\n';
  	
  	if(IN_STATUS!=0)out<<"Input check...successful.\n\n";
	
	out<<"Number of cells per dimension = "<< I <<'\n';
	out<<"Uniform distributed source strength = "<< q <<'\n';
	out<<"Thermal conductivity = "<< cnd <<'\n';
	out<<"Relative convergence criterion = "<< eps <<'\n';
	out<<"Number of iterations used during solution = "<< itr << '\n';
	out<<"Number of processors used during solution = "<< n_proc <<endl; 
	switch(opt){
	case 0:
	  out<<"Iterative method used = Point Jacobi \n";
	  break;
	case 1: 
	 out<<"Iterative method used = Gauss-Seidel \n";
	 break;
	case 2:
	 out<<"Iterative method used = SOR"<<endl;
	 out<<"Optimal SOR weight = "<< omega <<" \n";
	 break;
	}
}

void init (Element *node)
{
	int i;
	double size;
	size= (N1+2)*(N1+2)*(N1+2);
	for ( i=0; i<size; i++ )
	{
		node[i].T = 0.0;
		node[i].x = 0.0;
		node[i].y = 0.0;
		node[i].z = 0.0;
	}
}

void grid(Element *node)
{
	dx = xlength/I;
	dy = ylength/I;
	dz = zlength/I;
//	cout<< "dx used is "<< dx << " dy used is "<<dy <<'\n';
	int k;

	for(int l=1; l<=ny; l++)
	{
		for (int j=1; j<=nz; j++)
		{
			for (int i=1; i<=nx; i++)
			{
				//Calculating the node number corresponding to i,j
				k = i + (j-1)*nx + (l-1)*nx*nz;
				//cout<<k<<'\n';
				node[k].x = (i-1)*dx;
				node[k].y = (j-1)*dy; 
				node[k].z = (l-1)*dz; 
			}// End of i loop
		}//End of j loop
	} // End of l loop	
}
		

void create_array(void)
{
  long size, i;
  size = (N1+2)*(N1+2)*(N1+2);

  //dynamically allocating the size of Temp array T based on number of cells
  Tn = (double *) malloc ( sizeof(double)*( size ) ) ;
  oldT = (double *) malloc ( sizeof(double)*( size ) ) ;

  for( i=0; i<size; i++ )
    Tn[i] = oldT[i] = 0.0;

}


void free_array(void)
{
  free(Tn);
  free(oldT);

  free(proc);
  free(k_min);
  free(k_max);
  free(left_proc);
  free(right_proc);
}

void define_domains(void)
{
  int i, k, p;
  double eps, d, z_min, z_max;
  
  //Allocating arrays for process information
  proc = (int*) malloc( (N1 + 2) * sizeof(int) );
  k_min = (int*) malloc( n_proc * sizeof(int) );
  k_max = (int*) malloc( n_proc * sizeof(int));
  left_proc = (int*) malloc( n_proc * sizeof(int));
  right_proc = (int*) malloc( n_proc * sizeof(int));

  // divide zrange evenly among processes
  eps = 0.0001;
  d = (N1 - 1 + 2*eps) / n_proc;

  for (p = 0; p < n_proc; p++) {
    //process domain
    z_min = -eps + 1 + p*d;
    z_max = z_min + d;

    // identify zvertices belonging to the current process
    for (i = 1; i <= N1; i++)
      if (z_min <= i && i < z_max)proc[i] = p;
  }

  for (p = 0; p < n_proc; p++) {
    // find the smallest vertex index in domain
    for (i = 1; i <= N1; i++)
      if (proc[i] == p) break;
    k_min[p] = i;

    // find largest vertex index in domain
    for (i = N1; i >= 1; i--)
      if (proc[i] == p) break;
    k_max[p] = i;

    // find processes to left and right
    left_proc[p] = right_proc[p] = -1;
    if (proc[p] != -1) {
      if (k_min[p] > 1 && k_min[p] <= N1)left_proc[p] = proc[k_min[p] - 1];
      if (k_max[p] > 0 && k_max[p] < N1)right_proc[p] = proc[k_max[p] + 1];
    }
  }

  ofstream fout;
  //Writing the domain information (only for master)
  if( myrank == 0) {
    
    fout.open("domain.dat",ios::out);
    fout<<"Total number of processes is "<<n_proc<<endl;

    for( p=0; p<n_proc; p++)
      {
	fout<< "Process: "<< p << " k_min: " << k_min[p] << " k_max: "<< k_max[p] << " left_proc: "<< left_proc[p] <<"  right_proc: "<< right_proc[p]<<endl;
      }

    fout.close();
  }

}

void print_node(Element *node)
{
	int i;
	ofstream outfile;
  	outfile.open("WholeData.txt",ios::out);
	
	outfile << "i     X      Y       Z       T"<<'\n';
	for (i=1; i<=n; i++ )
	{
		outfile<<i<<'\t'<<node[i].x<<'\t'<<node[i].y<<'\t'<<node[i].z<<'\t';
		outfile<<node[i].T<<'\n';
	}

	outfile.close();
}

void plot_node(Element *node){

	for (int i=1; i<=n; i++ )
	{
		if(node[i].x== 0.5 && node[i].y==0.5 && node[i].z==0.5){
		         out.precision(15);
			 out<<"T(0.5,0.5,0.5) = "<< node[i].T<<'\n';
		}
       
		//for (int i=1; i<=n; i++ )
       
		if(node[i].x==0.25 && node[i].y==0.25 && node[i].z==0.25){
			out.precision(15); 
			out<<"T(0.25,0.25,0.25) = "<< node[i].T<<'\n';
		}
	       	
		//	for (int i=1; i<=n; i++ )
       
		if(node[i].x==0.75 && node[i].y==0.75 && node[i].z==0.75){
			out.precision(15);
		    out<<"T(0.75,0.75,0.75) = "<< node[i].T<<'\n';
		}
       
	}
	
//	ini.close();
}

void solve (void)
{
  long int i, j, k, cur, tcol, it;
  double curerror, maxerror, *swap;
  int IS_CONVERGE, CONVERGE; // is equal to zero if all domain has converged otherwisw 1

  //Iteration loop (PJ)
  for( it=1; it<itr; it++ )
    {
      //Error Initialization
      curerror = maxerror = maxRes = 0.0;
      IS_CONVERGE = CONVERGE = 1;

      if (opt== 0) PJ_mpi(); //for PJ
      else SOR_mpi();

      //Error calculation
      for (k = k_min[myrank]; k <= k_max[myrank]; k++)
	for (j = 1; j <= N1; j++)
	  for( i = 1; i <= N1; i++ )
	    {
	      cur = INDEX(i, j, k);
	      if( (oldT[cur]-0.0) > 1e-6 )curerror = fabs( (Tn[cur]/oldT[cur]) - 1.0 );
	      if( curerror > maxerror ) maxerror = curerror;
	    }

      if( (maxerror < eps) && (it>1) ) IS_CONVERGE = 0;
      if( maxRes < maxerror ) maxRes = maxerror;

      MPI_Allreduce(&IS_CONVERGE, &CONVERGE, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      if(CONVERGE == 0) break;

      //Update
      swap = oldT;
      oldT = Tn;
      Tn = swap;

    }//Iteration loop
  noit = it;
}


void PJ_mpi(void)
{

  int i, j, k;
  double h;
  int requests;
  MPI_Request request[4];
  MPI_Status status[4];

  h = 1.0 / (N1 + 1);// lattice spacing

  // update ghost layers using non-blocking send/receive
  requests = 0;
  if (left_proc[myrank] >= 0 && left_proc[myrank] < n_proc) {
    MPI_Irecv(oldT + INDEX(0, 0, k_min[myrank] - 1), (N1+2)*(N1+2), MPI_DOUBLE, left_proc[myrank], 0, MPI_COMM_WORLD, request + requests++);
    MPI_Isend(oldT + INDEX(0, 0, k_min[myrank]), (N1+2)*(N1+2), MPI_DOUBLE, left_proc[myrank], 1, MPI_COMM_WORLD, request + requests++);
  }
  if (right_proc[myrank] >= 0 && right_proc[myrank] < n_proc) {
    MPI_Irecv(oldT + INDEX(0, 0, k_max[myrank] + 1), (N1+2)*(N1+2), MPI_DOUBLE, right_proc[myrank], 1, MPI_COMM_WORLD, request + requests++);
    MPI_Isend(oldT + INDEX(0, 0, k_max[myrank]), (N1+2)*(N1+2), MPI_DOUBLE, right_proc[myrank], 0, MPI_COMM_WORLD, request + requests++);
  }

  // Jacobi update for internal vertices in my domain
  for (k = k_min[myrank] + 1; k <= k_max[myrank] - 1; k++)
    for (j = 1; j <= N1; j++)
      for (i = 1; i <= N1; i++)
	Tn[INDEX(i,j,k)] =  ( oldT[INDEX(i-1,j,k)] + oldT[INDEX(i+1,j,k)] + oldT[INDEX(i,j-1,k)] +
			     oldT[INDEX(i,j+1,k)] +  oldT[INDEX(i,j,k-1)] + oldT[INDEX(i,j,k+1)] + q/(I*I*cnd) )/6.0;

  // wait for all non-blocking communications to complete
  MPI_Waitall(requests, request, status);

  // Jacobi update for boundary vertices in my domain
  k = k_min[myrank];
  for (j = 1; j <= N1; j++)
    for (i = 1; i <= N1; i++)
      Tn[INDEX(i,j,k)] =  ( oldT[INDEX(i-1,j,k)] + oldT[INDEX(i+1,j,k)] + oldT[INDEX(i,j-1,k)] +
			   oldT[INDEX(i,j+1,k)] +  oldT[INDEX(i,j,k-1)] + oldT[INDEX(i,j,k+1)] + q/(I*I*cnd) )/6.0;

  k = k_max[myrank];
  if (k != k_min[myrank]) {
    for (j = 1; j <= N1; j++)
      for (i = 1; i <= N1; i++)
	Tn[INDEX(i,j,k)] =  ( oldT[INDEX(i-1,j,k)] + oldT[INDEX(i+1,j,k)] + oldT[INDEX(i,j-1,k)] +
			     oldT[INDEX(i,j+1,k)] +  oldT[INDEX(i,j,k-1)] + oldT[INDEX(i,j,k+1)] + q/(I*I*cnd) )/6.0;
  } //End of if

}


void SOR_mpi(void)
{

	 int i, j, k;
	 long cur;
	 double h;
	 int requests;
	 MPI_Request request[4];
	 MPI_Status status[4];


	 // update ghost layers using non-blocking send/receive
	 requests = 0;
	 if (left_proc[myrank] >= 0 && left_proc[myrank] < n_proc) {
		 MPI_Irecv(oldT + INDEX(0, 0, k_min[myrank] - 1), (N1+2)*(N1+2), MPI_DOUBLE, left_proc[myrank], 0, MPI_COMM_WORLD, request + requests++);
		 MPI_Isend(oldT + INDEX(0, 0, k_min[myrank]), (N1+2)*(N1+2), MPI_DOUBLE, left_proc[myrank], 1, MPI_COMM_WORLD, request + requests++);
	 }
	 if (right_proc[myrank] >= 0 && right_proc[myrank] < n_proc) {
		 MPI_Irecv(oldT + INDEX(0, 0, k_max[myrank] + 1), (N1+2)*(N1+2), MPI_DOUBLE, right_proc[myrank], 1, MPI_COMM_WORLD, request + requests++);
		 MPI_Isend(oldT + INDEX(0, 0, k_max[myrank]), (N1+2)*(N1+2), MPI_DOUBLE, right_proc[myrank], 0, MPI_COMM_WORLD, request + requests++);
	 }

	 // updating Red index from old Black for internal vertices in my domain
	 for (k = k_min[myrank] + 1; k <= k_max[myrank] - 1; k++)
	  for (j = 1; j <= N1; j++)
	   for (i = 1; i <= N1; i++)
	   {
		   cur = INDEX(i,j,k);
		   if( cur%2 == 0 ) {
			   Tn[INDEX(i,j,k)] =  ( oldT[INDEX(i-1,j,k)] + oldT[INDEX(i+1,j,k)] + oldT[INDEX(i,j-1,k)] +
			   oldT[INDEX(i,j+1,k)] +  oldT[INDEX(i,j,k-1)] + oldT[INDEX(i,j,k+1)] + q/(I*I*cnd) )/6.0;

			   Tn[cur] = oldT[cur] + omega*(Tn[cur] - oldT[cur]); //Relaxing the temperature values for RED index
		   }
	   }

	 // wait for all non-blocking communications to complete
	 MPI_Waitall(requests, request, status);

	 // updating Red index from old Black for boundary vertices in my domain
	 k = k_min[myrank];
	 for (j = 1; j <= N1; j++)
	  for (i = 1; i <= N1; i++)
	  {
		  cur = INDEX(i,j,k);
		  if( cur%2 == 0 ) {
			  Tn[INDEX(i,j,k)] =  ( oldT[INDEX(i-1,j,k)] + oldT[INDEX(i+1,j,k)] + oldT[INDEX(i,j-1,k)] +
			  oldT[INDEX(i,j+1,k)] +  oldT[INDEX(i,j,k-1)] + oldT[INDEX(i,j,k+1)] + q/(I*I*cnd) )/6.0;

			  Tn[cur] = oldT[cur] + omega*(Tn[cur] - oldT[cur]); //Relaxing the temperature values for RED index
		  }
	  }

	 k = k_max[myrank];
	 if (k != k_min[myrank]) {
	  for (j = 1; j <= N1; j++)
	   for (i = 1; i <= N1; i++)
		{
			  cur = INDEX(i,j,k);
			  if( cur%2 == 0 ) {
				  Tn[INDEX(i,j,k)] =  ( oldT[INDEX(i-1,j,k)] + oldT[INDEX(i+1,j,k)] + oldT[INDEX(i,j-1,k)] +
				  oldT[INDEX(i,j+1,k)] +  oldT[INDEX(i,j,k-1)] + oldT[INDEX(i,j,k+1)] + q/(I*I*cnd) )/6.0;

				  Tn[cur] = oldT[cur] + omega*(Tn[cur] - oldT[cur]); //Relaxing the temperature values for RED index
			  }
	  }//end of i for
	 } //End of if 'k'

	// For Black indexes

	// update ghost layers using non-blocking send/receive
	 requests = 0;
	 if (left_proc[myrank] >= 0 && left_proc[myrank] < n_proc) {
		 MPI_Irecv(Tn + INDEX(0, 0, k_min[myrank] - 1), (N1+2)*(N1+2), MPI_DOUBLE, left_proc[myrank], 0, MPI_COMM_WORLD, request + requests++);
		 MPI_Isend(Tn + INDEX(0, 0, k_min[myrank]), (N1+2)*(N1+2), MPI_DOUBLE, left_proc[myrank], 1, MPI_COMM_WORLD, request + requests++);
	 }
	 if (right_proc[myrank] >= 0 && right_proc[myrank] < n_proc) {
		 MPI_Irecv(Tn + INDEX(0, 0, k_max[myrank] + 1), (N1+2)*(N1+2), MPI_DOUBLE, right_proc[myrank], 1, MPI_COMM_WORLD, request + requests++);
		 MPI_Isend(Tn + INDEX(0, 0, k_max[myrank]), (N1+2)*(N1+2), MPI_DOUBLE, right_proc[myrank], 0, MPI_COMM_WORLD, request + requests++);
	 }

	 // updating Black index from new Red index for internal vertices in my domain
	 for (k = k_min[myrank] + 1; k <= k_max[myrank] - 1; k++)
	  for (j = 1; j <= N1; j++)
	   for (i = 1; i <= N1; i++)
	   {
		   cur = INDEX(i,j,k);
		   if( cur%2 == 1 ) {
			   Tn[INDEX(i,j,k)] =  (Tn[INDEX(i-1,j,k)] + Tn[INDEX(i+1,j,k)] + Tn[INDEX(i,j-1,k)] +
			   Tn[INDEX(i,j+1,k)] +  Tn[INDEX(i,j,k-1)] + Tn[INDEX(i,j,k+1)] + q/(I*I*cnd) )/6.0;

			   Tn[cur] = oldT[cur] + omega*(Tn[cur] - oldT[cur]); //Relaxing the temperature values for RED index
		   }
	   }

	 // wait for all non-blocking communications to complete
	 MPI_Waitall(requests, request, status);

	 // updating Red index from old Black for boundary vertices in my domain
	 k = k_min[myrank];
	 for (j = 1; j <= N1; j++)
	  for (i = 1; i <= N1; i++)
	  {
		  cur = INDEX(i,j,k);
		  if( cur%2 == 1 ) {
			  Tn[INDEX(i,j,k)] =  ( Tn[INDEX(i-1,j,k)] + Tn[INDEX(i+1,j,k)] + Tn[INDEX(i,j-1,k)] +
			  Tn[INDEX(i,j+1,k)] +  Tn[INDEX(i,j,k-1)] + Tn[INDEX(i,j,k+1)] + q/(I*I*cnd) )/6.0;

			  Tn[cur] = oldT[cur] + omega*(Tn[cur] - oldT[cur]); //Relaxing the temperature values for RED index
		  }
	  }

	 k = k_max[myrank];
	 if (k != k_min[myrank]) {
	  for (j = 1; j <= N1; j++)
	   for (i = 1; i <= N1; i++)
		{
			  cur = INDEX(i,j,k);
			  if( cur%2 == 1 ) {
				  Tn[INDEX(i,j,k)] =  ( Tn[INDEX(i-1,j,k)] + Tn[INDEX(i+1,j,k)] + Tn[INDEX(i,j-1,k)] +
				  Tn[INDEX(i,j+1,k)] +  Tn[INDEX(i,j,k-1)] + Tn[INDEX(i,j,k+1)] + q/(I*I*cnd) )/6.0;

				  Tn[cur] = oldT[cur] + omega*(Tn[cur] - oldT[cur]); //Relaxing the temperature values for RED index
			  }
	  }//end of i for
	 } //End of if 'k'



	 //synchronization point to make sure that all processes have finished calculation for Black index
	 //MPI_Barrier(MPI_COMM_WORLD);

}


void print_output(void)
{

  long int i, j, k, cur;
  double x, y, z;

  out<<"Number of iteration used  = "<<noit<<endl;
  out<<"Maximum iteration residual = "<<maxRes<<endl;

     // 0.25
    i = (I/4);
    j = (I/4);
    k = (I/4);
    z = (double) (k)/I ;
    y = (double) (j)/I;
    x = (double) (i)/I;

    cur = INDEX(i, j, k);
    out.precision(15);
    out<<"T(0.25, 0.25, 0.25) = "<<Tn[cur]<<endl;

    // 0.5
    i = (I/2);
    j = (I/2);
    k = (I/2);
    z = (double) (k)/I ;
    y = (double) (j)/I;
    x = (double) (i)/I;

    cur = INDEX(i, j, k);
    out.precision(15);
    out<<"T(0.5, 0.5, 0.5) = "<<Tn[cur]<<endl;

    // 0.75
    i = (3*I/4);
    j = (3*I/4);
    k = (3*I/4);
    z = (double) (k)/I ;
    y = (double) (j)/I;
    x = (double) (i)/I;

    cur = INDEX(i, j, k);
    out.precision(15);
    out<<"T(0.75, 0.75, 0.75) = "<<Tn[cur]<<endl;
 
}




void SetEqn(Element *node){
	long int i, j, k, n;
		
	//double **A, *B;   // Part of old code
	
	
	
	// Total number of grid points
	n = nx*ny*nz;
	
	C = new double[n+1];
	W = new double[n+1];
	N = new double[n+1];
	E = new double[n+1];
	S = new double[n+1];
	T = new double[n+1];
	B = new double[n+1];
	Q = new double[n+1];
	
	for( i=0; i<=n; i++) C[i] = W[i] = N[i] = E[i] = S[i] = T[i] = B[i] = Q[i]= 0.0;

	//Creation of band matrices
	for( int l = 1; l<=ny; l++)
	{
		for( j = 1; j<=nz; j++)
		{
			for( i = 1; i<=nx; i++)
			{
				k = i + (j-1)*nx + (l-1)*nx*ny;

				//Bdry Nodes
				if( i==1 || j==1 || l==1 || i==nx || j==nx || l==nx ) {
					C[k] = 1.0;
					Q[k] = 0.0 ;
					continue;
				}// End of if loop


				W[k] = E[k] = 1.0;
				N[k] = S[k] = 1.0;
				B[k] = T[k] = 1.0;

				C[k] = -6.0;

				Q[k] = -q/(I*I*cnd);
			}//end of i loop
		}// end og j loop
	}//end of k loop


	//Call to solver routine
 	SolveEqn();

	 for( i=0; i<(n+1); i++ )
  		node[i].T = Q[i];

 	delete[] C;
  	delete[] W;
  	delete[] E;
  	delete[] T;
  	delete[] B;
  	delete[] N;
  	delete[] S;
  	delete[] Q;
}
	
void SolveEqn(){
		
  double sum=0.0;
  double *Td, *Tn;        // Dummy array and a new array to store the intermediate solution
  int i,k;
  int n = nx*ny*nz;
  double error;
  
		Td = new double[n+1];
		Tn = new double[n+1];
		
		for (i=0; i<(n+1); i++) { Td[i] = Tn[i]= 0.0; }

		//ofstream red;
		//red.open("red_o.dat",ios::out);

		switch(opt){
	        case 0:
			
			for(k=1; k<=itr ; k++){
			  /*#pragma omp parallel private(i, TID, east, west, south, north, bot, top, sum) shared (nx, E, W, N, S, Td, Tn, T, B, C, TH, Q) num_threads(TH)
{  				
  TID = omp_get_thread_num();
  //  cout<<" I am thread "<< TID <<" of "<< TH <<" threads"<<endl;  	
  #pragma omp for schedule(static) */
 for (i=1; i<= n ; i++){
    
				west = i - 1;
				east = i + 1;

				south = i - nx;
				north = i + nx;

				bot = i - nx*nx;
				top = i + nx*nx;

				if( west < 0 || west >nx*nx*nx ) west =0;
				if( east < 0 || east >nx*nx*nx ) east =0;
				if( south < 0 || south >nx*nx*nx ) south =0;
				if( north < 0 || north >nx*nx*nx ) north =0;
				if( bot < 0 || bot >nx*nx*nx ) bot =0;
				if( top < 0 || top >nx*nx*nx ) top =0;
		        	
	   	 		sum = 0.0;
				sum = sum + E[i]*Td[east] + W[i]*Td[west] + N[i]*Td[north] + S[i]*Td[south] + T[i]*Td[top] + B[i]*Td[bot];
        			
				Tn[i] = (Q[i]-sum)/C[i];
 } // for loop ended
 //} // parallel ended
	   	  	  	error = err(Td,Tn);
	   	  	  
	   	  	  	if(error<=eps){
	   	  	  		out<<"Number of iterations required for the Point Jacobi = "<< k<<endl;
	   	  	  		break;
	   	  	  	}  
	   	  	  	
	   	  	  	if(error>= eps && k>=itr){
      				out<<"Convergence cannot be reached after  "<< k <<" iterations"<<endl;
      				out<<"Maximum error = "<<error<<endl;
      			}
				/*   			
				#pragma omp parallel private(i) shared(Td, Tn) num_threads(TH)
				{
				#pragma omp for schedule (static) */
				for(i=1; i<=n; i++)Td[i] = Tn[i];
				//} 
      		} // Iterations completed
		
			for(i=1; i<=n; i++)Q[i] = Td[i];
			break; // Point Jacobi Solver Finished 
	
	
		case 1:
			for(k=1; k<=itr ; k++){
		        for(i=1; i<= n ; i++){
		        	
		        	west = i - 1;
					east = i + 1;

					south = i - nx;
					north = i + nx;

					bot = i - nx*nx;
					top = i + nx*nx;
		        	
      		 		sum=0.0;
      		 		sum = sum + E[i]*Td[east] + W[i]*Tn[west] + N[i]*Td[north] + S[i]*Tn[south] + T[i]*Td[top] + B[i]*Tn[bot];
        
         			Tn[i] = (Q[i]-sum)/C[i];
           		}
      	
      			error = err(Td,Tn);
	   	  	  	
	   	  	  	if(error<=eps){
	   	  	  		out<<"Number of iterations required for the Gauss-Seidel = "<< k<<endl;
	   	  	  		break;
	   	  	  	}  
      	      	
      	      	if(error>= eps && k>=itr){
      				out<<"Convergence cannot be reached after  "<< k <<" iterations"<<endl;
      				out<<"Maximum error = "<<error<<endl;
      			}
				for(i=1; i<=n; i++)Td[i] = Tn[i];
      		} // Iterations completed
		
		
			for(i=1; i<=n; i++)Q[i] = Td[i];
			break; // Gauss-Seidel Solver Finished 
	
	
		case 2:
		
		  for(k=1; k<=itr ; k++){
		    /*		    
#pragma omp parallel private(i, TID, east, west, south, north, bot, top, sum) shared (nx, E, W, N, S, Td, Tn, T, B, C, TH, Q,omega,n) num_threads(TH)
		    {  				
		    #pragma omp for schedule(static) */
		      for(i=1; i<= n ; i=i+2){    // Updating Red index from Old Black index
			//	sum=0.0;
			
			west = i - 1;
			east = i + 1;
			      
			south = i - nx;
			north = i + nx;
			      
			bot = i - nx*nx;
			top = i + nx*nx;
		        	
			if( west < 0 || west >nx*nx*nx ) west =0;
			if( east < 0 || east >nx*nx*nx ) east =0;
			if( south < 0 || south >nx*nx*nx ) south =0;
			if( north < 0 || north >nx*nx*nx ) north =0;
			if( bot < 0 || bot >nx*nx*nx ) bot =0;
			if( top < 0 || top >nx*nx*nx ) top =0;
			
			sum=0.0;
			sum = sum + E[i]*Td[east] + W[i]*Td[west] + N[i]*Td[north] + S[i]*Td[south] + T[i]*Td[top] + B[i]*Td[bot];
			      
			Tn[i] = (Q[i]-sum)/C[i];
			//			red<<" I am starting red for " << k << " loop" << endl;
			//			red<<i<<endl;
			
			Tn[i]= Td[i] + omega*(Tn[i]-Td[i]);
		      } // End red iterations

		      //#pragma omp for schedule(static)
		      for(i=2; i<= n ; i=i+2){    // Updating black index from Old red index                                                                                
			//	sum=0.0;

			west = i - 1;
			east = i + 1;

			south = i - nx;
			north = i + nx;

			bot = i - nx*nx;
			top = i + nx*nx;

			if( west < 0 || west >nx*nx*nx ) west =0;
			if( east < 0 || east >nx*nx*nx ) east =0;
			if( south < 0 || south >nx*nx*nx ) south =0;
			if( north < 0 || north >nx*nx*nx ) north =0;
			if( bot < 0 || bot >nx*nx*nx ) bot =0;
			if( top < 0 || top >nx*nx*nx ) top =0;
			
			sum=0.0;
			sum = sum + E[i]*Tn[east] + W[i]*Tn[west] + N[i]*Tn[north] + S[i]*Tn[south] + T[i]*Tn[top] + B[i]*Tn[bot];

			Tn[i] = (Q[i]-sum)/C[i];

			Tn[i]= Td[i] + omega*(Tn[i]-Td[i]);
		      } // End black iterations
		      //} // End Parallel

		    error = err(Td,Tn);
	   	  	 
		    if(error<=eps){
		      out<<"Number of iterations required for the SOR = "<< k<<endl;
		      break;
		    }  
      	
		    if(error>= eps && k>=itr){
		      out<<"Convergence cannot be reached after  "<< k <<" iterations"<<endl;
		      out<<"Maximum error = "<<error<<endl;
		    }


		    /*#pragma omp parallel private(i) shared(Td, Tn) num_threads(TH)
		    {
		    #pragma omp for schedule (static) */
		      for(i=1; i<=n; i++)Td[i] = Tn[i];
		      //}

		  } // Iterations completed

//		  red.close();
		
		  for(i=1; i<=n; i++)Q[i] = Td[i];
		  break; // SOR Solver Finished 
		}	
		delete[] Td;
		delete[] Tn;
	
	}
	
	
inline	double err(double *Ta, double *Tb){
		
		double *e, maxerr=-32000;
		e = new double[n+1];
		for(int i=1;i<=n;i++)e[i]=fabs( (Tb[i]-Ta[i])/Ta[i] );
		
		for(int i=1;i<=n;i++){
			
		 	if (e[i]>maxerr)
	 			{
	   	 		maxerr=e[i];
	 			}
			}
	delete[] e;
	return maxerr;
	}
	
		
		
		

	
	
