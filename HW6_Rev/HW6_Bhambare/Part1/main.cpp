// Code written by Kamalesh Bhambare to solve heat equation with source in cube with banded matrices and a structure to store temperture 
// and location of grid point
#include<iostream>
#include<fstream>
#include<math.h>
#include<time.h>
#include<omp.h>

using namespace std;

// Defining the structure for storing variables at each grid point
struct Element
{	
	double T;          
	double x;
	double y;
	double z;
};

void check_in(void);    // Checks input for negative and non zero values as well as for input file existance
void echoin(void);       // Echos input
void grid(Element *);   // Calculate x and y for grid point
void init(Element *);   // Initialise the grid variables
void SetEqn(Element *); // Setting the coefficients of the matrix
void SolveEqn();
inline double err(double *, double *);

void print_node (Element *);  // For initial checking
void plot_node(Element *);
double *C, *W, *N, *E, *S, *T, *B, *Q; 
int nx, ny, nz, n, I;
double nxd;
double dx,dy,dz;
double xlength, ylength, zlength;
double cnd, q, Tb;
int info;
time_t date;
double exec_time;
int east, west, top, bot, north, south;
bool flag;
double eps,omega;              // Parameters used for iterative methods
int itr,opt,TH, TID; 

ofstream out;
ofstream outfile_seteq;

int main(){
	clock_t start = clock();
	check_in();
	echoin();
	
	Element *node;
	node = new Element[n+1];

	init(node); 
	grid(node);  // Generating grid with stored information of x,y,z of each node
		
	SetEqn(node); // Setting up the coefficients and solving for the temperature
	exec_time = ((double)clock() - start)/CLOCKS_PER_SEC;	
	plot_node(node); // Printing results at the specified cells
	
//	print_node(node); // Printing whole grid data 
		

	out<<"Execution time = "<<exec_time << " seconds \n";
	out.close();
	delete[] node;
	
} // main ends here

void check_in(){
	
 	ifstream inp;
  	inp.open("input.dat",ios::in);
  	
  	if(!inp){
    cout<<"Input file doesnot exist"<<'\n';
    exit(1);
  }	

  	inp >> nxd >> q >> cnd >> eps >> itr >> opt >> TH;
  	inp.close();
  	
  	I = (int) nxd;
  	nx = I + 1 ;
  
 	 if( I <= 1) {
		cout<<"Number of computational cells should be greater than one\nExiting....\n";
		exit (0);
	}	
    
  	if( q <= 0) {
		cout<<" Negative or zero source\nExiting....\n";
		exit (0);
	}

	if( cnd <= 0) {
			cout<<"Thermal condictivity can't be negative or zero\nExiting....\n";
			exit (0);
	}

	if(I%4 != 0) {
		cout<<"Number of cells should be multiple of 4 to place grid at 0.25, 0.5 and 0.75 location\n";
		cout<<"Still running the code for "<< nx << " cells \n";
	}
	
	if( itr <= 0) {
			cout<<"Number of iterations can't be negative or zero\nExiting....\n";
			exit (0);
	}
	
	if( opt >2 || opt< 0) {
			cout<<"Not a valid input for iterative method selection, Use 0 for Point Jacobi, \n";
			cout<<"1 for Gauss-Seidel, 2 for SOR ....\n";
			exit (0);
	}
	

	if(TH <= 0){
	  cout<<" Number of threads must be an integer and greater than zero \n ";
	}
	

	flag = true;
	ny = nx;
  	nz = nx; 	
  	n = nx*ny*nz;
  	xlength = 1; 
  	ylength = 1;
  	zlength = 1; 	
  	Tb = 273;	
  	omega = 2/(1+sin(3.14/I));	
}

void echoin()
{
	out.open("output.dat",ios::out);
  	time(&date);	
  	out<<"Code Name: "<< "main" <<'\n';
  	out<<"Version Number: "<< "1.0" <<'\n';
  	out<<"Author Name: "<< "Kamalesh Bhambare " <<'\n';
  	out<<"Date and time of execution: " << ctime(&date) <<'\n';
  	
  	if(flag==true)out<<"Input check...successful.\n\n";
	
	out<<"Number of cells per dimension = "<< I <<'\n';
	out<<"Uniform distributed source strength = "<< q <<'\n';
	out<<"Thermal conductivity = "<< cnd <<'\n';
	out<<"Relative convergence criterion = "<< eps <<'\n';
	out<<"Number of iterations used during solution = "<< itr << '\n';
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

	for ( i=0; i<=n; i++ )
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
		
		for (i=0; i<(n+1); i++) Td[i] = Tn[i]= 0.0;
		
		switch(opt){
	        case 0:
			
			for(k=1; k<=itr ; k++){
#pragma omp parallel private(i, TID, east, west, south, north, bot, top, sum) shared (nx, E, W, N, S, Td, Tn, T, B, C, TH, Q) num_threads(TH)
{  				
  TID = omp_get_thread_num();

#pragma omp for schedule(static)
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
} // parallel ended
	   	  	  	error = err(Td,Tn);
	   	  	  
	   	  	  	if(error<=eps){
	   	  	  		out<<"Number of iterations required for the Point Jacobi = "<< k<<endl;
	   	  	  		break;
	   	  	  	}  
	   	  	  	
	   	  	  	if(error>= eps && k>=itr){
      				out<<"Convergence cannot be reached after  "<< k <<" iterations"<<endl;
      				out<<"Maximum error = "<<error<<endl;
      			}
      			
				#pragma omp parallel private(i) shared(Td, Tn) num_threads(TH)
				{
				#pragma omp for schedule (static)
				for(i=1; i<=n; i++)Td[i] = Tn[i];
				}
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
		        for(i=1; i<= n ; i++){
      		 		sum=0.0;
      		 		
      		 		    	west = i - 1;
					east = i + 1;

					south = i - nx;
					north = i + nx;

					bot = i - nx*nx;
					top = i + nx*nx;
		        	
      		 		sum=0.0;
      		 		sum = sum + E[i]*Td[east] + W[i]*Tn[west] + N[i]*Td[north] + S[i]*Tn[south] + T[i]*Td[top] + B[i]*Tn[bot];
   
         			Tn[i] = (Q[i]-sum)/C[i];
 		 		
         			Tn[i]= Td[i] + omega*(Tn[i]-Td[i]);
           		}
      	
      			error = err(Td,Tn);
	   	  	 
	   	  	  	if(error<=eps){
	   	  	  		out<<"Number of iterations required for the SOR = "<< k<<endl;
	   	  	  		break;
	   	  	  	}  
      	
      			if(error>= eps && k>=itr){
      				out<<"Convergence cannot be reached after  "<< k <<" iterations"<<endl;
      				out<<"Maximum error = "<<error<<endl;
      			}
				for(i=1; i<=n; i++)Td[i] = Tn[i];
      		} // Iterations completed
		
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
	
		
		
		

	
	
