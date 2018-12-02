// Code written by Kamalesh Bhambare to solve heat equation with source in cube
#include<iostream>
#include<fstream>
#include<math.h>
#include<time.h>

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
void SolveEqn(double **A, double *B);
double err(double *, double *);

void print_node (Element *);  // For initial checking
void plot_node(Element *);

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
int itr,opt; 

ofstream out;
ofstream outfile_seteq;

//double *B;  // Defining B matrix
//double **A; // Defining A matrix to store


int main(){
	clock_t start = clock();
	check_in();
	echoin();
	
	Element *node;
	node = new Element[n+1];

	init(node); 
	grid(node);  // Generating grid with stored information of x,y,z of each node
		
	SetEqn(node); // Setting up the coefficients and solving for the temperature
	
	plot_node(node); // Printing results at the specified cells
	
	print_node(node); // Printing whole grid data 
		
	exec_time = ((double)clock() - start)/CLOCKS_PER_SEC;
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

  	inp >> nxd >> q >> cnd >> eps >> itr >> opt;
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
	if(opt==0)out<<"Iterative method used = Point Jacobi \n";
	if(opt==1)out<<"Iterative method used = Gauss-Seidel \n";
	if(opt==2){
			out<<"Iterative method used = SOR"<<endl;
			out<<"Optimal SOR weight = "<< omega <<" \n";
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
	
  //	ofstream ini;
	
	//	ini.open("Results.dat",ios::out);
	for (int i=1; i<=n; i++ )
	{
		if(node[i].x== 0.5 && node[i].y==0.5 && node[i].z==0.5){
			 out<<"T(0.5,0.5,0.5) = "<< node[i].T<<'\n';
		}
	}
		
	
	for (int i=1; i<=n; i++ )
	{
		if(node[i].x==0.25 && node[i].y==0.25 && node[i].z==0.25){
			 out<<"T(0.25,0.25,0.25) = "<< node[i].T<<'\n';
		}
	
	}	
		for (int i=1; i<=n; i++ )
	{
		if(node[i].x==0.75 && node[i].y==0.75 && node[i].z==0.75){
			 out<<"T(0.75,0.75,0.75) = "<< node[i].T<<'\n';
		}
	}
	
	
		//ini.close();
}

void SetEqn(Element *node){
	long int i, j, k, n;
	long int front, top, right,  left, bot, back;
	double **A, *B;

	// Total number of grid points
	n = nx*ny*nz;

	//dynamically allocation of Matrices
	A = new double * [n+1];
		for (i=0; i<(n+1); i++) A[i] = new double[n+1]; 
		
		B = new double[n+1];
		for (i=0; i<(n+1); i++) B[i] =  0.0;

  	//Initialization of the A, B, T matrices

  	for( i=0; i<(n+1); i++ )
  	{
		for( j=0; j<(n+1); j++)
		  A[i][j] = 0;
		  A[i][i] = 0;
	      B[i] = 0;
	 }

	//Creation of A and B matrices
	for( int l = 1; l<=ny; l++)
	{
		for( j = 1; j<=nz; j++)
		{
			for( i = 1; i<=nx; i++)
			{
				k = i + (j-1)*nx + (l-1)*nx*ny;

				//Bdry Nodes
				if( i==1 || j==1 || l==1 || i==nx || j==nx || l==nx ) {
					A[k][k] = 1.0;
					B[k] = 0.0 ;
					continue;
				}// End of if loop

				left = k - 1;
				right = k + 1;

				bot = k - nx;
				top = k + nx;

				back = k - nx*nx;
				front = k + nx*nx;

				A[k][left] = A[k][right] = 1.0;
				A[k][bot] = A[k][top] = 1.0;
				A[k][back] = A[k][front] = 1.0;

				A[k][k] = -6.0;

				B[k] = -q/(I*I*cnd);
			}//end of i loop
		}// end og j loop
	}//end of k loop

	//Call to solver routine
 	SolveEqn(A,B);

	for( i=0; i<(n+1); i++ )
  	node[i].T = B[i];

  	delete[] B;
  	delete[] A;
	}
	
	
	void SolveEqn(double **A, double *B){
		
		double sum=0.0;
		double *Td, *Tn;        // Dummy array and a new array to store the intermediate solution
		int i,j,k;
		int n = nx*ny*nz;
		double error;
			
		Td = new double[n+1];
		Tn = new double[n+1];
		
		for (i=0; i<(n+1); i++) Td[i] = Tn[i]= 0.0;
		
		
		if(opt==0){
			
			for(k=1; k<=itr ; k++){
				
		        for(i=1; i<= n ; i++){
      		 		sum = 0.0;
        			for(j=1; j<=n ;j++){
            		if(i!=j)sum = sum + A[i][j]*Td[j];
            		}
         
         			Tn[i] = (B[i]-sum)/A[i][i];
           		}
   
	   	  	  	error = err(Td,Tn);
	   	  	  
	   	  	  	if(error<=eps){
	   	  	  		out<<"Number of iterations required for the Point Jacobi = "<< k<<endl;
	   	  	  		out<<"Maximum error = "<<error<<endl;
					break;
	   	  	  	}  
	   	  	  	
	   	  	  	if(error>= eps && k>=itr){
      				out<<"Convergence cannot be reached after  "<< k <<" iterations"<<endl;
      				out<<"Maximum error = "<<error<<endl;
      			}
      			
				for(i=1; i<=n; i++)Td[i] = Tn[i];
			
      		} // Iterations completed
		
		for(i=1; i<=n; i++)B[i] = Td[i];
	} // Point Jacobi Solver Finished 
	
	
	if(opt==1){
			for(k=1; k<=itr ; k++){
		        for(i=1; i<= n ; i++){
      		 		sum=0.0;
        			for(j=1; j<=i-1 ;j++)sum = sum + A[i][j]*Tn[j];
        			for(j=i+1; j<=n ;j++)sum = sum + A[i][j]*Td[j];
         
         			Tn[i] = (B[i]-sum)/A[i][i];
           		}
      	
      			error = err(Td,Tn);
	   	  	  	
	   	  	  	if(error<=eps){
	   	  	  		out<<"Number of iterations required for the Gauss-Seidel = "<< k<<endl;
					out<<"Maximum error = "<<error<<endl;
	   	  	  		break;
	   	  	  	}  
      	      	
      	      	if(error>= eps && k>=itr){
      				out<<"Convergence cannot be reached after  "<< k <<" iterations"<<endl;
      				out<<"Maximum error = "<<error<<endl;
      			}
				for(i=1; i<=n; i++)Td[i] = Tn[i];
      		} // Iterations completed
		
		
		for(i=1; i<=n; i++)B[i] = Td[i];
	} // Gauss-Seidel Solver Finished 
	
	
	if(opt==2){
		
			for(k=1; k<=itr ; k++){
		        for(i=1; i<= n ; i++){
      		 		sum=0.0;
        			for(j=1; j<=i-1 ;j++)sum = sum + A[i][j]*Tn[j];
        			for(j=i+1; j<=n ;j++)sum = sum + A[i][j]*Td[j];
         
         			Tn[i] = (B[i]-sum)/A[i][i];
         			Tn[i]= Td[i] + omega*(Tn[i]-Td[i]);
           		}
      	
      			error = err(Td,Tn);
	   	  	 
	   	  	  	if(error<=eps){
	   	  	  		out<<"Number of iterations required for the SOR = "<< k<<endl;
					out<<"Maximum error = "<<error<<endl;
	   	  	  		break;
	   	  	  	}  
      	
      			if(error>= eps && k>=itr){
      				out<<"Convergence cannot be reached after  "<< k <<" iterations"<<endl;
      				out<<"Maximum error = "<<error<<endl;
      			}
				for(i=1; i<=n; i++)Td[i] = Tn[i];
      		} // Iterations completed
		
		for(i=1; i<=n; i++)B[i] = Td[i];
	} // SOR Solver Finished 
	
	delete[] Td;
  	delete[] Tn;
	
	}
	
	
	double err(double *Ta, double *Tb){
		int n1=nx*ny*nz;
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
	
		
		
		

	
	
