// This program is written by Kamalesh S. Bhambare 

#include<iostream>
#include<math.h>
#include<cstdlib>
#include<time.h>
#include<fstream>
#include<omp.h>
using namespace std;

// Variable declaration
double R, r1, r2;
double Vsp, Vsh;  // Volumes of the spheres inner and outer 

double dN;        
double ds;

int s;           // Seed for random number generator 
long int N;           // N- Number of iterations and w is number within the sphere
int TH;

time_t date;
double sum=0.0,sqrsum=0.0;
double q_avg, std_dev;
double sum2=0.0,sqrsum2=0.0;
double q_avg2, std_dev2;
double exec_time;


// Function declartion
double volume(int,double);   
double volume_sh(int, double, double);
void check_in();
void output();
void echoin();
void read();

ofstream out;
  
int main(){

  time_t starttime, endtime;
  starttime = time (NULL);

  out.open("output.dat",ios::out);
  //  out<<"Hello";  
  // Reading the inputs from file input.dat
  read();
  check_in();
  echoin();

  //  srand(s);   //Seed for psuedo random generator
  //  srand(time(0));  // Random number generator seeding from system time

  Vsp = volume(N,R);   // Calculates the volume of the sphere
  Vsh = volume_sh(N,r1,r2);  // Calculates the volume of shell 
  output();

  endtime = time (NULL);
  exec_time = difftime(endtime, starttime);

  out<<"Execution time = "<<exec_time << " seconds \n";
  out.close();
}



// Monte Carlo code to calculate the volume of the sphere
double volume(int n1, double rad){
  
  ofstream plot1;
  plot1.open("Plot1.dat",ios::app); // Output file for writting the standard deviation and mean volume for error bar plot
  
  double w = 0.0;
  double qi, vol;
  int  TID;
  long int i;

#pragma omp parallel private(TID) num_threads(TH) reduction( +: w) 
 {
    TID = omp_get_thread_num();
    cout<<" I am thread # " << TID<< " of " << TH <<endl;

#pragma omp for
    for (i = 0 ; i <= n1 ; i++) {
      double x = rad*rand()/((double)RAND_MAX);   // This generates random number in the range of 0-1.
      double y = rad*rand()/((double)RAND_MAX);
      double z = rad*rand()/((double)RAND_MAX);

      if (sqrt(x*x + y*y + z*z)< rad) {
	w = w + 1.0;
	//   qi = (w)/n1;
	//  sum = sum + qi;
	//sqrsum = sqrsum + pow(qi,2);
      }

    } // End parallel for loop
 } // End parallel

  // Calculating the standard deviation
  sqrsum = w;
  sqrsum = sqrsum/n1;
  q_avg = w/n1;
  std_dev = sqrt(( sqrsum - pow(q_avg,2) )/n1)*8*pow(rad,3);

  vol = w/n1*8*pow(rad,3);
  
  plot1<<n1<<'\t'<<vol<<'\t'<<std_dev<<'\t'<< 4.0/3.0*3.14159*pow(rad,3)<<'\n';
  plot1.close();

  return vol;
}


double volume_sh(int n2, double r_1, double r_2){

  ofstream plot2;
  plot2.open("Plot2.dat",ios::app);

  double w2 = 0.0;
  double vol2, qi2;
  int TID;
  long int i;

#pragma omp parallel private(TID) num_threads(TH) reduction( +: w2) 
{

    TID = omp_get_thread_num();
    cout<<" I am thread # " << TID<< " of " << TH <<endl;

#pragma omp for  
for (i = 0 ; i < n2 ; i++){
  
    double x2 = r_2*rand()/((double)RAND_MAX);
    double y2 = r_2*rand()/((double)RAND_MAX);
    double z2 = r_2*rand()/((double)RAND_MAX);
    if (sqrt(x2*x2 + y2*y2 + z2*z2)< r_2 && sqrt(x2*x2 + y2*y2 + z2*z2)>r_1){

      w2 = w2 + 1.0;
      //   qi2 = (w2)/n2;
      //   sum2= sum2 + qi2;
      // sqrsum2= sqrsum2 + pow(qi2,2);
    }
    
} // end parallel for

} // end parallel
  q_avg2 = w2/n2;
  sqrsum2= w2/n2;
  std_dev2 = sqrt((sqrsum2-pow(q_avg2,2))/n2)*8*pow(r2,3); 

  vol2 = w2/n2*8*pow(r2,3);
  
  plot2<<n2<<'\t'<<vol2<<'\t'<<std_dev2<<'\t'<< 4.0/3.0*3.14159*(pow(r2,3)-pow(r1,3))<<'\n';
  plot2.close();

  return vol2;

}

    
void echoin(){
  

  time(&date);
  out<<"Code Name: "<< "main" <<'\n';
  out<<"Version Number: "<< "1.0" <<'\n';
  out<<"Author Name: "<< "Kamalesh Bhambare " <<'\n';
  out<<"Date and time of execution: " << ctime(&date) <<'\n';
  out<<"Radius of sphere          = "<< R <<'\n' ;
  out<<"Inside radius of sphere   = "<< r1 <<'\n';
  out<<"Outside radius of sphere  = "<< r2 <<'\n' ;
  out<<"Number of histories       = "<< N <<'\n' ;
  out<<"Seed for random number    = "<< s <<'\n' ;

}
 
 
void output(){
  
  
  out<<"Volume of sphere          = "<< Vsp <<'\n';
  out<<"Standard deviation        = "<< std_dev<<'\n';
  out<<"Volume of spherical shell = "<< Vsh <<'\n';
  out<<"Standard deviation        = "<< std_dev2<<'\n';

}

void read(){
 ifstream inp;
  inp.open("input.dat",ios::in);
  if(!inp){
    cout<<"Input file doesnot exist"<<'\n';
    exit(1);
  }

  inp>>R>>r1>>r2>>dN>>ds>>TH;
  N = (int)dN;               // Type conversion from double to integer for reading any kind of input from user
  s = (int)ds;
  inp.close();
}


void check_in(void){

  //checking for positive radius value (Input check)
  if( R <= 0) {
    cout << "Radius of sphere should be greater than zero\nExiting....\n" ;
    exit (0);
  }

  if( r1 <= 0 || r2 <= 0) {
    cout<< "Radius of spherical shell should be greater than zero\nExiting....\n" ;
    exit (0);
  }

  if( r2 <= r1 ) {
    cout<<"Outer radius of spherical shell should be greater than inner radius\nExiting....\n";
    exit (0);
  }

  if( N <= 0) {
    cout<<"Number of histories should be greater than zero\nExiting....\n" ;
    exit (0);
  }

  if( s <= 0) {
    cout<<"seed should be greater than zero\nExiting....\n" ;
    exit (0);
  }

  if( TH <= 0) {
    cout<<"Number of threads should be greater than zero\nExiting....\n" ;
    exit (0);
  }

   out<< "Input check...successful.\n\n" ;
}

