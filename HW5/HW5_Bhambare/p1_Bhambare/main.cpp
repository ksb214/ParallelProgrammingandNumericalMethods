// Program written by Kamalesh S. Bhambare
#include<iostream>
#include<fstream>
#include<omp.h>
using namespace std;

int main(){
  int TH, TID;

  ifstream in;
  in.open("input.dat",ios::in);
  in>>TH;  
  in.close();

  ofstream out;
  out.open("output.dat",ios::out);
  
#pragma omp parallel \
num_threads(TH) \
private(TID)
  {
    TID = omp_get_thread_num();
    out<<" I am thread # "<<TID<<" of "<<TH<<" threads \n";

  } // End of parallel section

  out.close();
} // End of main program
