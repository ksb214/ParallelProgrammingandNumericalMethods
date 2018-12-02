
void gauss(double **a,double *b,long eqn);

/*
       To solve the linear equation by gauss elimination using row and column pivoting

a = coeffecient matrix
b = RHS and the solution will be obtained in the same 
eqn = no. of quations to be solved in the solver i.e. no. of rows in matrix a

Note: After calling of this function the a matrix will be change
*/

void gauss(double **a,double *b,long eqn)
{
  int k,row,r,s;
  int i,j;
  double sum ,c,*temp,btemp,large=-1e-15;
 
  temp = (double *) malloc ( sizeof(double)*(eqn+1) ) ;
  
  for ( k = 1 ; k<=(eqn-1) ; k++)
    { 
      for ( i = k+1 ; i<=eqn ; i++)
      {
	/* pivoting begins */
	if ( a[k][k]==0 )
	  { 
	    for( r = k+1; r<=eqn ; r++) 
	      {
		if( (a[r][k]>large) && (a[r][k]!=0) ) 
		  {
		    large = a[r][k];
		    row=r;
		  }
	      }
	    
	    btemp = b[k];
	    b[k] = b[row];
	    b[row] = btemp;
	    
	    for( s = 1 ; s<=eqn ; s++) 
	      {
		temp[s] = a[k][s];
		a[k][s] = a[row][s];
		a[row][s] = temp[s];
	      }
	  }
	
	/* end of pivoting */
	large=-1e-15;
	
	c = a[i][k]/a[k][k];
	
	for( j = k+1; j<=eqn ; j++ ) 
	  a[i][j] = a[i][j]- c*a[k][j] ;
	
	b[i] = b[i] - c * b[k] ;
	
      }
    }
  
  b[eqn]=b[eqn]/a[eqn][eqn];
  for ( i = eqn-1 ; i>0 ; i--)
    {
      k = i+1;
      sum = 0 ;
      for( j = k; j<=eqn; j++)
	sum += a[i][j]*b[j];
      
      b[i] = (b[i]- sum)/a[i][i];
    }
  free (temp);
}

