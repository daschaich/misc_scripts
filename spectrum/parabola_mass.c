/* for use with "fit.c".
* fit to form  b[0] + b[1]*(x-L/2)^2 + b[4]*(-1)**x*exp(-b[5]*x)
*             + b[2]*exp(-b[3]*x) + b[6]*(-1)**x*exp(-b[7]*x)
*             + "image"
*
* compile with cofit.c min.c linalg.c
*/
#include <stdio.h>
#include <math.h>
#include "minimize.h"

extern int npar;
int NT;

void f_init(FILE *fp)
{
  fprintf(stderr,"Kogut-Susskind parabola plus mass fitting\n");
  fprintf(stderr,"Enter time size of lattice:\n");
  fscanf(fp,"%d",&NT);
  npar = 8;
}

double f(double x, double *b)
/* b[0] = constant in parabola
*  b[1] = curvature of parabola
*  b[2] = coefficient 2 
*  b[3] = mass 2 (natural parity)
*  b[5] = coefficient 3 (alternating term)
*  b[6] = mass 3
*  b[7] = coefficient 4 (alternating term)
*  b[6] = mass 4
*/
{
	double z;
	int sign,ix;
	ix = (int)(x+0.5);
	sign = 1 - 2*(ix%2);
	z  = b[0] + b[1]*(NT/2-x)*(NT/2-x) + b[2]*exp(-b[3]*x);
	z +=                                 b[2]*exp(-b[3]*(NT-x));
	z += b[4]*sign*exp(-b[5]*x) + b[6]*sign*exp(-b[7]*x);
	z += b[4]*sign*exp(-b[5]*(NT-x)) + b[6]*sign*exp(-b[7]*(NT-x));
	return(z);
}

double df(double x, double *b, int i)
{
	int sign,ix;
	ix = (int)(x+0.5);
	sign = 1 - 2*(ix%2);
	if     (i==0)return( 1.0 );
	else if(i==1)return( (NT/2-x)*(NT/2-x) );
	else if(i==2)return( exp(-b[3]*x) + exp(-b[3]*(NT-x)) );
	else if(i==3)return( -b[2]*x*exp(-b[3]*x)
			     - b[2]*(NT-x)*exp(-b[3]*(NT-x)) );
	else if(i==4)return( sign*exp(-b[5]*x) + sign*exp(-b[5]*(NT-x)) );
	else if(i==5)return( -b[4]*sign*x*exp(-b[5]*x)
			     - b[4]*sign*(NT-x)*exp(-b[5]*(NT-x)) );
	else if(i==6)return( sign*exp(-b[7]*x) + sign*exp(-b[7]*(NT-x)) );
	else if(i==7)return( -b[6]*sign*x*exp(-b[7]*x)
			     - b[6]*sign*(NT-x)*exp(-b[7]*(NT-x)) );
	else return(0.0) ;
}

double ddf(double x, double *b, int i, int j)
{
	int k,sign,ix;

	if(i>j){k=i;i=j;j=k;}
	ix = (int)(x+0.5);
	sign = 1 - 2*(ix%2);

	if(i==3 && j == 3)return( b[2]*x*x*exp(-b[3]*x)
	    +b[2]*(NT-x)*(NT-x)*exp(-b[3]*(NT-x)) );
	else if(i==5 && j == 5)return( b[4]*sign*x*x*exp(-b[5]*x)
	    +b[4]*sign*(NT-x)*(NT-x)*exp(-b[5]*(NT-x)) );
	else if(i==7 && j == 7)return( b[6]*sign*x*x*exp(-b[7]*x)
	    +b[6]*sign*(NT-x)*(NT-x)*exp(-b[7]*(NT-x)) );
	else if(i==0 && j==1)return( -x*exp(-b[1]*x)
	    - (NT-x)*exp(-b[1]*(NT-x)) );
	else if(i==2 && j==3)return( -x*exp(-b[3]*x)
	    - (NT-x)*exp(-b[3]*(NT-x)) );
	else if(i==4 && j==5)return( -x*sign*exp(-b[5]*x)
	    - (NT-x)*sign*exp(-b[5]*(NT-x)) );
	else if(i==6 && j==7)return( -x*sign*exp(-b[7]*x)
	    - (NT-x)*sign*exp(-b[7]*(NT-x)) );
	else return(0.0);
}
