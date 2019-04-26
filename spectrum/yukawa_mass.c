/* for use with "fit.c".
* fit to form b[0]*exp(-b[1]*x)/x**b[3]
*             + "image"
*
* compile with cofit_n.c min.c linalg.c
*/
#include <stdio.h>
#include <math.h>
#include "minimize.h"

extern int npar;
int NT;

void f_init(FILE *fp)
{
  fprintf(stderr,"Kogut-Susskind meson mass fitting\n");
  fprintf(stderr,"Enter time size of lattice:\n");
  fscanf(fp,"%d",&NT);
  npar = 8;
}

double f(double x, double *b)
/* b[0] = coefficient 1
*  b[1] = mass 1 (natural parity)
*  b[2] = Yukawa exponent
*/
{
	double z;
	int sign,ix;
	ix = (int)(x+0.5);
	sign = 1 - 2*(ix%2);
	z  = b[0]*exp(-b[1]*x)*exp(-b[2]*log(x));
	z += b[0]*exp(-b[1]*(NT-x))*exp(-b[2]*log(NT-x));
	return(z);
}

double df(double x, double *b, int i)
{
	int sign,ix;
        double xx,xxt;
	ix = (int)(x+0.5);
	sign = 1 - 2*(ix%2);
        xx = exp(-b[1]*x)*exp(-b[2]*log(x));
        xxt = exp(-b[1]*(NT-x))*exp(-b[2]*log(NT-x));
	if     (i==0)return( xx + xxt ) ;
	else if(i==1)return( -b[0]*x*xx
			     - b[0]*(NT-x)*xxt );
	else if(i==2)return(  -b[0]*log(x)*xx
                             - b[0]*log(NT-x)*xxt );
	else if(i==3)return(0.0);
	else if(i==4)return( 0.0);
	else if(i==5)return( 0.0);
	else if(i==6)return( 0.0);
	else if(i==7)return( 0.0);
	else return(0.0) ;
}

double ddf(double x, double *b, int i, int j)
{
	int k,sign,ix;
        double xx,xxt;

	if(i>j){k=i;i=j;j=k;}
	ix = (int)(x+0.5);
	sign = 1 - 2*(ix%2);
        xx = exp(-b[1]*x)*exp(-b[2]*log(x));
        xxt = exp(-b[1]*(NT-x))*exp(-b[2]*log(NT-x));

	if     (i==1 && j == 1)return( b[0]*x*x*xx
	    +b[0]*(NT-x)*(NT-x)*xxt );
	else if(i==2 && j==2)return( b[0]*log(x)*log(x)*xx
                               +b[0]*log(NT-x)*log(NT-x)*xxt );
	else if(i==0 && j==1)return( - x*xx
                                - (NT-x)*xxt );
	else if(i==0 && j==2)return( - log(x)*xx
                               - log(NT-x)*xxt );
	else if(i==1 && j==2)return( b[0]*x*log(x)*xx
                               +b[0]*(NT-x)*log(NT-x)*xxt );
	else return(0.0);
}
