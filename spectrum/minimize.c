/*	Minimization using variants of Newtons		*/
/*	 and gradient descent methods			*/
/*							*/
/*	function to be minimized is 			*/
/*							*/
/*		phi( x )				*/
/*							*/
/*	whose gradient is given by			*/
/*							*/
/*		dphi( x, grad )				*/
/*							*/
/*	second partials are given by			*/
/*							*/
/*		ddphi( x , a )				*/
/*							*/
/*	minimization routine is				*/
/*							*/
/*	minimize( x, n , maxit , eps , display )	*/
/*							*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "linalg.h"
#include "minimize.h"

double minimize( double *x , int n , int maxit , double eps , int display ) {
    double **a,*b,*y;
    int	*p,*q,i,count,flag;
    double *u,*grad,*xnew,step,nstep,grad_size,tol,phi_size,new_phi_size;
    p=(int*)malloc(n*sizeof(int));
    q=(int*)malloc(n*sizeof(int));
    a=(double**)malloc(n*sizeof(double*));
    b=(double*)malloc(n*sizeof(double));
    y=(double*)malloc(n*sizeof(double));
    for(i=0;i<n;i++) a[i]=(double*)malloc(n*sizeof(double));
    grad=(double*)malloc(n*sizeof(double));
    u=(double*)malloc(n*sizeof(double));
    xnew=(double*)malloc(n*sizeof(double));
    count=0;
    nstep=1.0;
    step=1;
    tol=1.000000001;
    new_phi_size=phi_size=phi(x);
    dphi( x , grad );

    while( (grad_size=norm(grad,n)) > eps && count <= maxit  ){
	if(display){
	    printf("iteration %d\n\t|phi| =\t%e\t|grad(phi)| =\t%e\n"
		,count,phi_size,grad_size);
	    for(i=0;i<n;i++)printf("\tx[%d] =\t%e\tgrad[%d] =\t%e\n"
		,i,x[i],i,grad[i]);
	}
	ddphi( x , a );
	if(0 !=(flag=factor(a,p,q,n))) {
	    subst(a,grad,y,p,q,n);
	    if( nstep*dot( y , grad ,n) < 0.0 ) nstep = - nstep;
	    /* Newton iteration	*/
	    for(i=0;i<n;i++){
		xnew[i]= x[i] - nstep* y[i];
	    }
	    if((new_phi_size=phi(xnew))*tol <phi_size){
		/* Try to take a bigger step */
		do {
		    phi_size=new_phi_size;
		    nstep *= 1.8;
if(fabs(nstep) > 10){printf("BREAKING-BIG, phi = %e\n",phi_size); break;}
		    for(i=0;i<n;i++) xnew[i]= x[i] - nstep * y[i];
		} while( tol*phi_size > (new_phi_size=phi(xnew)) );
		new_phi_size = phi_size;
		nstep /= 1.8;
		for(i=0;i<n;i++) xnew[i]= x[i] - nstep * y[i];
	    }
	    else{
		while( !( phi_size*tol > (new_phi_size=phi(xnew)))){
		    /* if phi increased take a smaller step */
		    /* comparison is !(cond) so will take correct action
		      if new_phi_size is NaN, for which comparisons fail */
		    nstep /= 1.8;
if(fabs(nstep) < 1e-4){printf("BREAKING-SMALL, phi = %e\n",phi_size); break;}
		    for(i=0;i<n;i++) xnew[i]= x[i] - nstep * y[i];
		}
	    }

	    if( display ){
		 printf("\tNewton descent\t\tnstep =\t%e\n",nstep);
		 fflush(stdout);
	    }
	}
	else
	{	/* Gradient descent	(Last resort!!!)	*/
/**printf("min: hit a singular matrix\n");
exit(0);**/
	    for(i=0;i<n;i++){
		u[i] = -grad[i] ;
		xnew[i]= x[i] + step * u[i];
	    }
	    if( (new_phi_size=phi(xnew))*tol <phi_size){
		do {	/* try to take a bigger byte! */
printf("min: taking bigger bite\n");
		    phi_size=new_phi_size;
		    step *= 2;
		    for(i=0;i<n;i++) xnew[i]= x[i] + step * u[i];
		}
		while(  phi_size*tol > (new_phi_size=phi(xnew)) );
		new_phi_size = phi_size;
		step /=2;
		for(i=0;i<n;i++) xnew[i]= x[i] + step * u[i];
	    }
	    else {
		while(  !(tol*phi_size > (new_phi_size=phi(xnew))) && step > 0.001 ){
printf("min: taking smaller bite\n");
		    /* better take a smaller byte! */
		    step /= 2;
		    for(i=0;i<n;i++) xnew[i]= x[i] + step * u[i];
		}
	    }
	    if( display ){
		printf("\tgradient descent\t\tstep =\t%e\n",step);
		fflush(stdout);
	    }
	}
	for(i=0;i<n;i++) x[i]= xnew[i] ;
	phi_size=new_phi_size;
	dphi( x , grad );
	count++;
    }
    if(display){
	if(count >= maxit) printf("No ");
	printf("convergence after %d iterations\n|phi| = %e\t|grad(phi)| = %e\n"
	    ,count-1,phi_size,grad_size);
	printf("final parameters:\n");
	for(i=0;i<n;i++)printf("\tx[%d] =\t%e\tgrad[%d] =\t%e\n"
	    ,i,x[i],i,grad[i]);
	printf("\n");
    }

    return(new_phi_size);
}
double norm(double *x, int n)
{
    double sum;
    int i;
    sum=0;
    for(i=0;i<n;i++) sum+=x[i]*x[i];
    return( sqrt(sum) );
}
double dot(double *x, double *y, int n)
{
    double sum;
    int i;
    sum=0;
    for(i=0;i<n;i++) sum+=x[i]*y[i];
    return( sum );
}
