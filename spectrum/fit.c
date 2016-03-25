// least squares fitting for independent data points
//
// needs a file containing the function and its first and
// second derivatives
//
// least squares fitting for independent data points
//
// needs a file containing the function and its first and
// second derivatives
//
// to compile: cc -o somefunction somefunction.c fit.c linalg.c minimize.c -lm

#define NDATA	1000	/* Max number of data points */
#define NPAR	30	/* Max number of parameters  */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linalg.h"
#include "minimize.h"
double x[NPAR],xdata[NDATA],ydata[NDATA],weight[NDATA],xlow,xhigh;
int ndata,npar,fixpar,parlist[NPAR];
void graph() ;
double gammq(double a, double x) ;
void gser(double *gamser,double a,double x,double *gln) ;
void gcf(double *gammcf,double a,double x,double *gln) ;
double gammln(double xx) ;

int main(int argc,char *argv[])
{
	FILE    *file1,*file2 ;
	double **a,*b,dx[NPAR],yt[NPAR],eps,y,res,chi_square;
	int i,j,k,maxit,display,*p,*q;
	npar=NPAR;
	ndata=NDATA;
	fixpar = 0;
	b=(double*)malloc(npar*sizeof(double));
	a=(double**)malloc(npar*sizeof(double*));
	for(i=0;i<npar;i++)a[i]=(double*)malloc(npar*sizeof(double));
	p=(int*)malloc(npar*sizeof(int));
	q=(int*)malloc(npar*sizeof(int));


	if( argc!=2  &&  argc!=3) {
		printf("\tfit error\n\tcorrect format is:\n");
		printf("\tfit   datafile   [parameterfile]\n");
		exit(0);
	}
	file1 = fopen(argv[1],"r") ;
	if(file1 == NULL){
	    fprintf(stderr,"Cannot open data file %s\n",argv[1]);
	    exit(1);
	}
	if( argc==3 ) file2=fopen(argv[2],"r");
	else file2=stdin;

	/* do any initialization required by function */
	npar = -1;
	f_init(file2);

	/*read in eps */
	if(argc!=3)fprintf(stderr,"Enter epsilon, maximum iterations\n");
	fscanf(file2,"%le",&eps);

	/*read in maxit */
	fscanf(file2,"%d",&maxit);

	/*read in display flag */
	if(argc!=3)fprintf(stderr,"Enter 1 for verbose, 0 for terse\n");
	fscanf(file2,"%d",&display);

	/* read in range of x to be fitted		*/
	if(argc!=3)fprintf(stderr,"Enter range of x to be fitted\n");
	fscanf(file2,"%le%le",&xlow,&xhigh);

	/*  read in number of parameters */
	/* f_init may have already set npar */
	if(npar < 0){
	    if(argc!=3)fprintf(stderr,"Enter Number of parameters\n");
	    fscanf(file2,"%d",&npar);
	}
	else {
	    fprintf(stderr,"%d parameters\n",npar);
	}
	if(npar>NPAR){printf("To many parameters!\n");exit(0);}

	/*  read in number of parameters to be held fixed */
	if(argc!=3)fprintf(stderr,"Enter number of parameters held fixed\n");
	fscanf(file2,"%d",&fixpar);
	printf("%d of %d  parameters are held fixed\n",fixpar,npar);
	for(i=0;i<fixpar;i++){
	    if(argc!=3)fprintf(stderr,
	        "Which parameter is fixed? 0 = first parameter\n");
	    fscanf(file2,"%d",parlist+i);
	}

	/*  read in guess parameters */
	printf("\n\nInitial values:\n\n");
	for(i = 0 ; i < npar ; i++){
		if(argc!=3)fprintf(stderr,
		   "Enter initial value for parameter %d\n",i);
		if( 0 == fscanf(file2,"%le",x+i) ){
			printf("\nnot enough parameters\n");
			exit(1);
		}
		printf("x[%d]\t=\t%e\n",i,x[i]);
	}

	/*  read data	*/
	for(ndata=0;;ndata++){
		if(ndata > NDATA ){
			printf("too many data points ( >NDATA )\n") ;
			exit(1) ;
		}
		if(fscanf(file1,"%le %le %le"
		    ,xdata+ndata,ydata+ndata,weight+ndata) == EOF ) 
			goto endloop ;
		weight[ndata] = 1/(weight[ndata]*weight[ndata]) ;
		if(xlow > xdata[ndata] || xhigh < xdata[ndata])
		    ndata--;	/*reject datapoint*/
	}
endloop:
	fclose(file1) ;
	printf("\n\t read in %d data points\n",ndata);

	if(display){
		printf("data to be fitted: x,y,error\n");
		for(i=0;i<ndata;i++)printf("%e\t%e\t+/- %e\n",
		xdata[i],ydata[i],1.0/sqrt(weight[i]));
	}

	chi_square = minimize(x,npar,maxit,eps,display);

	if(display){
		/*	print differences		*/
		printf("\n\n\n");
		printf("x\t\tobserved\tpredicted\tdifference/sigma\n");
		for(i=0;i<ndata;i++){
			y=f(xdata[i],x);
			res=ydata[i]-y;
			printf("%6e\t%6e\t",xdata[i],ydata[i]);
			printf("%6e\t%f\n",y,res*sqrt(weight[i]));
		}
	}
	printf("\nchi-square with %d degrees of freedom is %e\n",
	ndata-npar+fixpar,chi_square);
	printf("confidence of fit is\t%e\n\n",
		gammq(0.5*(ndata-npar+fixpar),0.5*chi_square) );
	ddphi(x,a);
	for(i=0;i<npar;i++)dx[i]=0;
	if( 0 != factor(a,p,q,npar) ){
		for(i=0;i<ndata;i++){
			for(j=0;j<npar;j++){
				b[j]=0;
				for(k=0;k<fixpar;k++)if(j==parlist[k])goto bottom;
				b[j] = 2*weight[i] * df(xdata[i],x,j); 
bottom:	
				;
			}
			subst(a,b,yt,p,q,npar);
			for(j=0;j<npar;j++) dx[j] += yt[j]*(yt[j]/weight[i]);
		}
		for(i=0;i<npar;i++)dx[i]=sqrt(dx[i]);
	}
	else printf("Hessian is singular\tno error estimates\n");
	printf("\n\nFINAL PARAMETERS:\n\n");
	for(i = 0 ; i < npar ; i++){
	    for(j=0;j<fixpar;j++)if(parlist[j]==i)break;
	    if(j==fixpar)printf("x[%d]\t=\t%e\t+-\t%e\n",i,x[i],dx[i]);
	    else       printf("x[%d]\t=\t%e\t+-\t(fixed)\n",i,x[i]);
	}
	graph();
	return 0 ;
}

void graph(){
/* Routine to draw experimental data and compare to fit	*/
int i;
double tx,ty,z,xmin=1e30,xmax= -1e30;
FILE *fp ;
	fp=fopen("fit_graph.ax","w");
	for(i=0;i<ndata;i++){
		z = xdata[i];
		if( z < xmin )xmin = z;
		if( z > xmax )xmax = z;
	}
	z = (xmax-xmin)/(2*ndata);
	xmin -= z; xmax += z;
	tx = (xmax-xmin) / 100. ;
	/* data and error bars */
	for(i=0;i<ndata;i++){
		ty= sqrt( 1/weight[i]);
		fprintf(fp,"%e %e\n%e %e\n#\n",
			xdata[i]-tx,ydata[i]+ty,xdata[i]+tx,ydata[i]+ty);
		fprintf(fp,"%e %e\n%e %e\n#\n",
			xdata[i]-tx,ydata[i]-ty,xdata[i]+tx,ydata[i]-ty);
		fprintf(fp,"%e %e\n%e %e\n#\n",
			xdata[i],ydata[i]-ty,xdata[i],ydata[i]+ty);
	}
	for(i=0;i<101;i++){
		z=xmin + i*tx;
		fprintf(fp,"%e %e\n",z,f(z,x));
	}
	fclose(fp);
}

double phi(double *b)
{
	double sum,diff;
	int i;
	sum=0;
	for(i=0,sum=0 ;i<ndata;i++){
		diff= f( xdata[i] , b ) - ydata[i];
		sum += weight[i] * diff * diff ;
	}
	return(sum);
}

void dphi(double *b, double *grad)
{
	double diff;
	int i,j,k;
	for(j=0;j<npar;j++) {
		grad[j]=0;
		for(k=0;k<fixpar;k++)if(j==parlist[k])goto bottom;
		for(i=0;i<ndata;i++){
			diff= f( xdata[i] , b ) - ydata[i];
			grad[j] += 2*weight[i] * diff * df(xdata[i],b,j); 
		}
bottom:	
		;
	}
}

void ddphi(double *b, double **a)
{
	int i,j,k,l;
	for(j=0;j<npar;j++)for(k=0;k<=j;k++){
		a[j][k]=0;
		for(l=0;l<fixpar;l++)if(j==parlist[l] || k ==parlist[l] ){
			if(j==k)a[j][k]=1;
			goto bottom;
		}
		for(i=0;i<ndata;i++) a[j][k] += 2*weight[i]*
		    ( df(xdata[i],b,j) * df(xdata[i],b,k) + 
		    (f(xdata[i],b)-ydata[i]) * ddf(xdata[i],b,j,k) );
bottom:	
		;
		a[k][j] = a[j][k] ;
	}
}

/* To get confidence: gammq( 0.5*dof , 0.5*chisq ) */
/* From Numerical Recipes:  */
double gammq(double a, double x)
{
	double gamser,gammcf,gln;
	void gcf(),gser();

	if (x < 0.0 || a <= 0.0)
	    fprintf(stderr,"Invalid arguments in routine GAMMQ");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}

#define ITMAX 200
#define EPS 1.0e-10

void gser(double *gamser,double a,double x,double *gln)
{
	int n;
	double sum,del,ap;
	double gammln();

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) fprintf(stderr,"x less than 0 in routine GSER");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			ap += 1.0;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		fprintf(stderr,"a too large, ITMAX too small in routine GSER");
		return;
	}
}

void gcf(double *gammcf,double a,double x,double *gln)
{
	int n;
	double gold=0.0,g,fac=1.0,b1=1.0;
	double b0=0.0,anf,ana,an,a1,a0=1.0;
	double gammln();

	*gln=gammln(a);
	a1=x;
	for (n=1;n<=ITMAX;n++) {
		an=(double) n;
		ana=an-a;
		a0=(a1+a0*ana)*fac;
		b0=(b1+b0*ana)*fac;
		anf=an*fac;
		a1=x*a0+anf*a1;
		b1=x*b0+anf*b1;
		if (a1) {
			fac=1.0/a1;
			g=b1*fac;
			if (fabs((g-gold)/g) < EPS) {
				*gammcf=exp(-x+a*log(x)-(*gln))*g;
				return;
			}
			gold=g;
		}
	}
	fprintf(stderr,"a too large, ITMAX too small in routine GCF");
}

#undef ITMAX
#undef EPS

double gammln(double xx)
{
	double x,tmp,ser;
	static double cof[6]={76.18009173,-86.50532033,24.01409822,
		-1.231739516,0.120858003e-2,-0.536382e-5};
	int j;

	x=xx-1.0;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.0;
	for (j=0;j<=5;j++) {
		x += 1.0;
		ser += cof[j]/x;
	}
	return -tmp+log(2.50662827465*ser);
}
