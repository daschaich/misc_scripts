/* correlated least squares fit (designed for QCD hadron propagators).	*
*									*
* Fitting function and derivatives must be supplied as f(),df(),ddf().	*
* Read in average function and covariance matrix			*
* Then minimize quadratic form:						*
*  (y[i]-f(x[i])) * (inverse of covariance matrix[ij]) * (y[j]-f(x[j]))	*
*									*
* Form of input:							*
* "PROPAGATOR"	(or some other word)
* list of function points:
*	x y
* "COVARIANCE"  blocks   "blocks" (or some other words), where blocks
*     is the number of blocks used in computing the covariance matrix.
* list of covariance matrix elements:
*	x1 x2 covar12
* Missing covariance matrix elements are taken to be zero,
*  except diagonal ones must be specified.
* Extra covariance matrix elements are ignored.
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "linalg.h"
#include "minimize.h"
double Fdist_2( double chisq, int N, int DOF);
double conf_int( double chisq, int N, int DOF );

void dumpmat(char *label,char *format,double *matrix) ;
void graph() ;
double gammq(double a, double x) ;
void gser(double *gamser,double a,double x,double *gln) ;
void gcf(double *gammcf,double a, double x,double *gln) ;
double gammln(double xx) ;

double *par;		/* parameter array */
double *delpar,*wparmat1,*wparmat2;	/* errors on parameters, workspace */
double *xdata,*ydata;	/* x and y values, x are distances */
double *covar,*covarinv;	/* covariance matrix and inverse */
double xlow,xhigh;/* range of distance to include */
double *wdata1,*wdata2,*wdata3,*wdata4,*wpar1,*wpar2;
    /* working arrays, ndata & npar points*/
int nt,npar,fixpar,*parlist,ndata;
int nblocks;	/* number of blocks used in computing covariance matrix */
int ndof;	/* number of degrees of freedom in the fit */
double gammq();

int main(int argc,char *argv[]) {
        FILE    *file1,*file2;
	double **a,*b,*dpar,*yt,eps,y,res,chi_square,chi_square_adjusted;
	double xx,yy,zz;
	int i,j,k,l,maxit,display,*p,*q;
	fixpar = 0;

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

	/*  read in number of parameters, allocate storage */
	/* f_init may have already set npar */
	if(npar < 0){
	    if(argc!=3)fprintf(stderr,"Enter Number of parameters\n");
	    fscanf(file2,"%d",&npar);
	}
	else {
	    fprintf(stderr,"%d parameters\n",npar);
	}
	p=(int*)malloc(npar*sizeof(int));
	q=(int*)malloc(npar*sizeof(int));
	par=(double*)malloc(npar*sizeof(double));
	wpar1=(double*)malloc(npar*sizeof(double));
	wpar2=(double*)malloc(npar*sizeof(double));
	delpar=(double*)malloc(npar*npar*sizeof(double));
	wparmat1=(double*)malloc(npar*npar*sizeof(double));
	wparmat2=(double*)malloc(npar*npar*sizeof(double));
	parlist=(int*)malloc(npar*sizeof(int));
	dpar=(double*)malloc(npar*sizeof(double));
	yt=(double*)malloc(npar*sizeof(double));
	b=(double*)malloc(npar*sizeof(double));
	a=(double**)malloc(npar*sizeof(double*));
	for(i=0;i<npar;i++)a[i]=(double*)malloc(npar*sizeof(double));

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
		if( 0 == fscanf(file2,"%le",par+i) ){
			printf("\nnot enough parameters\n");
			exit(0);
		}
		printf("par[%d]\t=\t%e\n",i,par[i]);
	}

	/*  read data	*/
	fscanf(file1,"%*s");	/* flush header word */
	xdata = (double *)malloc(sizeof(double));
	ydata = (double *)malloc(sizeof(double));
	ndata=0;
	while(fscanf(file1,"%le %le",&xx,&yy) == 2 ){
	    if(xx < xlow || xx > xhigh ) continue;	/*reject datapoint*/
	    xdata = (double *)realloc(xdata,(ndata+1)*sizeof(double));
	    ydata = (double *)realloc(ydata,(ndata+1)*sizeof(double));
	    xdata[ndata] = xx;   ydata[ndata] = yy;
	    ndata++;
	}
	wdata1 = (double*)malloc(ndata*sizeof(double));
	wdata2 = (double*)malloc(ndata*sizeof(double));
	wdata3 = (double*)malloc(ndata*sizeof(double));
	wdata4 = (double*)malloc(ndata*sizeof(double));
	ndof = ndata-npar+fixpar;

	/* read covariances */
	fscanf(file1,"%*s%d%*s",&nblocks);	/* flush header words */
	printf("\n\t read in %d data points\n",ndata);
	if( nblocks <= ndata ){
	    fprintf(stderr,"Not enough blocks in covariance matrix %d  %d\n",ndata,nblocks);
	    exit(0);
	}
	/* initialize matrix */
	covar = (double*)malloc(ndata*ndata*sizeof(double));
	covarinv = (double*)malloc(ndata*ndata*sizeof(double));
	for(i=0;i<ndata;i++){
	    for(j=0;j<ndata;j++)covar[i*ndata+j] = 0.0;
	    covar[i*ndata+i] = -1.0;
	    /*diagonal can not be < 0, use to make sure it is specified */
	}
	while(fscanf(file1,"%le%le%le",&xx,&yy,&zz) == 3){
	    /* find indices of this element */
	    for(i=0;i<ndata;i++)if(xx == xdata[i])break;
	    for(j=0;j<ndata;j++)if(yy == xdata[j])break;
	    if(i < ndata && j < ndata)covar[i*ndata+j] = zz;
/****if(i != j )covar[i*ndata+j] = 0.0;	 temporary test */
	}
	fclose(file1);
	/* check to make sure diagonal elements were specified */
	for(i=0,j=0;i<ndata;i++){
	   if(covar[i*ndata+i] < 0.0){
	       fprintf(stderr,
		  "Diagonal covariance missing or negative: x = %e\n",xdata[i]);
	       j=1;
	    }
	}
	if(j)exit(1);

	/*invert covariance matrix */
	matinv(covar,covarinv,ndata);

	if(display){
		printf("data to be fitted: x,y,error\n");
		for(i=0;i<ndata;i++)printf("%e\t%e\t+/- %e\n",
		xdata[i],ydata[i],sqrt(covar[i*ndata+i]));
	}

	chi_square = minimize(par,npar,maxit,eps,display);

	if(display){
		/*	print differences		*/
		printf("\n\n\n");
		printf("x\t\tobserved\tpredicted\tdifference/naive_sigma\n");
		for(i=0;i<ndata;i++){
			y=f(xdata[i],par);	
			res=ydata[i]-y;
			printf("%6e\t%6e\t",xdata[i],ydata[i]);
			printf("%6e\t%f\n",y,res/sqrt(covar[i*ndata+i]));
		}
	}

	/* make error estimates on parameters */

	// 3/08 update the correction factors for finite sample size
	// see klingon:/home/statistics/largeNnotes.tex for derivation
	// Also, the mean and variance of chi^2 are changed
	// compute an adjusted chi^2 with the mean and variance expected
	// in uncorrelated fitting, and use that to compute an adjusted
	// confidence level

	/* first look at derivatives of chi-squared */
	/* essentially ask how they vary as chi-squared goes to one more
	   than its minimum value */
	printf("\nERROR ANALYSIS FROM CHI-SQUARED CONTOURS\n");
	ddphi(par,a);
	/* second derivatives of chi-square. Unit matrix for fixed params */
	for(i=0;i<npar;i++)for(j=0;j<npar;j++)wparmat1[i*npar+j]=a[i][j];
	matinv(wparmat1,delpar,npar);
	dumpmat("Hessian matrix","%e",wparmat1);
	dumpmat("Inverse Hessian matrix","%e",delpar);
	for(i=0;i<npar;i++){
	    for(j=0;j<npar;j++)wparmat2[i*npar+j] = 
		delpar[i*npar+j]/sqrt(delpar[i*npar+i]*delpar[j*npar+j]);
	}
	dumpmat("Scaled inverse Hessian matrix","%.3f",wparmat2);

	/* now look at dependence of fit parameters on data points */
	printf("\nERROR ANALYSIS FROM FULL EXPRESSION\n");
	/* wparmat1 = df/dpar * Cinv * df/dpar */
	for(i=0;i<npar;i++)for(j=0;j<npar;j++){
	    for(l=0;l<fixpar;l++)if(i==parlist[l] || j ==parlist[l] ){
		if(i==j)wparmat1[i*npar+j] = 1.0; else wparmat1[i*npar+j]=0.0;
		goto bottom2;
	    }
	    wparmat1[i*npar+j] = 0.0;
	    for(k=0;k<ndata;k++){
		wdata1[k] = df(xdata[k],par,i);
		wdata2[k] = df(xdata[k],par,j);
	    }
	    for(k=0;k<ndata;k++)for(l=0;l<ndata;l++)
		wparmat1[i*npar+j] += wdata1[k]*covarinv[k*ndata+l]*wdata2[l];
bottom2:    ;
	}
	/* wparmat2 = wparmat1 * delpar (delpar = inverse Hessian) */
	for(i=0;i<npar;i++)for(j=0;j<npar;j++)
	    for(wparmat2[i*npar+j]=0.0,k=0;k<npar;k++)
	    wparmat2[i*npar+j] += wparmat1[i*npar+k]*delpar[k*npar+j];
	/* wparmat1 = delpar * wparmat1 (delpar = inverse Hessian) */
	for(i=0;i<npar;i++)for(j=0;j<npar;j++)
	    for(wparmat1[i*npar+j]=0.0,k=0;k<npar;k++)
	    wparmat1[i*npar+j] += 4.0*delpar[i*npar+k]*wparmat2[k*npar+j];
	dumpmat("Parameter variance matrix","%e",wparmat1);
	for(i=0;i<npar;i++){
	    for(j=0;j<npar;j++)wparmat2[i*npar+j] = 
		wparmat1[i*npar+j]/sqrt(wparmat1[i*npar+i]*wparmat1[j*npar+j]);
	}
	dumpmat("Scaled parameter variance matrix","%.3f",wparmat2);

	printf("\nFINAL PARAMETERS, CHI_SQUARE ERRORS, FULL ERRORS\n");
	printf("chi-square with %d degrees of freedom is %e\n", ndof,chi_square);
	printf("confidence of fit is\t%e\n\n", gammq(0.5*(ndof),0.5*chi_square) );
	chi_square_adjusted = (ndof) + 
		(chi_square - (ndof)*(1+((double)ndof+2)/(double)nblocks) )
		* sqrt( ((double)nblocks)/((double)nblocks+3*ndof+6) );
	printf("adjusted-chi-square with %d degrees of freedom is %e\n",
	        ndof,chi_square_adjusted);
	printf("adjusted-confidence of fit is\t%e\n\n",
		gammq(0.5*(ndof),0.5*chi_square_adjusted) );
        printf("real-confidence of fit is \t%e\n\n",
                conf_int( chi_square, nblocks, ndof ) );

	for(i=0;i<npar;i++){
	    for(j=0;j<fixpar;j++)if(parlist[j]==i)break;
	    if(j!=fixpar) printf("par[%d] = %e  +-  (fixed)\n",i,par[i]);
	    else printf("par[%d] = %e  +-  %e \t%e\n",
		i,par[i],
		sqrt( (nblocks+ndof)/(double)(nblocks-(1.0+ndof)) * 2.0*delpar[i*npar+i] ),
		sqrt( (nblocks+ndof)/(double)(nblocks-(1.0+ndof)) * wparmat1[i*npar+i])        );
		//sqrt( (double)nblocks/(nblocks-ndata) * 2.0*delpar[i*npar+i] ),
		//sqrt( (double)nblocks/(nblocks-ndata) * wparmat1[i*npar+i])        );
		// see intro comments for source of this adjustment
	}

	graph();
}

/* dump a matrix in parameter space, skip fixed parameters */
void dumpmat(char *label,char *format,double *matrix)
{
int i,j,k;
    printf("\n%s\n",label);
    for(i=0;i<npar;i++){
	for(k=0;k<fixpar;k++)if(parlist[k]==i)break;if(k!=fixpar)continue;
	for(j=0;j<npar;j++){
	    for(k=0;k<fixpar;k++)if(parlist[k]==j)break;if(k!=fixpar)continue;
	    printf(format,matrix[i*npar+j]);printf(" \t");
	}
	printf("\n");
    }
}

void graph()
{
/* Routine to draw experimental data and compare to fit	*/
int i;
double tx,z,xmin=1e30,xmax= -1e30;
FILE *fp;
	fp=fopen("cofitpic","w");
	for(i=0;i<ndata;i++){
		z = xdata[i];
		if( z < xmin )xmin = z;
		if( z > xmax )xmax = z;
	}
	z = (xmax-xmin)/(2*ndata);
	xmin -= z; xmax += z;
	tx = (xmax-xmin) / 100. ;
	/* data and error bars */
	fprintf(fp,"# c \"\\oc\"\n");
	for(i=0;i<ndata;i++){
		fprintf(fp,"%e %e\n",xdata[i],ydata[i]);
	}
	fprintf(fp,"# c0\n");
	/* fitting function */
	for(i=0;i<101;i++){
		z=xmin + i*tx;
		fprintf(fp,"%e %e\n",z,f(z,par));
	}
	fclose(fp);
}

double phi(double *b)
{
    double sum;
    int i,j;
    sum=0.0;
    for(i=0,sum=0.0 ;i<ndata;i++)for(j=0;j<ndata;j++){
	sum += (f( xdata[i] , b ) - ydata[i]) * covarinv[i*ndata+j] 
	    * (f( xdata[j] , b ) - ydata[j]) ;
    }
    return(sum);
}

void dphi(double *b, double *grad)
{
    double deriv;
    int i,j,k;
    for(i=0;i<ndata;i++) wdata1[i] = f(xdata[i],b) - ydata[i];
	/* compute difference outside double loop */
    for(j=0;j<npar;j++) {
	grad[j]=0.0;
	for(k=0;k<fixpar;k++)if(j==parlist[k])goto bottom;
	for(i=0;i<ndata;i++){
	    deriv = df(xdata[i],b,j);
	    for(k=0;k<ndata;k++){
	        grad[j] += 2.0 * deriv * covarinv[i*ndata+k] * wdata1[k];
	    }
	}
bottom:	
	;
    }
}

void ddphi(double *b, double **a) 
{
    int i,j,k,l;
    for(i=0;i<ndata;i++)wdata4[i] = f(xdata[i],b) - ydata[i];
    for(j=0;j<npar;j++){
	for(i=0;i<ndata;i++)wdata1[i] = df(xdata[i],b,j);
	for(k=0;k<=j;k++){
	    a[j][k]=0.0;
	    for(l=0;l<fixpar;l++)if(j==parlist[l] || k ==parlist[l] ){
	        if(j==k)a[j][k]=1.0;
	        goto bottom;
	    }
	    for(i=0;i<ndata;i++){
		wdata2[i] = df(xdata[i],b,k);
		wdata3[i] = ddf(xdata[i],b,j,k);
	    }
	    for(i=0;i<ndata;i++){
	        for(l=0;l<ndata;l++){
	            a[j][k] += 2.0*(wdata1[i]*covarinv[i*ndata+l]*wdata2[l] + 
	            wdata4[i] * covarinv[i*ndata+l] * wdata3[l] );
	            /**a[j][k] += 2.0 *
	            ( df(xdata[i],b,j)*covarinv[i*ndata+l] * df(xdata[l],b,k) + 
	            (f(xdata[i],b)-ydata[i]) * covarinv[i*ndata+l] *
		    ddf(xdata[l],b,j,k) );**/
	        }
	    }
bottom:     ;
	    a[k][j] = a[j][k] ;
	}
    }
}

/* From Numerical Recipes:  */
/* To get confidence: gammq( 0.5*ndof , 0.5*chisq ) */
double gammq(double a, double x)
{
	double gamser,gammcf,gln;
	void gcf(),gser();

	if( a==0.0 )return(0.0);
	if (x < 0.0 || a <= 0.0)
	    fprintf(stderr,"Invalid arguments in routine GAMMQ %e %e\n",
		a,x);
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}

// confidence = 1.0 - integral from 0 to chisq of Fdist_2
double conf_int( double chisq, int N, int DOF ){
	double x,y,z, sum, eps;
	double chisq_t;
	int nsteps;

	    if(chisq <= 0.0) return(1.0); // chisq < 0 maybe should be an error
	    if(N<=DOF || DOF <= 0 )return(0.0);  // again, probably really an error

	    // if confidence is not obviously small, try integrating from zero
	    if( chisq < 3.0*DOF+2.0 ){
	        nsteps = (int)( 100.0*chisq/(2.0*sqrt(2.0*DOF)) ); // about 100 steps in peak
//fprintf(stderr,"CONF_INT_1: chisq = %e  N = %d   DOF = %d   nsteps = %d\n",chisq,N,DOF,nsteps);
	        if(nsteps<20)nsteps=20;
	        eps = chisq/(double)nsteps;
	        for( sum=0.0, chisq_t=0.5*eps; chisq_t<chisq; chisq_t+= eps ){
	            sum += Fdist_2( chisq_t, N, DOF );
	        }
	        if( sum*eps < 0.9 )return( 1.0 - sum*eps );
	    }

	    // if confidence is small, integrate upper part of distribution
	    nsteps=0;
	    eps = (2.0*sqrt(2.0*DOF))/100.0;
//fprintf(stderr,"CONF_INT_2: chisq = %e  N = %d   DOF = %d   eps = %e\n",chisq,N,DOF,eps);
	    y = Fdist_2( chisq, N, DOF );
	    for( nsteps=0,sum=0.0, chisq_t=chisq+0.5*eps;  ; chisq_t+= eps,nsteps++ ){
	        z = Fdist_2( chisq_t, N, DOF );
		if( nsteps>1000 || z < 0.0001*y )break;
	        sum += z;
	    }
	    return( sum*eps );
}

// probability density - F distribution   F_{m,n}
// m=DOF     n=(N-DOF)
// using the more natural (for us) variables N and DOF
double Fdist_2( double chisq, int N, int DOF){
    double Nd,Dd;
    Nd=(double)N; Dd=(double)DOF;
    return(  exp( gammln((Nd)/2.0) - gammln(Dd/2.0) - gammln((Nd-Dd)/2.0)
        + (-Dd/2.0)*log(Nd)
        + ((Dd-2.0)/2.0)*log(chisq)
        + (-Nd/2.0)*log(1.0+chisq/Nd )
        ) );
}
#define ITMAX 200
#define EPS 1.0e-10

void gser(double *gamser,double a,double x,double *gln)
{
	int n;
	double sum,del,ap;

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

void gcf(double *gammcf,double a, double x,double *gln)
{
	int n;
	double gold=0.0,g,fac=1.0,b1=1.0;
	double b0=0.0,anf,ana,an,a1,a0=1.0;

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
