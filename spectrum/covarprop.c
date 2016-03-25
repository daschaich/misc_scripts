/* read a list of propagators, compute average propagator and
   covariance matrix.
   Propagator is list of numbers
      0   value
      1   value
      2   value
      ...
      max value
  Output is
    PROPAGATOR
    i  average[i]
    COVARIANCES
    i j covar[i,j]
   Arguments are maximum distance of propagator and blocking factor
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void main(int argc,char *argv[]) {
int blocksize,nblocks,npoints,nprops,dist,i,j;
double *covar,*average,*current,*block,x;

    /* find size of propagator */
    if(argc < 3){fprintf(stderr,"Not enough arguments\n"); exit(1);}
    if(sscanf(argv[1],"%d",&npoints) != 1){
	fprintf(stderr,"Bad argument (maximum distance)\n"); exit(1);
    }
    if(sscanf(argv[2],"%d",&blocksize) != 1){
	fprintf(stderr,"Bad argument (blocksize)\n"); exit(1);
    }
    npoints++;	/* zero distance counts! */
    /* allocate storage */
    average = (double*)malloc(npoints*sizeof(double));
    current = (double*)malloc(npoints*sizeof(double));
    block = (double*)malloc(npoints*sizeof(double));
    covar = (double*)malloc(npoints*npoints*sizeof(double));
    /* initialize block, average and covariance */
    for(i=0;i<npoints;i++){
	block[i] = average[i]=0.0;
	for(j=0;j<npoints;j++)covar[i*npoints+j]=0.0;
    }

    /* read propagators, accumulate sums, check for proper format and order */
    for(nblocks=0,nprops=0;;){
        for(dist=0;dist<npoints;dist++){
	    if((j=scanf("%d%le",&i,current+dist)) != 2){
		if(dist == 0){	/* end of file, probably */
		    goto done;	/* break two loops */
		}
		else{
		    fprintf(stderr,"Bad data in propagator %d\n",nprops);
		    exit(1);
		}
	    }
	    if(i != dist){	/* unexpected entry - give no mercy */
		fprintf(stderr,"Bad data in propagator %d\n",nprops);
		exit(1);
	    }
	}
	for(i=0;i<npoints;i++){
	    block[i] += current[i];
	}
	nprops++;
	if(nprops % blocksize == 0){	/* we have finished a block */
	    for(i=0;i<npoints;i++)block[i] /= (double)blocksize;
	    for(i=0;i<npoints;i++){
	        average[i] += block[i];
	        for(j=0;j<npoints;j++)covar[i*npoints+j]
		    += block[i]*block[j];
	    }
	    for(i=0;i<npoints;i++)block[i] = 0.0;
	    nblocks++;
	}
    }
done:
    fprintf(stderr,"%d propagators read in, %d blocks\n",nprops,nblocks);

    /* normalize and subtract disconnected parts */
    /* covariance matrix is normalized for "error of mean" */
    /* no attempt is made to correct for difference between the
	sample average and the real average - subsequent fitting
	programs must take care of this. */
    x = 1.0/(double)nblocks;
    for(i=0;i<npoints;i++) average[i] *= x;
    for(i=0;i<npoints;i++)for(j=0;j<npoints;j++) covar[i*npoints+j] =
	(x*covar[i*npoints+j] - average[i]*average[j])/((double)(nblocks));

    /* print results */
    printf("PROPAGATOR\n");
    for(i=0;i<npoints;i++)printf("%d\t%.15e\n",i,average[i]);
    printf("COVARIANCE:  %d blocks\n",nblocks);
    for(i=0;i<npoints;i++)for(j=0;j<npoints;j++)
	printf("%d\t%d\t%.15e\n",i,j,covar[i*npoints+j]);
}

