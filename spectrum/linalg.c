
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "linalg.h"

/* compute inverse of dim by dim matrix at x, put result at y */
void matinv( double *x, double *y, int dim ){
  double **a,*b,*c;
  int *perm,*perm2,i,j;

  perm = (int *)malloc(dim*sizeof(int));
  perm2 = (int *)malloc(dim*sizeof(int));
  a = (double **)malloc(dim*sizeof(double *));
  b = (double *)malloc(dim*sizeof(double));
  c = (double *)malloc(dim*sizeof(double));
  a[0] = (double *)malloc(dim*dim*sizeof(double));
  for( i=0; i<dim; i++ ) a[i] = a[0] + (i*dim);

  for( i=0; i<dim*dim; i++ ) a[0][i] = x[i];
  if( factor(a,perm,perm2,dim) == 0 ){
    printf("singular matrix\n");
    exit(0);
  }
  for( i=0; i<dim; i++ ){
    for( j=0; j<dim; j++ ) b[j] = 0.0;
    b[i] = 1.0;
    subst(a,b,c,perm,perm2,dim);
    for( j=0; j<dim; j++ ) y[j*dim+i] = c[j];
  }

  free(a[0]); free(a); free(b); free(c); free(perm); free(perm2);
}


/* solve the linear equation mat*ans = vec */
void lineq( double *mat, double *vec, double *ans, int dim ){
  double **a;
  int *perm,*perm2,i;

  perm = (int *)malloc(dim*sizeof(int));
  perm2 = (int *)malloc(dim*sizeof(int));
  a = (double **)malloc(dim*sizeof(double *));
  a[0] = (double *)malloc(dim*dim*sizeof(double));
  for( i=0; i<dim; i++ ) a[i] = a[0] + (i*dim);

  for( i=0; i<dim*dim; i++ ) a[0][i] = mat[i];
  if( factor(a,perm,perm2,dim) == 0 ){
    printf("singular matrix\n");
    exit(0);
  }
  subst(a,vec,ans,perm,perm2,dim);

  free(a[0]); free(a); free(perm); free(perm2);
}


/*      Gaussian elimination with scaled TOTAL pivoting */
#include        <stdio.h>
#define         max(a,b)        ( ((a)>(b))?(a) :(b) )
#define         abs(x)          ( ((x)>0)? (x) : -(x) )

int factor( double **a, int *p, int *q, int n ){
  int i,j,k,l,m,temp;
  double  big,x,*d,eps;

  eps = n*1e-16;

  /* get space for row norms */
  if( NULL == ( d=(double*)malloc(n*sizeof(double)) ) ){
    printf("no space!!\n");
    return(0);
  }
  for( i=0; i<n; i++ ){
    p[i] = q[i] = i;    /*      initialize identity permutation */
    big = 0.0;
    for( j=0; j<n; j++ ) big = max( big, abs(a[i][j]) );
    if( big == 0.0 ) return(0);
    d[i] = big;       /*      store row norms */
  }
  for( k=0; k<n-1; k++ ){
    j = k; m = k;
    big = 0.0;
    for( i=k; i<n; i++ ) for( l=k; l<n; l++ ){
      x = abs(a[p[i]][q[l]]/d[p[i]]);
      if( x > big ){
	j = i;            /* find pivoting element*/
	m = l;
	big = x;
      }
    }
    if( big == 0.0 )return(0);
    temp = p[k];
    p[k] = p[j];      /* keep track of permutation    */
    p[j] = temp;
    temp = q[k];
    q[k] = q[m];      /* keep track of permutation    */
    q[m] = temp;
    for( i=k+1; i<n; i++){
      /* store row multpliers in would-be zeroed location */
      x = (a[p[i]][q[k]] /= a[p[k]][q[k]]);
      for( j=k+1; j<n; j++) a[p[i]][q[j]] -= x*a[p[k]][q[j]];
    }
  }
  free(d);
  if( abs(a[p[n-1]][q[n-1]]/d[p[n-1]]) < eps ) return(0);
  else return(1);
}


void subst( double **a, double *b, double *x, int *p, int *q, int n ){
  int k,j;
  double sum;

  for( k=0; k<n; k++ ){       /* foward elimination */
    sum = 0.0;
    for( j=0; j<k; j++ ) sum += a[p[k]][q[j]] * x[q[j]];
    x[q[k]] = b[p[k]] - sum;
  }
  for( k=n-1; k>=0; k-- ){    /* back substitution */
    sum = 0.0;
    for( j=k+1; j<n; j++ ) sum += a[p[k]][q[j]] * x[q[j]];
    x[q[k]] = ( x[q[k]] - sum )/a[p[k]][q[k]];
  }
}


void trisolve( double *lower, double *diag, double *upper, double *vec,
	       double *ans, int dim ){
  /* solves a tridiagonal system of linear eqs.:

     / diag[0]   upper[0]                                \ / ans[0]     \
     | lower[1]  diag[1]   upper[1]                      | | ans[1]     |
     |                     ...                           | |  ...       |  =
     |           lower[dim-2]  diag[dim-2]   upper[dim-2]| | ans[dim-2] |
     \                         lower[dim-1]  diag[dim-1] / \ ans[dim-1] /

     / vec[0]     \
     | vec[1]     |
     |  ...       |
     | vec[dim-2] |
     \ vec[dim-1] /

  */

  int i;
  double *beta,*rho;

  beta = (double *)malloc(dim*sizeof(double));
  rho = (double *)malloc(dim*sizeof(double));

  if( diag[0] == 0.0 ){
    fprintf(stderr,"zero diag element in trisolve\n");
    exit(1);
  }

  /* forward elimination */
  beta[0] = diag[0];
  rho[0] = vec[0];
  for( i=1; i<dim; i++ ){
    beta[i] = diag[i] - lower[i]*upper[i-1]/beta[i-1];
    rho[i] = vec[i] - lower[i]*rho[i-1]/beta[i-1];
    if( beta[i] == 0.0 ){
      fprintf(stderr,"zero diag element in trisolve\n");
      exit(1);
    }
  }

  /* back substitution */
  ans[dim-1] = rho[dim-1]/beta[dim-1];
  for( i=2; i<=dim; i++ ){
    ans[dim-i] = (rho[dim-i] - upper[dim-i]*ans[dim-i+1])/beta[dim-i];
  }

  free(beta); free(rho);
}
