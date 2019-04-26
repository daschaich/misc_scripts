
void lineq( double *mat, double *vec, double *ans, int dim );
void matinv( double *x, double *y, int dim );
int factor( double **a, int *p, int *q, int n );
void subst( double **a, double *b, double *x, int *p, int *q, int n );
void trisolve( double *lower, double *diag, double *upper, double *vec,
	       double *ans, int dim );
