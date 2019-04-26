void f_init(FILE *fp) ;
double f(double x, double *b) ;
double df(double x, double *b, int i);
double ddf(double x, double *b, int i, int j) ;

double phi(double *b) ;
void dphi(double *b, double *grad) ;
void ddphi(double *b, double **a)  ;


double minimize( double *x , int n , int maxit , double eps , int display ) ;
double dot(double *x, double *y, int n) ;
double norm(double *x, int n) ;
