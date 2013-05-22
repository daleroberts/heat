/* Solve 1D Heat Equation Using Fully Implicit Scheme
 * 
 * Dale Roberts <dale.o.roberts@gmail.com>
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <Accelerate/Accelerate.h>

#define minX -5
#define maxX 5
#define minT 0
#define maxT 2
#define J 10
#define M 100

double dl[J];
double d[J+1];
double du[J];
double B[J+1];

void print_vec(double *vec, int n)
{
  for (int i = 0; i < n; i++) {
    if (i%11 == 0 && i != 0)
      putchar('\n');

    printf("%8.2lf", vec[i]);
  }
    putchar('\n');
}

void lapack_dgtsv(long N, long nRHS, double *dl, double *d, double *du, double *b, long ldb)
{
  long info;
  dgtsv_(&N, &nRHS, dl, d, du, b, &ldb, &info);
  assert( info == 0 );
}

void alloc_vec(double *vec, int n, double val)
{
  vec = malloc(n * sizeof(double));
  assert(vec != NULL);
  
  for (int i = 0; i < n; i++)
    vec[i] = val;
}

void init_vec(double *vec, int n, double val)
{
  for (int i = 0; i < n; i++)
    vec[i] = val;
}

void set_initial(double *vec, double dx, int n)
{
  double x;

  for (int i = 0; i < n; i++) {
    x = minX + i*dx;
    vec[i] = exp(-x*x);
  }
}

int main() 
{
  double dx = (maxX - minX) / (double)J;
  double dt = (maxT - minT) / (double)M;
  double mu = dt / (dx*dx);

  init_vec(dl, J, -mu);
  init_vec(d, J+1, (1+2*mu));
  init_vec(du, J, -mu);
  init_vec(B, J+1, 0);

  set_initial(B, dx, J+1);

  
  print_vec(B, J+1);

  for (int m = 0; m < M; m++) {
    lapack_dgtsv(J+1, 1, dl, d, du, B, J+1);
    print_vec(B, J+1);
  }

}
