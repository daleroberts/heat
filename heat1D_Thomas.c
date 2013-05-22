/* Heat equation using the Thomas algorithm
 *
 * Dale Roberts <dale.o.roberts@gmail.com>
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N 100

double U[N+1];
double b[N+1];

double minX = 0.0;
double maxX = 1.0;

double du[N];
double dc[N+1];
double dl[N];

float dx, dt, mu, theta;

void dumpSolution(int n) {
  char fname[20];
  int s;
  double x;
  FILE *fd;
  
  snprintf(fname, sizeof(fname), "data.%i", n);
  fd = fopen(fname, "w");

  for (s = 0; s < N+1; s++) {
    x = minX + s * dx;

    fprintf(fd, "%10.5lf %10.5lf\n", x, U[s]);
  }
  
  fclose(fd);
}

void dumpRHS(int n) {
  char fname[20];
  int s;
  double x;
  FILE *fd;
  
  snprintf(fname, sizeof(fname), "rhs.%i", n);
  fd = fopen(fname, "w");

  for (s = 0; s < N+1; s++) {
    x = minX + s * dx;

    fprintf(fd, "%10.5lf\n", b[s]);
  }
  
  fclose(fd);
}

void lapack_dgtsv(long dim, long nRHS, double *dl, double *d, double *du, double *b, long ldb) {
  long info;
  dgtsv_(&dim, &nRHS, dl, d, du, b, &ldb, &info);

  if (info != 0)
    printf("sgtsv = %li\n", info);
}

void initialDensity(void) {
  int s;
  double x;
  
  for (s = 0; s < N+1; s++) {
    x = minX + s * dx;
    if (x <= 0.5)
      U[s] = 2*x;
    else
      U[s] = 2-2*x;
  }
}

void initDiagonals(void) {
  int s;

  for (s = 0; s < N; s++) {
    dc[s] = (1+2*theta*mu);
    du[s] = -0.5 * mu;
    dl[s] = -0.5 * mu;
  }
  
  /* boundary conditions */
  dc[0]   = 1.0;
  du[0]   = 0.0;

  dl[N-1] = 0.0;
  dc[N]   = 1.0;
}

void createRHS(void) {
  int s;

  for (s = 1; s < N; s++)
    b[s] = U[s] + (1-theta)*mu*(U[s-1] - 2 * U[s] + U[s+1]);

  /* boundary conditions */
  b[0] = 0.0;
  b[N] = 0.0;
}


void stepSolver(int t) {
  int s;

  createRHS();
  dumpRHS(t);
  lapack_dgtsv(N+1, 1, dl, dc, du, b, N+1);

  for (s = 0; s < N+1; s++)
    U[s] = b[s];
}

void initSolver(void) {
  dx = (maxX - minX) / (double) N;
  dt = 0.001;

  theta = 0.6;

  mu = dt / (dx*dx);

  printf("\ndx=%.5lf dt=%.5lf mu=%.5lf theta=%.5lf\n", dx, dt, mu, theta);

  if (0.0 <= theta && theta < 0.5 && mu > 0.5/(1-2*theta))
    puts("unstable parameters");
  else
    puts("stable parameters");

  initDiagonals();

  initialDensity();
}

int main(void) {
  int t;

  initSolver();
  dumpSolution(0);

  for (t = 1; t < 20; t++) {
    stepSolver(t);
    dumpSolution(t);
  }

  return 0;
}
