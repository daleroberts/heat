/* 2D Heat equation using the Peaceman-Rachford algorithm
 *
 * Dale Roberts <dale.o.roberts@gmail.com>
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define R 64

#define MAX(A, B) A < B? B: A

/*
 * Set initial condition where m = rows, n = cols, dr = step
 */
void setInitialCondition(double *U, int m, int n, double dr) {
  double x, y;
  int r, c;

  for (r = 0; r < m; r++) {
    y = r * dr;
    for (c = 0; c < n; c++) {
      x = c * dr;
      /*      U[r*n+c] = MAX(1.0/3.0 - sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)), 0.0); */
      U[r*n+c] = MAX(0.01, 1.0 - sqrt(20*(x-0.1)*20*(x-0.1)+20*(y-0.1)*20*(y-0.1)));
      U[r*n+c] = MAX(U[r*n+c], 1.0 - sqrt(20*(x-0.9)*20*(x-0.9)+20*(y-0.9)*20*(y-0.9)));
    }
  }
}

/*
 * Force boundary conditions
 */
void forceBoundaryConditions(double *U, int m, int n, double x) {
  catlas_dset((R+1), x, U, 1);
  catlas_dset((R+1), x, &U[R*(R+1)], 1);
  catlas_dset((R+1), x, U, (R+1));
  catlas_dset((R+1), x, &U[R], (R+1));
}

/*
 * Construct RHSs - x sweep
 */
void setRHSx(double *U, double *b, int m, int n, double mu) {
  int r, c;

  for (r = 0; r < m; r++) {
    b[r*n] = 0.0;
    for (c=1; c<(n-1); c++)
      b[r*n+c] = mu*U[r*n+(c-1)] + (1-2*mu)*U[r*n+c] + mu*U[r*n+(c+1)];
    b[r*n+(n-1)] = 1.0;
  }
}

/*
 * Construct RHSs - y sweep
 */
void setRHSy(double *U, double *b, int m, int n, double mu) {
  int r, c;

  for (c = 0; c < n; c++) {
    b[c] = 0.0;
    for (r=1; r<(m-1); r++)
      b[r*n+c] = mu*U[(r-1)*n+c] + (1-2*mu)*U[r*n+c] + mu*U[(r+1)*n+c];
    b[(m-1)*n+c] = 1.0;
  }
}

void setA(double *du, double *dc, double *dl, int m, int n, double mu) {
  catlas_dset(m-1, -mu, du, 1);
  catlas_dset(m, (1+2*mu), dc, 1);
  catlas_dset(m-1, -mu, dl, 1);

  dc[0] = 1;
  du[0] = 0;
  
  dc[m-1] = 1;
  dl[m-2] = -1;
}

void printMatrix(double *A, int m, int n) {
  int r, c;

  puts("");

  for (r = 0; r < m; r++) {
    for (c = 0; c < n; c++)
      printf("%7.3lf", A[c*m+r]);
    putchar('\n');
  }
}

void printVector(double *v, int m) {
  int r;

  puts("");

  for (r = 0; r < m; r++)
      printf("%6.2lf", v[r]);

  putchar('\n');
}

void lapack_dgtsv(long dim, long nRHS, double *dl, double *d, double *du, double *b, long ldb) {
  long info;
  dgtsv_(&dim, &nRHS, dl, d, du, b, &ldb, &info);

  if (info != 0)
    printf("dgtsv = %li\n", info);
}

void diffuse(double *U, double *b, int m, int n, double *du, double *dc, double *dl, double mu) {
  setRHSx(U, b, (R+1), (R+1), mu);
  lapack_dgtsv((R+1), (R+1), dl, dc, du, b, (R+1));

  cblas_dcopy((R+1)*(R+1), b, 1, U, 1);

  setRHSy(U, b, (R+1), (R+1), mu);
  lapack_dgtsv((R+1), (R+1), dl, dc, du, b, (R+1));

  cblas_dcopy((R+1)*(R+1), b, 1, U, 1);

/*   forceBoundaryConditions(U, (R+1), (R+1), 0.0); */
}

void dumpMatrix(double *A, int m, int n, double dr, int no) {
  double x, y;
  int r, c;
  char fname[20];
  FILE *fd;

  snprintf(fname, sizeof(fname), "data.%i", no);
  fd = fopen(fname, "w");

  for (c = 0; c < n; c++) {
    x = c * dr;
    for (r = 0; r < m; r++) {
      y = r * dr;
      fprintf(fd, "%10.5lf %10.5lf %10.5lf\n", x, y, A[r*n+c]);
    }
  }

  fclose(fd);
}

int main() {
  double dr = 1.0 / (double) R;
  double dt = 0.0001;
  double mu = 0.5 * dt / (dr*dr);
  int i = 0;
  
  double U[(R+1)*(R+1)];
  double b[(R+1)*(R+1)];

  double du[R];
  double dc[R+1];
  double dl[R];

  printf("mu = %.5lf\n", mu);
  
  setInitialCondition(U, (R+1), (R+1), dr);
  setA(du, dc, dl, (R+1), (R+1), mu);
  dumpMatrix(U, (R+1), (R+1), dr, 0);

  while (i++ < 5) {
    diffuse(U, b, (R+1), (R+1), du, dc, dl, mu);
    dumpMatrix(U, (R+1), (R+1), dr, i);
  }
  
}
