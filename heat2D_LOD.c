/* 2D Heat equation using the LOD method
 *
 * Dale Roberts <dale.o.roberts@gmail.com>
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "config.h"

extern float U[N+1][N+1];
extern float V[N+1][N+1];
extern float D[N+1][N+1];

float minX = -5.0;
float maxX =  5.0;
float minY = -5.0;
float maxY =  5.0;

float b[N+1];

float du[N];
float dc[N];
float dl[N];

float dx, dy, dt, mux, muy;

extern int showVectorField;


void print_vec(float *vec, int n)
{
  for (int i = 0; i < n; i++) {
    if (i%8 == 0 && i != 0)
      putchar('\n');

    printf("%10.2Ef", vec[i]);
  }
    putchar('\n');
}

void lapack_sgtsv(long dim, long nRHS, float *dl, float *d, float *du, float *b, long ldb)
{
  long info;
  sgtsv_(&dim, &nRHS, dl, d, du, b, &ldb, &info);
  assert( info == 0 );
}

inline float max(float a, float b) {
  return a<b? a: b;
}

int inBound(int coord) {
  if (coord < 0)
    return 0;
  else if (coord > N)
    return N;
  else
    return coord;
}

void initialDensity(void) {
  float x, y;
  
  for (int s = 0; s < N+1; s++) {
    x = minX + s * dx;
    for (int r = 0; r < N+1; r++) {
      y = minY + r * dy;

      D[s][r] = exp(-x*x-y*y);
    }
  }
  
}

void initialVectorField(void) {
  float x, y;

  for (int s = 0; s < N+1; s++) {
    x = minX + s * dx;
    for (int r = 0; r < N+1; r++) {
      y = minY + r * dy;
      
      U[s][r] = sin(x*y);
      V[s][r] = cos(x*y);
    }
  }
}


void initDiagonals(void) {
  dc[0] = 1;
  du[0] = 0;

  for (int s = 1; s < N; s++) {
    dc[s] = (1+mux);
    du[s] = -0.5 * mux;
    dl[s-1] = -0.5 * mux;
  }
  
  dc[N]   = 1;
  dl[N-1] = 0;
}


void stepSolver(void) {
  // Generate U^{n+*}
  for (int r = 0; r < N+1; r++) {
    b[0] = D[0][r];

    for (int s = 1; s < N; s++)
      b[s] = 0.5 * mux * D[s-1][r] + (1 - mux) * D[s][r] + 0.5 * mux * D[s+1][r];
   
    b[N] = D[N][r];

    lapack_sgtsv(N+1, 1, dl, dc, du, b, N+1);

    for (int s = 0; s < N+1; s++)
      D[s][r] = b[s];
  }

  // Generate U^n
  for (int s = 0; s < N+1; s++) {
    b[0] = D[s][0];

    for (int r = 1; r < N; r++)
      b[r] = 0.5 * mux * D[s][r-1] + (1 - mux) * D[s][r] + 0.5 * mux * D[s][r+1];
   
    b[N] = D[s][N];

    lapack_sgtsv(N+1, 1, dl, dc, du, b, N+1);

    for (int r = 0; r < N+1; r++)
      D[s][r] = b[r];
  }

}

void initSolver(void) {
  dx = 10.0 / (float)N;
  dy = 10.0 / (float)N;
  dt = 0.001;

  mux = dt / (dx*dx);
  muy = dt / (dy*dy);

  //printf("dx=%.3f dy=%.3f dt=%.3f mu_x=%.3f mu_y=%.3f\n", dx, dy, dt, mux, muy);

  initDiagonals();

  initialVectorField();
  initialDensity();

  showVectorField = 0;
}

