#include <stdlib.h>
#include <math.h>
#include "matrix.h"

int fwdsub(unsigned long n, double alpha, double **R, double *b);
double* step1(vector_t c, const matrix_t *R, matrix_t *C, unsigned long k);

int tri_sylvester_solve(const matrix_t *R, matrix_t *C) {
  if (R == NULL || C == NULL || R->A == NULL || C->A == NULL){return -2;} // Check if input is null
  if (R->m != R->n || C->m != C->n || R->m != C->m){return -2;} // Check that dimensions match

  for (int k = 0; k < C->m; k++){
    vector_t c_k;
    c_k.n = C->n;
    c_k.v = (double*) calloc(c_k.n, sizeof(double));
    if (c_k.v == NULL){return -2;}
    
    // First step
    c_k.v = step1(c_k, R, C, k);
    // Second step
    if (fwdsub(C->m, R->A[k][k], R->A, c_k.v) == -1){return -1;} // use fwdsub where n is the number of rows, alpha is the diagonal element of R, R is the R matrix and b is the kth row of C

    for (int i = 0; i < c_k.n; i++){C->A[k][i] = c_k.v[i];}
    free(c_k.v);
  }
  
  return 0;
}

double* step1(vector_t c_k, const matrix_t *R, matrix_t * C, unsigned long k){
  for (int i = 0; i < c_k.n; i++){
      double sum = 0;
      c_k.v[i] = C->A[k][i];
      for (int j = 0; j < k; j++){
        sum += C->A[j][i] * R->A[j][k]; // cant use ck as we need to acces to columns instead of rows
      }
      c_k.v[i] -= sum;
    }
  return c_k.v;
}


int fwdsub(
  unsigned long n,
  double alpha,
  double **R,  /* two-dimensional array, row-major */
  double *b    /* one-dimensional array */
  ){
  for(int k = 0; k < n; k++){
    double div = (alpha+R[k][k]);
    if (div == 0){return -1;}

    double sum = 0;
    for(int i = 0; i < k; i++){
      sum += b[i]*R[i][k];
    }
    
    b[k] = (b[k]-sum)/div;
  }

  return 0;
}

