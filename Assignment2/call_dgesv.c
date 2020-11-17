#include "matrix_io.h"

/* C prototype for LAPACK routine DGESV */
void dgesv_(const int *n,    /* columns/rows in A          */
            const int *nrhs, /* number of right-hand sides */
            double *A,       /* array A                    */
            const int *lda,  /* leading dimension of A     */
            int *ipiv,       /* pivoting array             */
            double *B,       /* array B                    */
            const int *ldb,  /* leading dimension of B     */
            int *info        /* status code                */
            );

double* flatten(matrix_t *A);

/* call_dgesv : wrapper for LAPACK's DGESV routine

Purpose:
Solves system of equations A*x=b using LAPACK's DGESV routine
Upon exit, the input vector b is overwriten by the solution x.

Return value:
The function returns the output `info` from DGESV with the
following exceptions: the return value is

   -9 if the input A is NULL and/or the input B is NULL
   -10 if A is not a square matrix 
   -11 if the dimensions of A and b are incompatible
   -12 in case of memory allocation errors.
*/
int call_dgesv(matrix_t * A, vector_t * b) {
  // Check input
  if (A == NULL || b == NULL){fprintf(stderr, "Input matrix was NULL");return -9;}
  if (A->n != A->m){fprintf(stderr, "A is not a quadratic matrix");return -10;}
  if (A->n != b->n){fprintf(stderr, "Leading dimnesions of A and b did not match");return -11;}

  // Allocate memory
  int n = A->n;
  int nrhs = 1; /* number of right-hand sides */
  int lda = A->n;  /* leading dimension of A     */
  int *ipiv;      /* pivoting array             */
  ipiv = malloc(sizeof(int)*A->n);
  int ldb = b->n;  /* leading dimension of B     */
  int info;
  double *A_flat = flatten(A);

  // Check for allocation error
  if (ipiv == NULL){
    fprintf(stderr, "Allocation failed");
    free(ipiv);
    return -12;
  }
  if (A_flat == NULL){
    fprintf(stderr, "Allocation failed");
    return -12;
  }

  // call the function
  dgesv_(&n, &nrhs, A_flat, &lda, ipiv, b->v, &ldb, &info);
 
  // free memory
  free(ipiv); 
  free(A_flat);

  // return result
  return info;
}


double* flatten(matrix_t *A){
  // Function to flatten the matrix A to a vector A_flat, this is done columnwise as LAPACK uses column-major order
  double *A_flat = calloc(A->m*A->n, sizeof(double));
  if (A_flat == NULL){
    fprintf(stderr, "Allocation failed");
    free(A_flat);
    return NULL;
  }

  for (int i = 0; i < A->m; i++){
    for (int j = 0; j < A->n; j++){
      A_flat[i*A->n+j] = A->A[j][i];
    }
  }
  return A_flat;
}


// int main(void){
//   matrix_t *A = malloc_matrix(2, 2);
//   vector_t *b = malloc_vector(2);

//   A->A[0][0] = 1;
//   A->A[0][1] = 0;
//   A->A[1][0] = 0;
//   A->A[1][1] = 1;


//   b->v[0] = 2;
//   b->v[1] = 3;

//   print_matrix(A);
//   print_vector(b);

//   int info = call_dgesv(A, b);
  
//   printf("Outcome: %d\n", info);
//   // free_vector(b); free_matrix(A);

//   printf("Solution: %lf, %lf\n", b->v[0], b->v[1]);
//   free_matrix(A); free_vector(b);

//   return 0;
// }