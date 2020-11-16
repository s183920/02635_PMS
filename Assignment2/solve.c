#include <stdlib.h>
#include <stdio.h>
#include "matrix_io.h"

int call_dgesv(matrix_t * A, vector_t * b);

int main(int argc, char *argv[]) {

  if (argc != 4) {
    fprintf(stderr,"Usage: %s A b x\n", argv[0]);
    return EXIT_FAILURE;  
  }

  // read data
  matrix_t *A = read_matrix(argv[1]);
  if (!A){return EXIT_FAILURE;}
  vector_t *b = read_vector(argv[2]);
  if (!b){return EXIT_FAILURE;}

  // solve system
  int info = call_dgesv(A, b);
  if (info < 0){
    fprintf(stderr, "The %d th value had an illegal value", info);
    return EXIT_FAILURE;
  } else if (info > 0){
    fprintf(stderr, "Solution could not be computed, as matrix was singular");
    return EXIT_FAILURE;
  }

  // write solution
  write_vector(argv[3], b);

  free_matrix(A); free_vector(b);

  return EXIT_SUCCESS;
}
