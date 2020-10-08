#include <stdlib.h>
#include <math.h>
int fwdsub(
  unsigned long n,
  double alpha,
  double **R,  /* two-dimensional array, row-major */
  double *b    /* one-dimensional array */
  ){
  for(int k = 0; k < n; k++){
    double sum = 0;
    for(int i = 0; i < k; i++){
      sum += b[i]*R[i][k];
    }
    double div = (alpha+R[k][k]);
    if (div == 0){return -1;}
    b[k] = (b[k]-sum)/div;
  }

  return 0;
}
