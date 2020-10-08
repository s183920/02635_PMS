#ifndef MATRIX_H
#define MATRIX_H

/* Structure representing a vector */
typedef struct vector {
  unsigned long n;   /* length of vector                     */
  double * v;        /* pointer to array of length n         */
} vector_t;

/* Structure representing a matrix */
typedef struct matrix {
  unsigned long m;   /* number of rows                       */
  unsigned long n;   /* number of columns                    */
  double ** A;       /* pointer to two-dimensional array     */
} matrix_t;

#endif
