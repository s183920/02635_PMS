CC = gcc
LDLIBS = -lopenblas
CFLAGS=-Wall -std=c99 -g
LDFLAGS=-L.

ifeq ($(shell uname), Darwin)
	# Link against system default BLAS/LAPACK library on macOS
	LDLIBS=-llapack -lblas
endif

lssolve: matrix_io.o call_dgesv.o solve.o

.PHONY: run
run: call_dgesv clean

matrix_io.o: matrix_io.c

call_dgesv: matrix_io.o

.PHONY: clean
clean:
	-rm *.o
