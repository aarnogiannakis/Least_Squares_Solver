# Least Squares Solver Using LAPACK's DGELS Routine

This project implements a command-line tool to solve the least-squares problem. Given an over-determined system of linear equations (more equations than unknowns), the least-squares method finds the best-fitting solution by minimizing the sum of the squares of the residuals (the differences between the observed and computed values). This tool leverages the LAPACK library's DGELS routine to efficiently compute the least-squares solution.

### Prerequisites

- GCC Compiler
- LAPACK and BLAS libraries

### Building the Project

To build the project, run the following commands:

```bash
make clean
make all

###
in your command line type ./lssolve Data/A1.txt Data/b1.txt Data/x1.txt


## Acknowledgements

This project was developed as part of [02365 Mathematical Software Programming] at [DTU]. The Makefile and initial project structure were provided by Professor Martin Skovgaard Andersen. All code in `call_dgels.c` and `lssolve.c` was written by me as part of the course assignment.
