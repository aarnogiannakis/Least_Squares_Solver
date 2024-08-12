#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "msptools.h"
#include <string.h>

/* C prototype for LAPACK routine DGELS */
void dgels_(const char *trans, /* 'N' or 'T'             */
	const int *m,	   /* rows in A              */
	const int *n,	   /* cols in A              */
	const int *nrhs,   /* cols in B              */
	double *A,		   /* array A                */
	const int *lda,	   /* leading dimension of A */
	double *B,		   /* array B                */
	const int *ldb,	   /* leading dimension of B */
	double *work,	   /* workspace array        */
	int *lwork,		   /* workspace size         */
	int *info		   /* status code            */
);

double dcemv(array_t * b); // prototype for centering matrix function

int call_dgels(array2d_t * A, array_t * b, double * resnorm, double * rsquared)
{
	// null pointers
	if (A==NULL || b==NULL){
		return -12;
	}

	int m = A->shape[0]; // number of rows of A
	int n = A->shape[1]; // number of col of A

	// if A has fewer rows than columns 
	if (m < n){ // it must m>=n)
		return -13;
	}

	// incompatible dimensions of A and b
	if (m != b->len){
		return -14;
	}

	#define min(a,b) ( (a) < (b) ? (a) : (b) )
	#define max(a,b) ( (a) > (b) ? (a) : (b) )
	int nrhs = 1; // b is a vector
	int lwork = max(1, min(m,n) + max(min(m,n),nrhs));
	double *work = (double *)malloc(lwork * sizeof(double));

	// in case of memory allocation erros
	if (work == NULL){
		return -15;
	}
	int info;

	// Make a copy of b for later use in rsquared calculation
	double *original_b = malloc(b->len*sizeof(double)); 
	if (!original_b){
		free(work);
		return -15;
	}
	memcpy(original_b, b->val, b->len*sizeof(double));
	array_t original_b_array = {.len = b->len, .val = original_b};
	double computed_resnorm = 0.0;

    if (A->order == RowMajor){
		int lda = A->shape[1]; // if you use shape[1] you fail test1
		dgels_("T", &n, &m, &nrhs, A->val, &lda, b->val, &m, work, &lwork, &info);


		if (info != 0){
			fprintf(stderr,"Error: The matrix A is rank-deficient at column %d.\n", info);
			free(work);
			free(original_b);
			return info; //-16; 
		}
		b->len = min(n,m); // will be n
		if (m > n) {
			double sum_of_squares = 0.0;
			for (size_t i = n; i < m; i++) {
            	sum_of_squares += b->val[i] * b->val[i];
        	}
			computed_resnorm = sqrt(sum_of_squares);
			// printf("sum %f\n", sum_of_squares);
			// printf("residual %f\n", *resnorm);
			// exactly determined
		} else if (m ==n){
			computed_resnorm = 0.0; // // No residuals for a square system
		}

    } else {
        // When A is column-major, use the original m and n
        dgels_("N", &m, &n, &nrhs, A->val, &m, b->val, &m, work, &lwork, &info); 
		// Immediately check and print the value of 'info'
		// printf("metal cjeck DGELS info: %d\n", info);

		if (info != 0){
			fprintf(stderr,"Error: The matrix A is rank-deficient at column %d.\n", info);
			free(work);
			free(original_b);
			return info; 
		}
	
		
    	if (m > n) { // mayve this is redundant, cause we know that m>n
			double sum_of_squares = 0.0;
			for (size_t i = n; i < m; i++) {
				sum_of_squares += b->val[i] * b->val[i];
			}
			computed_resnorm = sqrt(sum_of_squares);
		} else if (m ==n) { // the system has no residuals
		//  This shouldn't happen since m should be >= n for column-major
			computed_resnorm = 0.0;
		}
    } // end of A->order check

    // Calculate rsquared using the centered original_b
	// now this makes the problem b is now changed it is not robust, you need the b len of the original
    
	if (resnorm != NULL){
		*resnorm = computed_resnorm;
	}

	// If both are NULL, your function would effectively just be solving the linear system without providing the additional statistics.
    // Calculate rsquared
    if (rsquared != NULL) {
		double total_sum_of_squares = dcemv(&original_b_array);
        *rsquared = 1 - (computed_resnorm * computed_resnorm) / total_sum_of_squares;
    }

	// Adjust the length of b to match the solution size
		// The first min(m, n) entries of b are the solution
	b->len = min(n,m);

    // Free allocated memory
    free(work);
    free(original_b);
    return info; // Successful exit, not 0 
}

// calls the copy of the original b, which it then overwrites
double dcemv(array_t *original_b_array) {
    double sum = 0.0;
	double sum_of_squares = 0.0;
    double avg = 0.0;

    //Error checks
    if (original_b_array == NULL || original_b_array->len ==0 || original_b_array->val == NULL ) {
    	return -1.0; // Return a negative value to indicate an error
    }

    // Calculate the sum of the elements in b->val
    for (size_t i=0; i< original_b_array->len; i++) { //i<5;
		// printf("%f\n", original_b_array->val[i]);
		
    	sum += original_b_array->val[i];
		// printf("%f\n",sum);
    }
	// printf("%f", original_b_array->len);
    // Calculate the avg
    avg = sum / original_b_array->len;// /5 
	// printf("mean is %f\n",avg);

	// with 5 its correct i need to fix the len of origianl b

	// total sum of squares 44.140948 should 68.9

    // Overwrite the vector b->val with centered values and calculate the sum of squares
    for (size_t i=0; i<original_b_array->len; i++) { //dia 5
    	original_b_array->val[i] -= avg; // Center the value
		sum_of_squares += original_b_array->val[i] * original_b_array->val[i]; // Add the square of the centered value
    }
	// // Return the sum of squares of the centered vector
    return sum_of_squares;
}


