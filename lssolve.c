#include <stdlib.h>
#include <stdio.h>
#include "msptools.h"

#include "array.h" // --> array_from_file(), array_to_file()
#include "array2d.h" // --> array2d_from_file()

int call_dgels(array2d_t * A, array_t * b, double * resnorm, double * rsquared);

/* prototypes 
int array_to_file(const char *filename, const array_t *a);
array_t *array_from_file(const char *filename);
array2d_t *array2d_from_file(const char *filename);
*/

int main(int argc, char *argv[])
{

  if (argc != 4)
  {
    fprintf(stderr, "Usage: %s A b x\n", argv[0]);
    return EXIT_FAILURE;
  }

  // Read the matrix from the file
  array2d_t *A = array2d_from_file(argv[1]);
  
  if (A==NULL){
    // Handle error:print message to stdrer and exit
    fprintf(stderr, "Error reading matrix A from file %s.\n", argv[1]);
    return EXIT_FAILURE; // No need to free A here, as it's NULL
  }

  // Error Handling
  if (A->shape[0]<A->shape[1]){
    fprintf(stderr,"Error: The matrix A has more columns than rows (m < n), which is not supported.\n");
    free(A->val);
    free(A);
    return EXIT_FAILURE;
  }


  // sos kati gnt edo 
  // A->order = ColMajor; //assume colmajor

  // Read the vector b from the file
  array_t *b = array_from_file(argv[2]);
  if (b == NULL){
    // Handle error:print message to stdrer and exit
    fprintf(stderr, "Error reading vector b from file %s.\n", argv[2]);
    // Free the previously allocated memory
    free(A->val); // Free the memory allocated for A's values
    free(A);
    return EXIT_FAILURE;
  }

  // Process Data
  double resnorm, rsquared;
  
  int result = call_dgels(A, b, &resnorm, &rsquared); //returns 0 if successful
  if (result!=0) {
    switch (result){
      case -12:
          fprintf(stderr, "Error: One of the input pointers is NULL.\n");
          break;
      case -13:
          fprintf(stderr, "Error: The matrix A has fewer rows than columns.\n");
          break;
      case -14:
          fprintf(stderr, "Error: Incompatible dimensions of A and b.\n");
          break;
      case -15:
          fprintf(stderr, "Error: Memory allocation failer.\n");
          break;
      // case -16:
      //     fprintf(stderr,"Error: The matrix A is rank-deficient and cannot be processed.\n");
      //     break;
      default:
          fprintf(stderr, "Error: DGELS returned non-zero info value.\n");
          break;
    }
    free(A);
    free(A->val);
    free(b);
    free(b->val);
    return EXIT_FAILURE;
  }

  // fprintf(stderr, "DGELS failed with info = %d.\n", result);
   // If successful, print the residual norm and the coefficient of determination

  
  printf("Residual norm: %f\n", resnorm);
  printf("Coefficient of determination: %f\n", rsquared);


  // write solution to text
  result = array_to_file(argv[3], b);
  if (result != 0) {
    fprintf(stderr, "Error writing the solution to file %s.\n", argv[3]);
    // No need to free the contents of b as array_to_file() should not modify it on failure
  }else {
    printf("The solution has been written to %s.\n", argv[3]); // Confirm the success to the user
  }

  // Free the allocated memory
  free(A->val);
  free(A);
  free(b->val);
  free(b);
  // this is unconditional --> return EXIT_SUCCESS;
  return result == 0? EXIT_SUCCESS : EXIT_FAILURE;
}
