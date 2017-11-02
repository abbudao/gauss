#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include "minunit.h"
#include "gauss_func.c"

int tests_run = 0;

static char *test_intializeMatrix() {
double *A,order=10;

initializeMatrix(&A, order);
for(int i=0;i<10;i++){
  for(int j=0;j< 10;j++) {
   A[i] =0;
  }
}
for(int i=0;i<10;i++){
  for(int j=0;j< 10;j++) {
   A[i] =0;
   mu_assert("intializeMatrix failed test: Couldn't assign value to new matrix", A[i]== 0);
  }
}

return 0;
}

static char * all_tests() {
     mu_run_test(test_intializeMatrix);
     return 0;
}

int main(int argc, char **argv) {
     char *result = all_tests();
     if (result != 0) {
         printf("%s\n", result);
     }
     else {
         printf("ALL TESTS PASSED\n");
     }
     printf("Tests run: %d\n", tests_run);
 
     return result != 0;
 }
