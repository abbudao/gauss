void printMatrix(double *A, int order){
  int i, j;
  for(i=0;i<order;i++){
    for(j=0;j<order;j++){
      printf("%lf ", A[j*order+i]);
    }
    printf("\n" );
  }
  printf("-------------------\n");
}

