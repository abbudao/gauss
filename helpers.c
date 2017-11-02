void printMatrix(double *A, int order){
  int i, j;
  for(i=0;i<order;i++){
    for(j=0;j<order;j++){
      if(j==0){
        printf("|");
      }
      printf(" %-6.1f ", A[j*order+i]);
      if(j==order-1){
        printf("|");
      }
    }
      printf("\n");
  }
  for(i=0;i<order;i++){
    printf("----");
  }
      printf("\n");
}
