#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include "gauss_func.c"
#include "helpers.c"

int main(int argc, char **argv){

  MPI_Init(&argc, &argv);
  /* order is for n rows, not for collumns */
  int master, row, col, myrank, np, index, i, order;
  double *A, *array;
  order = 8;
  master = 0;
  col = 0;
  double *div = (double*)calloc(order,sizeof(double));
  row = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm comm;
  comm = MPI_COMM_WORLD;

  if (myrank==0){
    initializeMatrix(&A, order);
    printMatrix(A, order);
  }else{
    A = (double*)calloc(order, sizeof(double));
  }
  separateMatrix(&array, A, order, np, myrank);

while(master<np){

  if(myrank==master){
    findPivot(array, order*order/np ,order, row, col, &index);
  }
  changeRow(&array, order*order/np, order, row, col,index, master, myrank, np,comm);
  if(myrank==master)
    findDivisors(array, order*order/np, order, row, col, &div);
  sendDivisors(&div, order, master, myrank, comm);
  subtract(&array, div, order*order/np, order, np, row, col, myrank, master);
  getMatrixTogether(array, &A, order, myrank, np, comm);
  nextPivot(&master, &row, &col, myrank, order, np, &comm);
}
if(myrank==0)
  printMatrix(A, order);

/*while test*/
/*    if(myrank==master){
    findPivot(array, order*order/np ,order, row, col, &index);
  }
  if (myrank == master)
    printf("index = %d\n", index);
  changeRow(&array, order*order/np, order, row, col,index, master, myrank, np,comm);
  if(myrank==master)
    findDivisors(array, order*order/np, order, row, col, &div);
  sendDivisors(&div, order, master, myrank, comm);
  if(myrank==master){
    for(i=0;i<order;i++)
      printf("%lf ", div[i]);
  }
  printf("\n" );
  subtract(&array, div, order*order/np, order, np, row, col, myrank, master);
  getMatrixTogether(array, &A, order, myrank, np, comm);
  nextPivot(&master, &row, &col, myrank, order, np, &comm);
  nextPivot(&master, &row, &col, myrank, order, np, &comm);
//}
if(myrank==0)
  printMatrix(A, order);

  /*while test2*/
/*    if(myrank==master){
      findPivot(array, order*order/np ,order, row, col, &index);
    }
    if (myrank == master)
      printf("index = %d\n", index);
    changeRow(&array, order*order/np, order, row, col,index, master, myrank, np,comm);
    if(myrank==master)
      findDivisors(array, order*order/np, order, row, col, &div);
    sendDivisors(&div, order, master, myrank, comm);
    if(myrank==master){
      for(i=0;i<order;i++)
        printf("%lf ", div[i]);
    }
    printf("\n" );
    subtract(&array, div, order*order/np, order, np, row, col, myrank, master);
    getMatrixTogether(array, &A, order, myrank, np, comm);
    nextPivot(&master, &row, &col, myrank, order, np, &comm);
  //}

  /*while test*/
  /*  if(myrank==master){
      findPivot(array, order*order/np ,order, row, col, &index);
    }
    if (myrank == master)
      printf("index = %d\n", index);
    changeRow(&array, order*order/np, order, row, col,index, master, myrank, np,comm);
    if(myrank==master)
      findDivisors(array, order*order/np, order, row, col, &div);
    sendDivisors(&div, order, master, myrank, comm);
//    if(myrank==master){
//      for(i=0;i<order;i++)
//        printf("%lf ", div[i]);
//    }
    printf("\n" );
    subtract(&array, div, order*order/np, order, np, row, col, myrank, master);
    getMatrixTogether(array, &A, order, myrank, np, comm);
//    nextPivot(&master, &row, &col, myrank, order, np, &comm);
  //}
/*  if(myrank==0)
    printMatrix(A, order);
*/
  if(myrank==0)
    printMatrix(A, order);


  MPI_Finalize();
 return 0;
}
