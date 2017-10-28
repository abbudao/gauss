#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int findPivot(double *array, int len, int begin,int *index){
/*
  @param array, array of the process
  @param len, length of the array
  @param begin, begining index of the array
  @param index, index of the pivot

  Finds the row with the Pivot from an array.
  Returns by parameters index of the array with the
  Pivot.
  Returns 0 if sucessful and 1 if the array is full of
  zeros
  -----------------
  (test case)
  double A[10]={0,2,3,4,5,6,7,8,9,10};
  int index;
  findPivot(A, 10, 3, &index);
  printf("%lf\n", index);
  ----------------
*/
  //finds the first not zero
  int i, min, iaux;
  for(i=begin;i<len;i++){
    if(array[i]!=0){
      min = array[i];
      *index = i;
      iaux = i;
      break;
    }
  }
  //if not found exits with error
  if(array[iaux]==0){
  }
    return 1;
  //else, checks for the minimum
  for (i=iaux;i<len;i++){
    if(array[i]<min && array[i]!=0){
      *index = i;
      min = array[i];
    }
  }
  return 0;
}

void changeRow(double **array, int len, int order, int pivot,
  int row, int master, int myrank, MPI_Comm comm){
/*
  @param array, array of the process
  @param len, the length of the array
  @param order, the order of the matrix
  @param pivot, the position where there will be the pivot
  @param row, the index of the row to be changed
  @param master, the process with Pivot
  @param myrank, rank of the process
  @param comm, comunicator

  Master sends the index of the row that will be changed.
  All the processes changes the two rows. If the process
  has more than one collum, it also works using len and
  order correctly.

  ------------------
  (test case)
  double *A;
  A = (double*)calloc(4, sizeof(double));
  int index;
  int myrank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if(myrank == 0){
    A[0]=0;
    A[1]=1;
    A[2]=3;
    A[3]=2;
    A[4]=15;
    A[5]=14;
    findPivot(A, 6, 0, &index); //remember this condition
  }
  else{
    A[0]=5;
    A[1]=3;
    A[2]=6;
    A[3]=7;
    A[4]=9;
    A[5]=11;
  }
  changeRow(&A,6,3,0,index,0, myrank, MPI_COMM_WORLD);
  int i;
  for(i=0;i<6;i++){
    printf("%lf and %d\n", A[i], myrank);
  }
   -----------------------
*/

  MPI_Bcast(&row,1, MPI_INT,master, comm);
  int aux, i, stride;

  for(i=0;i<len/order-pivot/order;i++){
    stride = i*order;
    aux = (*array)[pivot+stride];
    (*array)[pivot+stride] = (*array)[row+stride];
    (*array)[row+stride] = aux;
  }

}

void findDivisors(double *array, int len, int order, int row,
  int col, double **div){
/*
  @array, array of the process
  @len, length of the array
  @order, order of the matrix
  @row, index of the row where the pivot is
  @col, index of the collum where the pivot is
  @div, divisors of each row

  Finds the divisors of the pivot's collum and returns them as param div.

  -------------------
  (test case)
  double *A;
  A = (double*)calloc(4, sizeof(double));
  int index;
  int myrank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if(myrank == 0){
    int pivotRow = 1;
    int pivotCol = 0;
    int order = 2;
    double *div;
    div = (double*)calloc(order-pivotRow, sizeof(double));
    A[0]=0;
    A[1]=1;
    A[2]=3;
    A[3]=2;
    A[4]=15;
    A[5]=14;
    findDivisors(A, 6, order, pivotRow, pivotCol, &div);
    int i;
    for(i=0;i<order-pivotRow;i++){
      printf("%lf\n", div[i]);

    }
  }
  else{
    A[0]=5;
    A[1]=3;
    A[2]=6;
    A[3]=7;
    A[4]=9;
    A[5]=11;
  }
  -------------------
*/
  int i;
  for(i=0;i<len/order-row;i++){
    int stride = col*len/order;
    (*div)[i]=array[stride+row+1+i]/array[stride+row];
  }
}

void sendDivisors(double **div, int len, int master, int myrank, MPI_Comm comm){
/*
  @param div, array of divisors to be sent
  @param len, length of the array
  @param master, number of the process sending
  @param myrank, number of the process
  @param comm, communicator

  Sends the divisors to all process
  --------------
  (test case)
  MPI_Init(&argc, &argv);
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  double *div;
  div = (double*)calloc(6, sizeof(double));
  div[0] = 1;
  div[1] = 2;
  div[2] = 3;
  div[3] = 4;
  div[4] = 5;
  div[5] = 6;
  sendDivisors(&div, 6, 0, myrank, MPI_COMM_WORLD);
  int i;
  for(i=0;i<6;i++){
    printf("%lf\n", div[i]);
  }
  --------------
*/
  double *buf = (double*)calloc(len, sizeof(double));
  if (myrank==master){
    buf = *div;
    }
  MPI_Bcast(buf, len, MPI_DOUBLE , master, comm);
  *div = buf;
}

void subtract(double **array, double *div, int len, int order, int row){
/*
  @param array, array of the process
  @param div, divisors from the pivot process
  @param len, length of the array
  @param order, order of the matrix
  @param row, index of the pivot's row

  Multiply the row of the pivot and subtracts it in the collunm
  ----------------
  (test case)
  MPI_Init(&argc, &argv);
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  double *div;
  div = (double*)calloc(4, sizeof(double));
  double *A = (double*)calloc(6, sizeof(double));
  div[0] = 1;
  div[1] = 2;
  div[2] = 3;
  div[3] = 4;
  A[0]=5;
  A[1]=3;
  A[2]=6;
  A[3]=7;
  A[4]=9;
  A[5]=11;
  subtract(&A, div, 6, 2, 0);
  int i;
  for(i=0;i<6;i++){
    printf("%lf\n", A[i]);
  }
  MPI_Finalize();
 return 0;
 ----------------
*/
  int i, stride, k;
  double aux;
  for(i=0;i<order;i++){
    stride = i*len/order;
    for (k=0;k<(len/order-row-1);k++){
      aux = (*array)[stride+row]*div[stride+k-i*(row+1)];
      (*array)[stride+row+1+k] = (*array)[stride+row+k+1] - aux;
    }
  }
}

void nextPivot(int *master, int *row, int *col, int myrank, int order, int np,
  MPI_Comm *comm){
/*
  Pass the responsability to the next pivot because the last one has
  finished its job or because there are no other rows different of zero in
  that collumn

  @param master, rannk of the process with pivot
  @param row, row of the pivot in the process
  @param col, col of the pivot in the process
  @param myrank, rank of the process
  @param order, order of the matrix
  @param np, number of processes
  @param comm, comunicator
  ---------------
  (test case)
  MPI_Init(&argc, &argv);
  int master,row,col;
  master = row = col = 0;
  int myrank;
  MPI_Comm comm;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  comm = MPI_COMM_WORLD;
  int test= 0;
  MPI_Comm newC;
  nextPivot(&master, &row, &col, myrank, 4, 2, &comm);
  printf("myrank = %d master = %d row = %d col = %d\n", myrank, master, row, col);
  //test++;
  //MPI_Comm newC2;
  nextPivot(&master, &row, &col, myrank, 4, 2, &comm);
  printf("myrank = %d master = %d row = %d col = %d\n",myrank, master, row, col);
  nextPivot(&master, &row, &col, myrank, 4, 2, &comm);
  printf("myrank = %d master = %d row = %d col = %d\n",myrank, master, row, col);

  MPI_Finalize();
 return 0;
 ----------------
*/
  if(order/np==1||*row!=0 && *row%order/np==0){
    int color = 0;
    if(myrank==*master){
      color = MPI_UNDEFINED;
    }
    MPI_Comm_split(*comm, color, myrank, comm);
    *master = *master+1;
    *row = *col = 0;
  }
  else{
    (*row)++;
    (*col)++;
  }

}


int main(int argc, char **argv){
  MPI_Init(&argc, &argv);
  int master,row,col;
  master = row = col = 0;
  int myrank;
  MPI_Comm comm;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  comm = MPI_COMM_WORLD;
  
  MPI_Finalize();
 return 0;
}
