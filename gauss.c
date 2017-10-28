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
  @array, the array of the process
  @len, the length of the array
  @row, the index of the row where the pivot is
  @col, the index of the collum where the pivot is
  @div, the divisors of each row

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



int main(int argc, char **argv){



  MPI_Finalize();
 return 0;
}
