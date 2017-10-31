#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>

int findPivot(double *array, int len, int order, int row, int col, int *index){
/*
  @param array, array of the process
  @param len, length of the array
  @param order, order of the matrix
  @param row, begining index of the array or the pivot's row
  @param col, the collumn of the pivot
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
  int i, iaux, stride;
  double min;
  stride = col*order;
  for(i=row;i<order;i++){
    if(array[i+stride]!=0){
      min = array[i+stride];
      *index = i+stride;
      iaux = i+stride;
      break;
    }
  }
  //if not found exits with error
  if(array[iaux]==0){
      return 1;
  }
  //else, checks for the minimum
  for (i=iaux;i<order+stride;i++){
    if(fabs(array[i])<fabs(min) && array[i]!=0){
      *index = i;
      min = array[i];
  //    printf("minimo = %lf\n", min);
    }
  }
  return 0;
}

void changeRow(double **array, int len, int order, int row, int col,
  int index, int master, int myrank, int np, MPI_Comm comm){
/*
  @param array, array of the process
  @param len, the length of the array
  @param order, the order of the matrix
  @param row, the position where there will be the pivot
  @param index, the index of the row to be changed
  @param master, the process with Pivot
  @param myrank, rank of the process
  @param comm, comunicator

  Master sends the index of the row that will be changed.
  All the processes changes the two rows. If the process
  has more than one collum, it also works using len and
  order correctly.

  ------------------
  (test case old)
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
//  printf("%d\n", master)  ;
  //MPI_Comm_rank(comm, &myrank);
  MPI_Bcast(&index, 1, MPI_INT, master, comm);
  if(myrank>=master){
    int  i, stride;
    double aux;
    //index = 2;
    for(i=col;i<len/order;i++){
  //    printf("iiiiiiii = %d\n", i);
      stride = i*order;
    //  printf("stride = %d\n", stride);
      aux = (*array)[row+stride];
      (*array)[row+stride] = (*array)[index];
    //  printf("[index+stride] %lf\n", (*array)[index]);
      (*array)[index] = aux;
    //  printf("array = %lf \n", aux);
    }
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
  int stride = col*order;
  (*div)[row] = array[stride+row];
  for(i=0;i<order;i++){
    if(i!=row)
      (*div)[i]=array[stride+i]/array[stride+row];
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

void subtract(double **array, double *div, int len, int order, int np ,int row,
  int col, int myrank, int master)
{
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
  if(myrank>=master){
    int i, stride, k;
    double aux;
    for(i=col;i<order/np;i++){
    //  printf("subtracttttttttttt %d\n", i );
      stride = i*order;
      for (k=0;k<order;k++){
        if(k!=row){
          aux = (*array)[stride+row]*div[k];
          (*array)[stride+k] = (*array)[stride+k] - aux;
        }
      }
    //  printf("row-----=%lf\n", div[row]);
      (*array)[stride+row] = (*array)[stride+row]/div[row];
    }
  //  if(myrank == master){
    //  (*array)[row+col*order] = 1;
  //  }

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
//order=4, np =2, order/np =2
  MPI_Comm newCom;
  if(order/np==1||(*row)!=0 && ((*row)+1)%(order/np)==0){
    int color = 0;
    if(myrank==*master){
      color = MPI_UNDEFINED;
    //  color = 1;
      //printf("here\n" );
    }
  ///  MPI_Comm_split(*comm, color, myrank, &newCom);
    *master = *master+1;
    *col = 0;
    (*row)++;
  //  *comm = newCom;
    int nRank;
  //  MPI_Comm_rank(*comm, &nRank);
  //  printf("old %d  new %d\n", myrank, nRank);
  }
  else{
    if(*row==order-1){
      *master=*master+1;
    }
    (*row)++;
    (*col)++;

  }

//  int i;
}

void initializeMatrix(double **A, int order){
  int i;
  *A = (double*)malloc(order*order*sizeof(double));
  for(i=0;i<order*order;i++)
      (*A)[i] = rand()%50;
}

void printMatrix(double *A, int order){
  int i, j;
  //printf("-------------------\n");
  for(i=0;i<order;i++){
    for(j=0;j<order;j++){
      printf("%lf ", A[j*order+i]);
    }
    printf("\n" );
  }
  printf("-------------------\n");
}

void separateMatrix(double **array, double *matrix, int order, int np, int myrank){
    /* case order can't be divided by number of processes */
    if(order%np!=0){
      int count;
      int nDiv = np - (order%np);
      int n = (nDiv+order)/np;
      int nNormal = order/n;
      int nRest = order%n;
      printf("nDiv %d n %d nNOrmal %d  nRest %d \n",nDiv, n, nNormal, nRest );
      if(myrank==nNormal){
        *array = (double*)malloc(order*nRest*sizeof(double));
        count = order*nRest;
      }
      else if(myrank<nNormal){
        *array = (double*)malloc(order*n*sizeof(double));
        count = order*n;
      }
      else{
        *array = (double*)malloc(1*sizeof(double));
        count = 0;
      }
      int *displs = (int *)calloc(np,sizeof(int));
      int *elements = (int *)calloc(np,sizeof(int));
      int i;
      for (i=0; i<nNormal; ++i) {
        displs[i] = i*order*n;
        elements[i] = order*n;
    //    printf("elements = %d\n", elements[0]);
      }
      displs[nNormal] = (nNormal)*order*n;
      elements[nNormal] = order*nRest;
      printf("elem %d order = %d\n", elements[0], order*n);

      MPI_Scatterv(matrix, elements, displs, MPI_DOUBLE, *array, count, MPI_DOUBLE,
                                                        0, MPI_COMM_WORLD);

    //int i;
//    for(i=0;i<count;i++){
//      printf("%lf ", (*array)[i]);
//    }
//    printf("\n" );
    }
    /* normal case */
    else{
    *array = (double*)malloc(order*order/np*sizeof(double));
    MPI_Scatter (matrix, order*order/np, MPI_DOUBLE, *array, order*order/np,
      MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
}

void getMatrixTogether(double *array, double **A, int order, int myrank,
  int np, MPI_Comm comm){

  double *matrix = (double*)calloc(order*order,sizeof(double));
  if(order%np!=0){
    int count;
    int nDiv = np - (order%np);
    int n = (nDiv+order)/np;
    int nNormal = order/n;
    int nRest = order%n;
    if(myrank==nNormal){
      count = nRest;
    }
    else if(myrank<nNormal){
      count = n;
    }
    else{
      count = 0;
    }
    int *displs = (int *)calloc(np,sizeof(int));
    int *elements = (int *)calloc(np,sizeof(int));
    int i;
    for (i=0; i<np  ; ++i) {
        displs[i] = i*order*n;
        //if(count==0)
          //displs[i] = 0;
        if(i<nNormal)
          elements[i] = n*order;
        if(i==nNormal)
          elements[i] = nRest*order;
    }
    MPI_Gatherv(array, count*order, MPI_DOUBLE, matrix, elements, displs, MPI_DOUBLE,
                                                               0, MPI_COMM_WORLD);
  }
  else{
    MPI_Gather (array, order*order/np, MPI_DOUBLE, matrix, order*order/np, MPI_DOUBLE,
    0, MPI_COMM_WORLD);
  }
  *A = matrix;
}


int main(int argc, char **argv){

  MPI_Init(&argc, &argv);
  /* order is for n rows, not for collumns */
  int master, row, col, myrank, np, index, i, order;
  double *A, *array;
  order = 9;
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


  MPI_Finalize();
 return 0;
}
