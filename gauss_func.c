#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>

void getProcess(int *type, int *size, int order ,int np){
/*
  Returns 0 for normal process
          1 for last process
          2 for the rest
*/
  if(order%np!=0){
    int myrank;
    int nDiv = np - (order%np);
    int n = (nDiv+order)/np;
    int nNormal = order/n;
    int nRest = order%n;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if(myrank==nNormal){
      (*type) = 1;
      (*size) = nRest;
    }
    else if(myrank<nNormal){
      (*type) = 0;
      (*size) = n;
    }
    else{
      (*type) = 2;
      (*size) = 0;
    }
  }else{
    (*type) = 0;
    (*size) = order/np;
  }
}

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
*/
  MPI_Bcast(&index, 1, MPI_INT, master, comm);
  if(myrank>=master){
    int  i, stride;
    double aux;
    int type, size;
    getProcess(&type, &size, order ,np);
    for(i=col;i<size;i++){
      stride = i*order;
      aux = (*array)[row+stride];
      (*array)[row+stride] = (*array)[(index%order)+stride];
      (*array)[(index%order)+stride] = aux;
      printf("aux = %lf array row+stride = %lf i = %d\n",aux,  (*array)[row+stride], i);
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
*/
  printf("function ---\n" );
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
*/
  if(myrank>=master){
    int i, stride, k;
    double aux;
    int type, size;
    getProcess(&type, &size, order ,np);
    for(i=col;i<size;i++){
      stride = i*order;
      for (k=0;k<order;k++){
        if(k!=row){
          aux = (*array)[stride+row]*div[k];
          (*array)[stride+k] = (*array)[stride+k] - aux;
        }
      }
      (*array)[stride+row] = (*array)[stride+row]/div[row];
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
*/
  MPI_Comm newCom;
  int type, size;
//  getProcess(&type, &size, order ,np);
  if(myrank ==0)
  printf("%d row = %d col = %d \n",*master,  *row, *col);
  int nDiv = np - (order%np);
  int n = (nDiv+order)/np;
  int nNormal = order/n;
  int nRest = order%n;
  if((order/np==1 && order%np==0) || ( (*row)!=0 && ((*col+1))%(n)==0 && *master<nNormal-2)
    ||((*row)!=0 && ((*col+1))%(nRest)==0 && *master==nNormal-2)){
    *master = *master+1;
    *col = 0;
    (*row)++;
//    printf("hereeeeee, %d  \n",myrank );
  }
  else{
  //  if(*row==n-1 && master<nNormal-2 || *row==nRest && Master){
  //    *master=*master+1;
  //  }
    (*row)++;
    (*col)++;
  }
  if(myrank==0)
  printf("master = %d row = %d col = %d\n", *master, *row, *col );
}

void initializeMatrix(double **A, int order){
/*
  @param A, matrix pointer to be intialized
  @param order, order of matrix
Intialize a random square matrix [order x order]
*/
  int i;
  *A = (double*)malloc(order*order*sizeof(double));
  for(i=0;i<order*order;i++)
      (*A)[i] = rand()%50;
}


void separateMatrix(double **array, double *matrix, int order, int np, int myrank){
/*
  @param array, array of the collumns of each process
  @param matrix, matrix to be separated
  @param order, order of the matrix
  @param np, number of processes created
  @param myrank, number of the process

  Sepatares matrix between the processes
*/
    /* case order can't be divided by number of processes */
    if(order%np!=0){
      /*
        Logic behind the process (variables have bad names, they should be changed):
          order%np = collumns that can't be divided equally
          nDiv = what should be added to get to the next number greater than
            order that can be divided by np
          nDiv+order = next number that can be divided by np
          n = number of collumns distributed between the processes
          nNormal = number of processes that have 'n' number of collumns
          nRest = the process that have the rest collumns

          Example:
            order: 9 ; np = 5

            order%np => 9%5 = 4
            nDiv =>     5-4 = 1
            n =>        10/5 = 2
            nNormal =>  9/2 = 4
            nRest =>    9%2 = 1

        This logic should be revised too, it works based on faith.

      */
      int count;
      int nDiv = np - (order%np);
      int n = (nDiv+order)/np;
      int nNormal = order/n;
      int nRest = order%n;
      /*
        There are 3 kinds of processes from the logic:
          Type 0) they have n collumns and there are nNormal of them (myrank<nNormal)
          Type 1) they have nRest collumns and there is only one (myrank=nNOrmal)
          Type 2) they have 0 collumns and do nothing, processes created in excess
      */
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
      /*
        MPI directive needs the number of elements of each array indexed by the
          number of process. It also needs the position of the first element
          for each process (displs=displacement).

      */
      int *displs = (int *)calloc(np,sizeof(int));
      int *elements = (int *)calloc(np,sizeof(int));
      int i;
      for (i=0; i<nNormal; ++i) {
        displs[i] = i*order*n;
        elements[i] = order*n;
      }
      displs[nNormal] = (nNormal)*order*n;
      elements[nNormal] = order*nRest;
      printf("elem %d order = %d\n", elements[0], order*n);

      MPI_Scatterv(matrix, elements, displs, MPI_DOUBLE, *array, count, MPI_DOUBLE,
                                                        0, MPI_COMM_WORLD);
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
  /*
    Takes all the collumns distributed between the processes and unifies them
    into a unique matrix in the process 0.
    The logic here works similar to the function 'separateMatrix'. 
  */

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
