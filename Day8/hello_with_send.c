/****************************************************************
 ****
 **** This program file is part of the book 
 **** `Parallel programming with MPI and OpenMP'
 **** by Victor Eijkhout, eijkhout@tacc.utexas.edu
 ****
 **** copyright Victor Eijkhout 2012-7
 ****
 **** MPI Exercise
 ****
 ****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int main(int argc,char **argv) {

/**** your code here ****/
  int number_of_processes;
  int my_rank;
  int mpi_error_code;

  mpi_error_code = MPI_Init(&argc, &argv);
  mpi_error_code = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  mpi_error_code = MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

  if (){
	  /will be run by server process
  }else {
	  /will be run by client process
  }

  printf("Hello, world! My rank is %d among %d\n", my_rank, number_of_processes);
  mpi_error_code = MPI_Finalize();
  
  return 0;
}
