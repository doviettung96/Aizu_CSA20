#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]);

int main(int argc, char *argv[])
{
  int rank;
  int size;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int N = atoi(argv[1]);
  int T = atoi(argv[2]);

  // Load balancing
  int counts[size];
  int displs[size];
  int quotient = N / size;
  int remainder = N % size;
  for (int i = 0; i < size; ++i)
  {
    counts[i] = quotient;
    if (remainder > 0)
    {
      counts[i]++;
      remainder--;
    }

    displs[0] = 0;
    if (i > 0)
    {
      displs[i] = displs[i - 1] + counts[i - 1];
    }
  }

  double *global_h = NULL;
  FILE *h_file;
  if (rank == 0)
  {
    // init global data values
    global_h = malloc(N * sizeof(double));
    for (int i = 0; i < N; ++i)
    {
      global_h[i] = 50.0;
    }

    // write the first result to file
    h_file = fopen("h_data_mpi.txt", "w");
    for (int i = 0; i < N; ++i)
    {
      fprintf(h_file, " %f", global_h[i]);
    }
    fprintf(h_file, "\n");
  }

  // run the main algorithm for each process
  MPI_Status status;
  int tag;

  double *h;
  double *h_new;
  int n = counts[rank];
  double cfl = 1.0;

  /*
  Set the values of H at the initial time.
*/
  h = (double *)malloc((n + 2) * sizeof(double));
  h_new = (double *)malloc((n + 2) * sizeof(double));
  h[0] = 0.0;
  for (int i = 1; i <= n; i++)
  {
    h[i] = 50.0;
  }
  h[n + 1] = 0.0;

  /*
  Compute the values of H at the next time, based on current data.
*/
  for (int j = 1; j <= T; j++)
  {
    if (rank > 0)
    {
      tag = 1;
      MPI_Send(&h[1], 1, MPI_DOUBLE, rank - 1, tag, MPI_COMM_WORLD);
    }
    if (rank < size - 1)
    {
      tag = 1;
      MPI_Recv(&h[n + 1], 1, MPI_DOUBLE, rank + 1, tag, MPI_COMM_WORLD, &status);
    }
    if (rank < size - 1)
    {
      tag = 2;
      MPI_Send(&h[n], 1, MPI_DOUBLE, rank + 1, tag, MPI_COMM_WORLD);
    }
    if (rank > 0)
    {
      tag = 2;
      MPI_Recv(&h[0], 1, MPI_DOUBLE, rank - 1, tag, MPI_COMM_WORLD, &status);
    }

    /*
  Update the temperature based on the four point stencil.
*/
    for (int i = 1; i <= n; i++)
    {
      h_new[i] = h[i] + cfl * (h[i - 1] - 2.0 * h[i] + h[i + 1]);
    }

    /*
  H at the extreme left and right boundaries was incorrectly computed
  using the differential equation.  Replace that calculation by
  the boundary conditions.
*/
    if (rank == 0)
    {
      h_new[1] = 90.0;
    }
    if (rank == size - 1)
    {
      h_new[n] = 70.0;
    }
    /*
  Update temperature.
*/
    for (int i = 1; i <= n; i++)
    {
      h[i] = h_new[i];
    }
    /*
    Gather the results from child processes (from index 1 to n) to global_h to print to file
*/
    MPI_Gatherv(h + 1, n, MPI_DOUBLE, global_h, counts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
      for (int i = 0; i < N; i++)
      {
        fprintf(h_file, "  %f", global_h[i]);
      }
      fprintf(h_file, "\n");
    }
  }

  if (rank == 0)
  {
    fclose(h_file);
  }

  free(h);
  free(h_new);
  if (rank == 0)
    free(global_h);

  MPI_Finalize();

  return 0;
}
