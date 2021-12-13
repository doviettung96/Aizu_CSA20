#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define N 10

int main(int argc, char **argv)
{
    int size, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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

    // printf("Counts: ");
    // for (int i = 0; i < size; ++i)
    // {
    //     printf("%d ", counts[i]);
    // }
    // printf("Displs: ");
    // for (int i = 0; i < size; ++i)
    // {
    //     printf("%d ", displs[i]);
    // }

    // Init global data values
    int *globaldata = NULL;
    if (rank == 0)
    {
        globaldata = malloc(N * sizeof(int));
        for (int i = 0; i < N; i++)
            globaldata[i] = 2 * i + 1;

        printf("Processor %d has data: ", rank);
        for (int i = 0; i < N; i++)
            printf("%d ", globaldata[i]);
        printf("\n");
    }

    int *localdata = (int *)malloc(counts[rank] * sizeof(int));

    MPI_Scatterv(globaldata, counts, displs, MPI_INT, localdata, counts[rank], MPI_INT, 0, MPI_COMM_WORLD);
    printf("Processor %d has data: ", rank);
    for (int i = 0; i < counts[rank]; ++i)
    {
        printf("%d ", localdata[i]);
        localdata[i] *= 2;
    }
    printf("\n");

    MPI_Gatherv(localdata, counts[rank], MPI_INT, globaldata, counts, displs, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        printf("Processor %d has data: ", rank);
        for (int i = 0; i < N; i++)
            printf("%d ", globaldata[i]);
        printf("\n");
    }

    if (rank == 0)
        free(globaldata);

    MPI_Finalize();
    return 0;
}