#include <stdio.h>
#include <stdlib.h>
int main()
{
    double *h;
    double *h_new;
    FILE *h_file;
    int i;
    int j;
    int j_max = 10;
    int n = 10;
    double x_delta = 1.0;
    double cfl = 1.0;
    /*
    Set the values of H at the initial time.
    */
    h = (double *)malloc((n + 2) * sizeof(double));
    h_new = (double *)malloc((n + 2) * sizeof(double));
    h[0] = 0.0;
    for (i = 1; i <= n; i++)
    {
        h[i] = 50.0;
    }
    h[n + 1] = 0.0;

    /*
    Save the results to a file at timestep=0
    */

    h_file = fopen("h_data.txt", "w");

    for (i = 1; i <= n; i++)
    {
        fprintf(h_file, "  %f", h[i]);
    }
    fprintf(h_file, "\n");

    for (j = 1; j <= j_max; j++)
    {
        /*
    Update the temperature based on the four point stencil.
    */
        for (i = 1; i <= n; i++)
        {
            h_new[i] = h[i] + cfl * (h[i - 1] - 2.0 * h[i] + h[i + 1]);
        }

        /*
    H at the extreme left and right boundaries was incorrectly computed
    using the differential equation.  Replace that calculation by
    the boundary conditions.
    */
        h_new[1] = 90.0;
        h_new[n] = 70.0;
        /*
    Update temperature.
    */
        for (i = 1; i <= n; i++)
        {
            h[i] = h_new[i];
        }

        for (i = 1; i <= n; i++)
        {
            fprintf(h_file, "  %f", h[i]);
        }
        fprintf(h_file, "\n");
    }

    fclose(h_file);
    free(h);
    free(h_new);

    return 0;
}