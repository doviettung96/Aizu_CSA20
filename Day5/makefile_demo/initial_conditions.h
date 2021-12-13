void initial_conditions ( int nx, int nt, double x[], double t, double h[], 
  double uh[] );
/******************************************************************************/
/*
  Purpose:

    INITIAL_CONDITIONS sets the initial conditions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 July 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NX, the number of spatial nodes.

    Input, int NT, the number of times steps.

    Input, double X[NX], the coordinates of the nodes.

    Input, double T, the current time.

    Output, double H[NX], the initial height for all space.

    Output, double UH[NX], the initial mass velocity for all space.
*/
