void boundary_conditions ( int nx, int nt, double x[], double t, double h[], 
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

    Input/output, double H[NX], the height, with H(1) and H(NX) 
    adjusted for boundary conditions.

    Input/output, double UH[NX], the mass velocity, with UH(1) 
    and UH(NX) adjusted for boundary conditions.
*/
