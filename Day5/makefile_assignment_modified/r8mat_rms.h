#ifndef H_R8MAT_RMS_H
#define H_R8MAT_RMS_H
#define NX 11
#define NY 11
double r8mat_rms(int nx, int ny, double a[NX][NY]);
#endif
/******************************************************************************/
/*
  Purpose:

    R8MAT_RMS returns the RMS norm of a vector stored as a matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 March 2003

  Author:

    John Burkardt

  Parameters:

    Input, int NX, NY, the number of rows and columns in A.

    Input, double A[NX][NY], the vector.

    Output, double R8MAT_RMS, the root mean square of the entries of A.
*/
/******************************************************************************/