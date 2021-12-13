#include <stdio.h>
#include <stdlib.h>

#include "r8vec_write.h"
/******************************************************************************/

void r8vec_write ( char *output_filename, int n, double x[] )
{
  int j;
  FILE *output;
/*
  Open the file.
*/
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output file.\n" );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    fprintf ( output, "  %24.16g\n", x[j] );
  }
/*
  Close the file.
*/
  fclose ( output );

  return;
}
