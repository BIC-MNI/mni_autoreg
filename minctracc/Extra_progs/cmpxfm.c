/*
 * Compare transforms contained in two XFM files, element-by-element.
 * Exit with success if no element differs from the corresponding element
 * by more than a specified amount.
 */

#include <config.h>
#include <float.h>
#include <math.h>
#include <Proglib.h>
#include <ParseArgv.h>
#include <volume_io.h>

#include "local_macros.h"

static char* ProgName;
static double linear_tolerance;        /* ... for linear matrix */
static double translation_tolerance;   /* ... for translation vector */
static int show_max = 0;               /* show maximum deviation */
static int show_all = 0;               /* show all deviations */


static ArgvInfo argTable[] = {
    { "-linear_tolerance", ARGV_FLOAT, NULL, (char*)&linear_tolerance,
      "Absolute tolerance for transformation matrix elements"},
    { "-translation_tolerance", ARGV_FLOAT, NULL, (char*)&translation_tolerance,
      "Absolute tolerance for transformation matrix elements"},
    { "-show_max", ARGV_CONSTANT, (char*)1, (char*)&show_max,
      "Show the maximum deviation"},
    { "-show_all", ARGV_CONSTANT, (char*)1, (char*)&show_all,
      "Show all deviations"},
    { NULL, ARGV_END, NULL, NULL, NULL } 
};


/* Given a file name, check that it is a linear transform,
   return pointer to (linear) VIO_Transform structure, if so.
*/
VIO_Transform* input_linear_transform( char* filename )
{
    VIO_General_transform* gt;

    ALLOC( gt, 1 );

    if ( input_transform_file( filename, gt ) != OK ) {
        fprintf( stderr, "%s: cannot read file %s\n",
                 ProgName, filename );
        return NULL;
    }

    if ( get_transform_type(gt) != LINEAR ) {
        fprintf( stderr, "%s: transform in %s is not LINEAR\n",
                 ProgName, filename );
        return NULL;
    }

    return get_linear_transform_ptr(gt);
}


int do_compare( VIO_Transform* xfm1, VIO_Transform* xfm2,
                int i_min, int i_cnt, 
                int j_min, int j_cnt,
                double tolerance )
{
    int i, j;
    int ret = 0;

    for( i = i_min; i < i_min+i_cnt; ++i ) {
        for( j = j_min; j < j_min+j_cnt; ++j ) {
            double diff = fabs( Transform_elem( *xfm1, i, j ) 
                                - Transform_elem( *xfm2, i, j ) );
            if ( diff > tolerance ) {
                ret = 1;
                if ( show_all ) {
                    fprintf( stderr, "%s: difference at (%d,%d) is %f.\n",
                             ProgName, i, j, diff );
                } else {
                    fprintf( stderr, "%s: transforms differ at element (%d,%d)\n", 
                             ProgName, i, j );
                    return 1;
                }
            }
        }
    }
    return ret;
}


int main( int ac, char* av[] ) {
    VIO_Transform* xfm1;
    VIO_Transform* xfm2;

    /* Set defaults before parsing arguments */
    ProgName = av[0];

    /* On POSIX systems, the smallest value that can be added to 1.0 to
       give a distinct float number is contained in the symbol FLT_EPSILON.
       This is in header <float.h>, which we obtain from <config.h>.
    
       Default to about half the number of bits of precision 
    */
    linear_tolerance = translation_tolerance = sqrt( FLT_EPSILON );

    if ( ParseArgv( &ac, av, argTable, 0 ) || (ac != 3 )) {
        fprintf( stderr, "\nUsage: %s [options] xfm1 xfm2\n",
                 ProgName );
        fprintf( stderr, "       %s -help\n\n", ProgName );
        return 1;
    }

    xfm1 = input_linear_transform(av[1]);
    xfm2 = input_linear_transform(av[2]);

    if ( xfm1 == NULL || xfm2 == NULL )
        return 1;

    return do_compare( xfm1, xfm2, 0, 3, 0, 3, linear_tolerance )
        || do_compare( xfm1, xfm2, 0, 3, 3, 1, translation_tolerance );
}
