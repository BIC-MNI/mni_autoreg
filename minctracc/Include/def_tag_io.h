#ifndef  DEF_TAG_IO
#define  DEF_TAG_IO

#include <stdio.h>

#include "../make_xfm_file/def_main.h"

int  input_tag_points(
    FILE      *file,
    int       *n_volumes,
    int       *n_tag_points,
    double    ***tags_volume1,
    double    ***tags_volume2,
    char      ***labels );

void  free_tag_points(
    int       n_volumes,
    int       n_tag_points,
    double    **tags_volume1,
    double    **tags_volume2,
    char      **labels );

int  output_transform(
    FILE      *file,
    char      comments[],
    double    transform[3][4] );

int  input_tag_comments( 
    FILE *file , 
    global_data_struct *main );

int  output_non_linear_transform(
    FILE      *file,
    char      comments[],
    int num_rows, int dim,
    double     **transform );

int  output_original_coords(FILE  *file,
			    double **points,
			    int num_rows, int dim);

int  input_tps_transform(
    FILE      *file,
    int       *n_tag_points,
    int       *dim,
    double    ***deformations,
    double    ***coords);

#ifndef DEBUG_PRINT
#   define DEBUG_PRINT(str) if (debug) (void) fprintf (stderr,  str  );
#   define DEBUG_PRINT1(str,a1) if (debug) (void) fprintf (stderr,  str ,a1 );
#   define DEBUG_PRINT2(str,a1,a2) if (debug) (void) fprintf (stderr,  str ,a1,a2 );
#   define DEBUG_PRINT3(str,a1,a2,a3) if (debug) (void) fprintf (stderr,  str ,a1,a2,a3 );
#   define DEBUG_PRINT4(str,a1,a2,a3,a4) if (debug) (void) fprintf (stderr,  str ,a1,a2,a3,a4 );
#   define DEBUG_PRINT5(str,a1,a2,a3,a4,a5) if (debug) (void) fprintf (stderr,  str ,a1,a2,a3,a4,a5 );
#   define DEBUG_PRINT6(str,a1,a2,a3,a4,a5,a6) if (debug) (void) fprintf (stderr,  str ,a1,a2,a3,a4,a5,a6 );
#   define DEBUG_PRINT7(str,a1,a2,a3,a4,a5,a6,a7) if (debug) (void) fprintf (stderr,  str ,a1,a2,a3,a4,a5,a6,a7 );

#endif

#endif




