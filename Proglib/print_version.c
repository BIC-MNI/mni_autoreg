#include  <stdio.h>

extern   char *prog_name;

/* ----------------------------- MNI Header -----------------------------------
@NAME       : print_version_info
@INPUT      : version info string
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: prints out program version information and exits
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1996 louis colling
@MODIFIED   : 

---------------------------------------------------------------------------- */

void  print_version_info( char *version_string )
{

    (void) printf( "The program <%s> was built from:\n", prog_name);
   
    (void) printf( "%s\n", version_string );

    exit(-1);
}
