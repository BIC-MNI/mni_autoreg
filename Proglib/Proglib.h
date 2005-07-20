/* ----------------------------- MNI Header -----------------------------------
@NAME       : Proglib.h
@DESCRIPTION: Header file for Proglib.a
@GLOBALS    : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */

/* a few macros for historical reasons */
#define VOXEL_DATA(vol) ((vol)->array.data)

#define MNI_AUTOREG_LONG_VERSION \
   "Package " PACKAGE_STRING \
   ", compiled by " MNI_AUTOREG_COMPILE_USER \
   "@" MNI_AUTOREG_COMPILE_HOST \
   " (" MNI_AUTOREG_COMPILE_SYSTEM ")" \
   " on " MNI_AUTOREG_COMPILE_DATETIME


void  print_error_and_line_num( char format[], char *name, int line, ... );
void  print_version_info( char *version_string);


/*
 *  * Generate a history string consisting of the output of time_stamp(),
 *   * (which ends with a newline) followed by "(MNI_AUTOREG_LONG_VERSION)\n".
 *    */
char* history_string( int ac, char* av[] );

