#include <stdarg.h>
#include <stdio.h>

/* void  (*print_function) ( char [] ); */
extern   char *prog_name;

/* ----------------------------- MNI Header -----------------------------------
@NAME       : print_error_and_line_num
@INPUT      : exactly same arguments as printf
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: prints the arguments to a temporary string buffer, then either
              printf's the or calls the user installed function to output
              the string.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : Mon Sep 13 13:56:49 EST 1993 Louis Collins
   modified dave print() function for print error.
Wed Apr 26 14:55:51 MET DST 1995
   changed name from print_error to print_error_and_line_num

---------------------------------------------------------------------------- */

/* VARARGS */
void  print_error_and_line_num( char format[], char *name, int line, ... )
{
    va_list  ap;

    fprintf( stderr, "Error in %s in file %s, line %d\n", prog_name, name, line );

    va_start( ap, line );
    vfprintf( stderr, format, ap );
    va_end( ap );

    exit(-1);
}
