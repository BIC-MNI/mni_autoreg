#include  <stdarg.h>

/* private  void  (*print_function) ( char [] ); */
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
    char     print_buffer[10000];


				/* do the standard error message: */
    sprintf (print_buffer,"Error in %s in file %s, line %d\n",prog_name,name,line);

/*
    if( print_function != NULL )
        (*print_function) ( print_buffer );
    else
*/
        (void) printf( "%s", print_buffer );



				/* now to the variable argument part! */
    va_start( ap, line );
    (void) vsprintf( print_buffer, format, ap );
    va_end( ap );

/*
    if( print_function != NULL )
        (*print_function) ( print_buffer );
    else
*/
        (void) printf( "%s\n", print_buffer );

    exit(-1);
}
