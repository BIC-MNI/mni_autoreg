#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "time_stamp.h"

/* ----------------------------- MNI Header -----------------------------------
@NAME       : time_stamp
@INPUT      : argc - number of arguments
              argv - list of arguments
@OUTPUT     : 
@RETURNS    : pointer to string containing time stamp.
@DESCRIPTION: Function to produce a time stamp string for a program.
              Returns a string of the form "date > command". The command
              is simply the concatenation of argv elements.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : February 1, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
char *time_stamp(int argc, char *argv[])
{
   char *str, *the_time;
   int length, i, last;
   static char separator[]={">>>"};
   time_t timer;

   /* Get the time, overwriting newline */
   timer = time(NULL);
   the_time = ctime(&timer);

   /* Get the total length of the string and allocate space */
   length=strlen(the_time) + strlen(separator) + 2;
   for(i=0; i<argc; i++) {
      length += strlen(argv[i]) + 1;
   }
   str = malloc(length);

   /* Copy the time and separator */
   (void) strcpy(str, the_time);
   str[strlen(str)-1]='\0';
   (void) strcat(str, separator);

   /* Copy the program name and arguments */
   for (i=0; i<argc; i++) {
      last = strlen(str);
      str[last]=' ';
      str[last+1]='\0';
      (void) strcat(str, argv[i]);
   }

   /* Add a terminating newline */
   last = strlen(str);
   str[last]='\n';
   str[last+1]='\0';


   return str;
}
