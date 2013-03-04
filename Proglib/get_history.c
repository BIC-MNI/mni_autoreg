#include <config.h>
#include <string.h>
#include <stdlib.h>

#include <minc.h>
#include <time_stamp.h>

char* history_string( int ac, char* av[] )
{
    char* stamp = (char *)time_stamp( ac, av );
    int len = strlen(stamp) + strlen(MNI_AUTOREG_LONG_VERSION) + 4;
    char* s = malloc(len);

    if ( s == NULL )
	print_error_and_line_num( "cannot malloc(%d)", __FILE__, __LINE__, len );

    s[0] = 0;
    strcat(s, stamp );
    free( stamp );

    strcat(s, "(" );
    strcat(s, MNI_AUTOREG_LONG_VERSION);
    strcat(s, ")\n" );
    
    return s;
}
