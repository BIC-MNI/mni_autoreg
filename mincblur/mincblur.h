/* ----------------------------- MNI Header -----------------------------------
@NAME       : mincblur.h
@DESCRIPTION: header file for mincblur.c
@CREATED    : Wed Jun 23 09:04:34 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */

public Status blur3D_volume(Volume data,
			    float kernel1, 
			    char *outfile, 
			    int ndim, int kernel_type, char *history);

public Status gradient3D_volume(FILE *ifd, 
				Volume data, 
				char *outfile, 
				int ndim,
				char *history);

public void calc_gradient_magnitude(char *infilename, char *history);

#include<print_error.h>

#define INTERNAL_X  2
#define INTERNAL_Y  1
#define INTERNAL_Z  0
