/* ----------------------------- MNI Header -----------------------------------
@NAME       : mincblur.h
@DESCRIPTION: header file for mincblur.c
@CREATED    : Wed Jun 23 09:04:34 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */

public Status blur3D_volume(Volume data,
			    float kernel1, 
			    char *outfile, 
			    int ndim);

public Status gradient3D_volume(FILE *ifd, 
				Volume data, 
				char *outfile, 
				int ndim);

