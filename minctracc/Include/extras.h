/*
------------------------------ MNI Header ----------------------------------
@NAME       : extras.h
@DESCRIPTION: prototypes for Optimize/extras.c              
@CREATED    : Mon Nov  3, 1997 , Louis Collins
@MODIFIED   : not yet!
@VERSION    : $Id: extras.h,v 1.2 2004-02-12 05:54:16 rotor Exp $
#-----------------------------------------------------------------------------
*/


void report_time(long start_time, STRING text);

void init_the_volume_to_zero(Volume volume);

Real get_volume_maximum_real_value(Volume volume);

void save_data(char *basename, int i, int j,
		      General_transform *transform);
