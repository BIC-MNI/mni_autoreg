/*
------------------------------ MNI Header ----------------------------------
@NAME       : extras.h
@DESCRIPTION: prototypes for Optimize/extras.c              
@CREATED    : Mon Nov  3, 1997 , Louis Collins
@MODIFIED   : not yet!
@VERSION    : $Id: extras.h,v 1.1 1997-11-03 19:52:55 louis Exp $
#-----------------------------------------------------------------------------
*/


public void report_time(long start_time, STRING text);

public void init_the_volume_to_zero(Volume volume);

public Real get_volume_maximum_real_value(Volume volume);

public void save_data(char *basename, int i, int j,
		      General_transform *transform);
