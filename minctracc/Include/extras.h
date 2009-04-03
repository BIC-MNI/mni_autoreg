/*
------------------------------ MNI Header ----------------------------------
@NAME       : extras.h
@DESCRIPTION: prototypes for Optimize/extras.c              
@CREATED    : Mon Nov  3, 1997 , Louis Collins
@MODIFIED   : not yet!
@VERSION    : $Id: extras.h,v 1.5 2009-04-03 18:36:59 louis Exp $
#-----------------------------------------------------------------------------
*/


void report_time(long start_time, VIO_STR text);

void init_the_volume_to_zero(VIO_Volume volume);

VIO_Real get_volume_maximum_real_value(VIO_Volume volume);

void save_data(char *basename, int i, int j,
                      VIO_General_transform *transform);

void get_volume_minimum_maximum_real_value(VIO_Volume volume, VIO_Real *min, VIO_Real *max);
