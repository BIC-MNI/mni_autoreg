/*
------------------------------ MNI Header ----------------------------------
@NAME       : extras.h
@DESCRIPTION: prototypes for Optimize/extras.c              
@CREATED    : Mon Nov  3, 1997 , Louis Collins
@MODIFIED   : not yet!
@VERSION    : $Id: extras.h,v 1.4 2006-11-30 09:07:31 rotor Exp $
#-----------------------------------------------------------------------------
*/


void report_time(long start_time, VIO_STR text);

void init_the_volume_to_zero(VIO_Volume volume);

VIO_Real get_volume_maximum_real_value(VIO_Volume volume);

void save_data(char *basename, int i, int j,
                      VIO_General_transform *transform);
