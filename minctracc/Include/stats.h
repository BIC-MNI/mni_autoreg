/* ----------------------------- MNI Header -----------------------------------
@NAME       : stats.h
@DESCRIPTION: definition of structures and procedures used to
              calculate and report simple statistics.
@COPYRIGHT  :
              Copyright 1995 Louis Collins, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.

@CREATED    : February 23, 1996
@MODIFIED   : $Log: stats.h,v $
@MODIFIED   : Revision 1.4  2002-12-13 21:09:45  lenezet
@MODIFIED   : *** empty log message ***
@MODIFIED   :
@MODIFIED   : Revision 1.3  2002/03/26 14:15:37  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 1.2  2000/03/15 08:42:40  stever
@MODIFIED   : Code cleanup: all functions prototyped (except ParseArgs.c), no useless declarations, etc
@MODIFIED   :
@MODIFIED   : Revision 1.1  1997/11/03 19:52:55  louis
@MODIFIED   : Initial revision
@MODIFIED   :
---------------------------------------------------------------------------- */

#ifndef  DEF_STATS
#define  DEF_STATS

#include  <volume_io/basic.h>

typedef  struct
{
  char   name[512];
  Real   mean;
  Real   standard_deviation;
  Real   variance;
  Real   rms;
  Real   sum;
  Real   sum_squared;
  int    count;
  Real   max_val;
  Real   min_val;
} stats_struct;

public void init_stats(stats_struct *stat,
		  char         title[]);

public void tally_stats(stats_struct *stat,
		   Real         val);

public void report_stats(stats_struct *stat);

public void stat_title(void);


public Real stat_get_mean(stats_struct *stat);

public Real stat_get_rms(stats_struct *stat);

public Real stat_get_standard_deviation(stats_struct *stat);

public Real stat_get_variance(stats_struct *stat);

public int stat_get_count(stats_struct *stat);


#endif

