/* ----------------------------- MNI Header -----------------------------------
@NAME       : stats.c
@DESCRIPTION: collection of routines to simplify the calculation of 
              mean, standard deviation, rms, etc of variables.
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
@MODIFIED   : $Log: stats.c,v $
@MODIFIED   : Revision 1.1  1997-11-03 19:59:49  louis
@MODIFIED   : Initial revision
@MODIFIED   :
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Numerical/stats.c,v 1.1 1997-11-03 19:59:49 louis Exp $";
#endif

#include <internal_volume_io.h>
#include <config.h>
#include <stats.h>

public void init_stats(stats_struct *stat,
		  char         title[])
{

  ALLOC ( stat->name, strlen( title )+1 );
  (void) strcpy( stat->name, title);
  stat->mean	           = 0.0;
  stat->standard_deviation = 0.0;
  stat->rms	           = 0.0;
  stat->sum                = 0.0;
  stat->sum_squared	   = 0.0;
  stat->count	           = 0;
  stat->min_val            = DBL_MAX;
  stat->max_val            = -DBL_MAX;
}

public void tally_stats(stats_struct *stat,
		   Real         val)
{
  stat->count++;
  stat->sum         += val;
  stat->sum_squared += val*val;
  if (val>stat->max_val) stat->max_val = val;
  if (val<stat->min_val) stat->min_val = val;
}

private void calc_stats(stats_struct *stat)
{
  if (stat->count>0) {
    stat->mean = stat->sum / (Real)stat->count;
    stat->rms  = sqrt( stat->sum_squared / (Real)stat->count);
    stat->variance  = (stat->sum_squared * stat->count - stat->sum*stat->sum) / 
                 ((Real)stat->count * (Real)(stat->count-1.0));
    if (stat->variance >= 0.0)
      stat->standard_deviation = sqrt(stat->variance);
    else
      stat->standard_deviation = 0.0;
	
  }  
  else {
    print ("warning: calc_stats(%s) called with zero counter\n", stat->name);
  }
}

public void report_stats(stats_struct *stat)
{
  if (stat != (stats_struct *)NULL) {
    if (stat->count>0) {
      calc_stats(stat);
      print ("%14s %12f %12f %12f %12f %12f %12d\n", 
	     stat->name, 
	     stat->mean, 
	     stat->standard_deviation, 
	     stat->rms,
	     stat->min_val,
	     stat->max_val,
	     stat->count);
    }
    else {
      print ("warning: report_stats(%s) called with zero counter\n", stat->name);
    }
  }
  else {
    print ("warning: report_stats() called with NULL stat_struct\n");
  }
}

public void stat_title()
{
  print ("Statistics report:\n");
  print ("%14s %12s %12s %12s %12s %12s %12s\n","variable name","mean","std","rms","min","max","cnt");
}


public Real stat_get_mean(stats_struct *stat)
{
  calc_stats(stat);
  return(stat->mean);
}

public Real stat_get_rms(stats_struct *stat)
{
  calc_stats(stat);
  return(stat->rms);
}

public Real stat_get_standard_deviation(stats_struct *stat)
{
  calc_stats(stat);
  return(stat->standard_deviation);
}

public Real stat_get_variance(stats_struct *stat)
{
  calc_stats(stat);
  return(stat->variance);
}

public int stat_get_count(stats_struct *stat)
{
  return(stat->count);
}

