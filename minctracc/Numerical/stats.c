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
@MODIFIED   : Revision 1.8  2006-11-30 09:07:32  rotor
@MODIFIED   :  * many more changes for clean minc 2.0 build
@MODIFIED   :
@MODIFIED   : Revision 1.7  2006/11/29 09:09:33  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 1.6  2005/07/20 20:45:49  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 1.5  2004/02/13 00:17:15  rotor
@MODIFIED   :  * removed /static defs
@MODIFIED   :
@MODIFIED   : Revision 1.4  2002/12/13 21:16:31  lenezet
@MODIFIED   : nonlinear in 2D has changed. The option -2D-non-lin is no more necessary. The grid transform has been adapted to feet on the target volume whatever is size. The Optimization is done on the dimensions for which "count" is greater than 1.
@MODIFIED   :
@MODIFIED   : Revision 1.3  2002/03/26 14:15:42  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 1.2  2000/03/15 08:42:44  stever
@MODIFIED   : Code cleanup: all functions prototyped (except ParseArgs.c), no useless declarations, etc
@MODIFIED   :
@MODIFIED   : Revision 1.1  1997/11/03 19:59:49  louis
@MODIFIED   : Initial revision
@MODIFIED   :
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/minctracc/Numerical/stats.c,v 1.8 2006-11-30 09:07:32 rotor Exp $";
#endif

#include <config.h>
#include <float.h>
#include <volume_io.h>
#include <stats.h>

#include "local_macros.h"

void init_stats(stats_struct *stat,
                  char         title[])
{

  /*  ALLOC ( stat->name, strlen( title )+1 ); */

  (void) strcpy( stat->name, title);

  stat->mean                   = 0.0;
  stat->standard_deviation = 0.0;
  stat->rms                   = 0.0;
  stat->sum                = 0.0;
  stat->sum_squared           = 0.0;
  stat->count                   = 0;
  stat->min_val            = DBL_MAX;
  stat->max_val            = -DBL_MAX;
}

void deinit_stats( stats_struct *stat )
{
  FREE (stat->name);
}


void tally_stats(stats_struct *stat,
                   VIO_Real         val)
{
  stat->count++;
  stat->sum         += val;
  stat->sum_squared += val*val;
  if (val>stat->max_val) stat->max_val = val;
  if (val<stat->min_val) stat->min_val = val;
}

static void calc_stats(stats_struct *stat)
{
  if (stat->count>0) {
    stat->mean = stat->sum / (VIO_Real)stat->count;
    stat->rms  = sqrt( stat->sum_squared / (VIO_Real)stat->count);
    stat->variance  = (stat->sum_squared * stat->count - stat->sum*stat->sum) / 
                 ((VIO_Real)stat->count * (VIO_Real)(stat->count-1.0));
    if (stat->variance >= 0.0)
      stat->standard_deviation = sqrt(stat->variance);
    else
      stat->standard_deviation = 0.0;
        
  }  
  else {
    print ("warning: calc_stats(%s) called with zero counter\n", stat->name);
  }
}

void report_stats(stats_struct *stat)
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

void stat_title(void)
{
  print ("Statistics report:\n");
  print ("%14s %12s %12s %12s %12s %12s %12s\n","variable name","mean","std","rms","min","max","cnt");
}


VIO_Real stat_get_mean(stats_struct *stat)
{
  calc_stats(stat);
  return(stat->mean);
}

VIO_Real stat_get_rms(stats_struct *stat)
{
  calc_stats(stat);
  return(stat->rms);
}

VIO_Real stat_get_standard_deviation(stats_struct *stat)
{
  calc_stats(stat);
  return(stat->standard_deviation);
}

VIO_Real stat_get_variance(stats_struct *stat)
{
  calc_stats(stat);
  return(stat->variance);
}

int stat_get_count(stats_struct *stat)
{
  return(stat->count);
}

