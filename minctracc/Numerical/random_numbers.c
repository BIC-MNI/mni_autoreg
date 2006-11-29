/* ----------------------------- MNI Header -----------------------------------
@NAME       : random_numbers.c
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: collection of routines to manipulate random numbers using the
              builtin function drand48.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@COPYRIGHT  :
              Copyright 1996 Louis Collins, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.

@CREATED    : Jan 8, 1996 LC
@MODIFIED   : $Log: random_numbers.c,v $
@MODIFIED   : Revision 1.4  2006-11-29 09:09:33  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 1.3  2005/07/20 20:45:49  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 1.2  2002/03/26 14:15:41  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 1.1  1997/11/03 19:59:49  louis
@MODIFIED   : Initial revision
@MODIFIED   :
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Numerical/random_numbers.c,v 1.4 2006-11-29 09:09:33 rotor Exp $";
#endif


#include <volume_io.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>

/*
   return a random number, from a gaussian distribution with unit
   standard deviation.
*/

VIO_Real gaussian_random_w_std(VIO_Real sigma)
{
  static int iset=0;
  static VIO_Real gset;
  VIO_Real fac,r,v1,v2;
  
  if  (iset == 0) {
    do {
      v1=2.0*drand48()-1.0;
      v2=2.0*drand48()-1.0;
      r=v1*v1+v2*v2;
    } while (r >= 1.0 || r == 0.0);
    fac=sqrt(-2.0*log(r)/r);
    gset=v1*fac;
    iset=1;
    return sigma * v2*fac;
  } else {
    iset=0;
    return sigma * gset;
  }
}


/*
   return a random number, uniformly distributed in the range [0, 1[
*/
VIO_Real uniform_random_0_1()
{
  VIO_Real t;
  t = drand48();
  return(t);
}


/*
   return a random number, uniformly distributed in the range [min, max[
*/
VIO_Real uniform_random_in_range(VIO_Real min, VIO_Real max)
{
  VIO_Real t,r;

  r = uniform_random_0_1();

  if (min>max) {
    t = max + r * (min-max);
  }
  else {
    t = min + r * (max-min);
  }
  return(t);
}


void init_random()
{
  union
    {
      long   l;
      char   c[4];
    } seedval;
   
  time_t t;
  char tmp;
                                /* initialize drand function */
  t = 2*time(NULL);
  seedval.l = t; 
  
  tmp = seedval.c[0]; seedval.c[0] = seedval.c[3]; seedval.c[3] = tmp; 
  tmp = seedval.c[1]; seedval.c[1] = seedval.c[2]; seedval.c[2] = tmp;
  
  srand48(seedval.l);
}
