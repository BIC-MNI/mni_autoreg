/* ----------------------------- MNI Header -----------------------------------
@NAME       : blur_support.c
@DESCRIPTION: prototypes for the support routines used by the blurring and 
              gradient procedures.
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

@CREATED    : Wed Jun 23 09:04:34 EST 1993 Louis Collins

@MODIFIED   : $Log: blur_support.h,v $
@MODIFIED   : Revision 96.2  2004-02-12 05:53:48  rotor
@MODIFIED   :  * removed public/private defs
@MODIFIED   :
@MODIFIED   : Revision 96.1  2000/01/27 18:03:51  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :
 * Revision 96.0  1996/08/21  18:22:24  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:22:16  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:16:28  louis
 * Never released version 0.95
 *
 * Revision 1.10  1996/08/12  14:16:19  louis
 * Pre-release
 *
---------------------------------------------------------------------------- */

void  muli_vects(float *r, float *s1, float *s2, int n);
int   next_power_of_two(int x);
float normal_dist(float c, float fwhm, float mu, float x);
float rect_dist(float c, float fwhm, float mu, float x);
void  make_kernel_FT(float *kern, int size, float vsize);
void  make_kernel(float *kern, float vsize, float fwhm, int size, int type);


/*
  I redefine X, Y and Z since I use them differently than is defined in
  volume_io.  
*/
#ifdef X
# undef X
#endif
#define X  2

#ifdef Y
# undef Y
#endif
#define Y  1

#ifdef Z
# undef Z
#endif
#define Z  0

#include "kernel.h"


