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
@MODIFIED   : Revision 1.10  1996-08-12 14:16:19  louis
@MODIFIED   : Pre-release
@MODIFIED   :
---------------------------------------------------------------------------- */

#ifndef public
   #define public
   #define private static
#endif

public void  muli_vects(float *r, float *s1, float *s2, int n);
public int   next_power_of_two(int x);
public float normal_dist(float c, float fwhm, float mu, float x);
public float rect_dist(float c, float fwhm, float mu, float x);
public void  make_kernel_FT(float *kern, int size, float vsize);
public void  make_kernel(float *kern, float vsize, float fwhm, int size, int type);


/*
  I redefine X, Y and Z since I use them differently than is defined in
  volume_io.  
*/
#ifdef X
# undef X
#endif
#define X  2
#ifdef y
# undef y
#endif
#define y  1
#ifdef z
# undef z
#endif
#define z  0

#include "kernel.h"


