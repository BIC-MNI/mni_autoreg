/* ----------------------------- MNI Header -----------------------------------
@NAME       : blur_support.c
@DESCRIPTION: prototypes for the support routines used by the blurring and 
              gradient procedures.
@CREATED    : Wed Jun 23 09:04:34 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */

#ifndef public
   #define public
   #define private static
#endif

public void  muli_vects(float *r, float *s1, float *s2, int n);
public int   next_power_of_two(int x);
public float normal_dist(float c, float fwhm, float mu, float x);
public void  make_kernel_FT(float *kern, int size);
public void  make_kernel(float *kern, float vsize, float fwhm, int size);

