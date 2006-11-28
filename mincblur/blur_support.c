/* ----------------------------- MNI Header -----------------------------------
@NAME       : blur_support.c
@DESCRIPTION: prototypes for the support routines used by the blurring and 
              gradient procedures.
@CREATED    : Wed Jun 23 09:04:34 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */
#include <math.h>
#include <stdio.h>
#include "kernel.h"
#include <string.h>

#define PI 3.1415927

/* ----------------------------- MNI Header -----------------------------------
@NAME       : muli_vects
@INPUT      : s1 - a zero offset array containing real,imag,real,imag
              s2 - a zero offset array containing real,imag,real,imag
              n  - the number of complex pairs to be multiplied
@OUTPUT     : r  - the result of the multiplication, a zero offset array 
                   containing real,imag,real,imag
@RETURNS    : nothing
@DESCRIPTION: 
              for c = a*b, all real:

                      c.r=a.r*b.r-a.i*b.i;a
                      c.i=a.i*b.r+a.r*b.i;
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Wed Jun 23 09:04:34 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */
void muli_vects(float *r, float *s1, float *s2, int n)
{
  int i;
  float r1,i1,r2,i2;

  r++; s1++; s2++; /* all arrays start at r[1],s1[1] and s2[1], where the real
                      part is r[1] and the imag in r[2], and so on... */

  if (r!=s1 && r!=s2) {                /* if separate arrays */
    for (i=0; i< n; ++i) { 
      *r = (*(s1) * *(s2))   - (*(s1+1) * *(s2+1)); 
      r++;
      *r = (*(s1+1) * *(s2)) + (*(s1) * *(s2+1)); 
      r++;
      s1++;s1++;
      s2++;s2++;
    } 
  }
  else {                        /* if one array=result array */
    for (i=0; i< n; ++i) { 
      r1= *(s1); i1= *(s1+1);
      r2= *(s2); i2= *(s2+1);
      *r = (r1 * r2)   - (i1 * i2); 
      r++;
      *r = (i1 * r2) + (r1 * i2);
      r++;
      s1++;s1++;
      s2++;s2++;
    } 
  }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : next_power_of_two
@INPUT      : x - integer number (positive)
@OUTPUT     : 
@RETURNS    : int - the next higher power of two
@DESCRIPTION: the routine returns the smallest number n > x, such that n = 2^?
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Wed Jun 23 09:04:34 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */
/************************************************************/
/* find the next highest power of 2, greater than or equal  */
/* to x, and return 2^n, when n is the required power.      */
/************************************************************/
int next_power_of_two(int x)
{
  int 
    n, power_of_two;
  
  power_of_two = 1;
  n = 0;
  
  while (power_of_two<x && n<32) {
    power_of_two *= 2;
    n++;
  }
  
  return(power_of_two);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : normal_dist
@INPUT      : c    - height of gaussian
              fwhm - full wifth half max of gaussian
              mu   - center of gaussian
              x    - value of x
@OUTPUT     : 
@RETURNS    : value of gaussian evaluated at x
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Wed Jun 23 09:04:34 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */
/************************************************************/
/* return the value of the normal dist at x, given c,sigma  */
/* and mu ----   all in mm                                  */
/************************************************************/
float normal_dist(float c, float fwhm, float mu, float x)
{
  float sigma,t1,t2,t3,f;
  
  sigma = fwhm/2.35482;
  
  if (sigma==0) {
    if (x==mu)
      f = c / (sqrt(2*PI) * sigma);      /* !!!!!!!!!!! */
    else
      f = 0;
  }
  else {
    t1 = c / (sqrt(2*PI) * sigma);
    t2 = (x-mu)*(x-mu)/(2*sigma*sigma);
    t3 = exp(-t2);
    f = t1*t3;
  }
  
  return(f);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : rect_dist
@INPUT      : c    - height of rect function
              fwhm - width of rect function
              mu   - center of rect function
              x    - value of x
@OUTPUT     : 
@RETURNS    : value of rect function evaluated at x
@CALLS      : 
@CREATED    : Wed Jun 23 09:04:34 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */
/************************************************************/
/* return the value of the normal dist at x, given c,sigma  */
/* and mu ----   all in mm                                  */
/************************************************************/
float rect_dist(float c, float fwhm, float mu, float x)
{
  
  float t;
  
  t = x-mu;
  
/*  t = t<0 ? -t : t ; */
  
  if ( t >= -fwhm/2  && t < fwhm/2 ) {
    return(  (float) (c/fwhm) );
  }
  else
    return ( (float) 0.0 );
  
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_kernel_FT
@INPUT      : kern - a zero offset array containing real,imag,real,imag
                     in which will be stored the kernel for the dirivitive
              size - the number of complex numbers in the kernel array
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Wed Jun 23 09:04:34 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */
/************************************************************/
void make_kernel_FT(float *kern, int size, float vsize)
{
  
  int 
    kindex,k;
  float
    factor,
    f_sample_size;
  
  (void)memset(kern,(int)0,(size_t)((2*size+1)*sizeof(float)));
  
  f_sample_size = 1.0/(vsize*size);
  factor = -2.0 * PI * f_sample_size;
  
  
  for ( k = -size/2; k<size/2; ++k) {
    kindex = ((k + size) % size)*2 +1; 
    kern[kindex+1] = factor*k;
  }
  
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_kernel
@INPUT      : kern - a zero offset array containing real,imag,real,imag
                     in which will be stored the Gaussian kernel for convolution
              vsize- the size (in mm) of the sample along the kern array
              fwhm - full-width-half-maximum of gaussian (in mm)
              size - the number of complex numbers in the kernel array
@OUTPUT     : kern - the Gaussian kernel used for convolution
@RETURNS    : nothing
@DESCRIPTION: 

   note that kern (the convolution kernel) goes into the array so that the
   peak is at kern[1] (remember unit offset), with the
   positive half of the kernel running from kern[1] to kern[n/2] and the
   negative half running from kern[array_size_pow2 - n/2] to kern[array_size_pow2]:

  -                            size / 2
   \                              V                                 /
  --|-----------------------------+--------------------------------|-
     \_/                                                        \_/ 
  ^                                                                 ^
  0                                                            size ^  

            this is not right? ---^
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Wed Jun 23 09:04:34 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */
void make_kernel(float *kern, float vsize, float fwhm, int size, int type)
{

  int kindex,k;
  float c,r,max;

  (void)memset(kern,(int)0,(size_t)((2*size+1)*sizeof(float)));
  
  for ( k = -size/2; k<size/2; ++k) {

    kindex = ((k + size) % size)*2 +1;


    switch (type) {
    case KERN_GAUSSIAN: kern[kindex] = normal_dist(1.0*vsize,fwhm,0.0,(float)(vsize*k)); break;
    case KERN_RECT:     kern[kindex] = rect_dist(1.0*vsize,fwhm,0.0,(float)(vsize*k)); break;
    default: 
      {
        (void) fprintf (stderr,"Illegal kernel type = %d\n",type);
        (void) fprintf (stderr,"Impossible error in make_kernel(), line %d of %s\n",
                        __LINE__,__FILE__);
        k = size/2;
      }
    }
  }

  if (type == KERN_RECT) {
    max = -1e37;
    for(k=1; k<=size; k++)
      if (kern[k]>max) max = kern[k];
    if (max != 0.0) {
      for(k=1; k<=size; k++)
        kern[k] = kern[k] / max;
      c = (int)( fwhm/vsize )==(int)(fwhm/vsize+0.5) ? (int)( fwhm/vsize )+1 : (int)( fwhm/vsize )+2;
      for(k=1; k<=size; k++)
        kern[k] = kern[k] / (fwhm/vsize + 1);
      r = 0.0;
      for(k=1; k<=size; k++)
        r +=kern[k];
      kern[ (int)(c/2 + 0.5) + 1] += (1.0 - r)/2;
      kern[ size - (int)(c/2 + 0.5) + 1] += (1.0 - r)/2;
      
    }
  }

}


