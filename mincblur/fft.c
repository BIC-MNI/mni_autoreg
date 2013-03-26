/* ----------------------------- MNI Header -----------------------------------
@NAME       : fft.c
@INPUT      : complex data signal stored in 
              real1 imag1 real2 imag2... realN imagN format, w_ith

              signal[1] = real1     <-- NOTE zero offset!
              signal[2] = imag1
              ...
              signal[2*N-1] = realN
              signal[2*N] = imagN

              where N is the number of points.
              
@OUTPUT     : forward Fourier transform if direction==1,
              inverse Fourier transform if direction==-1;
              data signal is overwritten with result.
@RETURNS    : nothing
@DESCRIPTION: 
@METHOD     : uses Cooley-Tukey type fast Fourier transform of
              data, constained to have length (stored in numpoints) equal
              to a power of two.
@GLOBALS    : none
@CALLS      : nothing
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

@CREATED    : Wed Sep  6 09:15:07 MET DST 1995
              this routine was inspired by pascal code written by
              D. CLARK in Doctor Dobbs 1984 and by  domaine
              software made available by R. Hellman 2/21/86 from:
                  qiclab.scn.rain.com:/pub/math
              in file fft.c

              the routine was modified (individual routines removed,
              recoded in single procedure, sin and cos calculated by a
              recurrence relation ) to increase speed.

@MODIFIED   : $Log: fft.c,v $
@MODIFIED   : Revision 96.3  2006-11-28 09:12:21  rotor
@MODIFIED   :  * fixes to allow clean compile against minc 2.0
@MODIFIED   :
@MODIFIED   : Revision 96.2  2004/02/12 05:53:48  rotor
@MODIFIED   :  * removed /static defs
@MODIFIED   :
@MODIFIED   : Revision 96.1  2000/01/28 16:21:35  stever
@MODIFIED   : Revamped configure process
@MODIFIED   :
@MODIFIED   : Revision 96.0  1996/08/21 18:22:24  louis
@MODIFIED   : Release of MNI_AutoReg version 0.96
@MODIFIED   :
 * Revision 9.6  1996/08/21  18:22:17  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:16:28  louis
 * Never released version 0.95
 *
 * Revision 1.3  1996/08/12  14:16:20  louis
 * Pre-release
 *
 * Revision 1.2  1995/09/18  06:45:42  collins
 * this file is a working version of mincblur.  All references to numerical
 * recipes routines have been removed.  This version is included in the
 * package mni_reg-0.1i
 *
 * Revision 1.2  1995/09/18  06:45:42  collins
 * this file is a working version of mincblur.  All references to numerical
 * recipes routines have been removed.  This version is included in the
 * package mni_reg-0.1i
 *
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/mincblur/fft.c,v 96.3 2006-11-28 09:12:21 rotor Exp $";
#endif


#include <config.h>
#include <math.h>

#define PI2 6.28318530717959

void fft1(float *signal, int numpoints, int direction)
{
  int 
    n, m, mmax, 
    i, j, istep;
  double 
    sin_a,
    wp_r,wp_i,
    w_r,w_i,
    angle;
  float 
    temp,
    temp_real,
    temp_imag;
                                /* scramble the entries into
                                   bit reversed order */
  n = numpoints << 1;
  j = 1;
  for (i = 1; i < n; i += 2) {
    if (j > i) {                /* swap entries i and j */
      temp=signal[j];   signal[j]=signal[i];     signal[i]=temp;
      temp=signal[j+1]; signal[j+1]=signal[i+1]; signal[i+1]=temp;
    }

    m = n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;

  }

                                /* calculate the butterflies, but
                                   leave normalization to calling
                                   program if this is an inverse
                                   xform.  */
  mmax = 2;
  while (n > mmax) {
    angle = PI2/(direction*mmax); /* set up trig recurrance */
    sin_a = sin(0.5*angle);        
    wp_r = -2.0*sin_a*sin_a;        
    wp_i = sin(angle);

    istep = mmax<<1;                /* increment on i      */
    w_r = 1.0;                        /* initial sin and cos */
    w_i = 0.0;

    for (m = 1; m<mmax; m+=2) {        
      for (i = m; i<=n; i+=istep) {
        j = i+mmax;
        temp_real    = w_r*signal[j]   - w_i*signal[j+1]; /* mult */
        temp_imag    = w_r*signal[j+1] + w_i*signal[j]; 
        signal[j]    = signal[i]       - temp_real;       /* sub */
        signal[j+1]  = signal[i+1]     - temp_imag;
        signal[i]   += temp_real;                         /* add */
        signal[i+1] += temp_imag;
      }
      sin_a = w_r;
      w_r = sin_a*wp_r - w_i*wp_i   + w_r; /* trig recurrence  */
      w_i = w_i*wp_r   + sin_a*wp_i + w_i;
    }
    mmax = istep;
  }
}

