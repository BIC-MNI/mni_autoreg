/* ----------------------------- MNI Header -----------------------------------
@NAME       : rotmat_to_ang.c
@INPUT      : rot      - rotation matrix (3x3 in zero offset form) calculated
                         by the calling program.
@OUTPUT     : ang      - vector giving rx,ry and rz rotation angles (in 
                         radians). This vector must be defined by the 
                         calling routine.
@RETURNS    : TRUE if ok, FALSE if error.
@DESCRIPTION: this routine extracts the rotation angles from the rotation
              matrix.  The rotation matrix is assumed to be a 3x3 matrix in
              zero offset form [1..3][1..3].  It is locally copied into a 
              4x4 homogeneous matrix for manipulation.

              we assume that the matrix rotation center is (0,0,0).
              we assume that the application of this matrix to a vector
                 is done with rot_mat*vec = premultiplication by matrix

              the resulting angles rx=ang[1],ry=ang[2],rz=ang[3], follow
              the following assumptions:

              -all rotations are assumed to be in the range -pi/2..pi/2
               routine returns FALSE is this is found not to hold
              -rotations are applied 1 - rx, 2 - ry and 3 - rz
              -applying these rotations to an identity matrix will
               result in a matrix equal to `rot' (the input matrix)
              -positive rotations are counter-clockwise when looking
               down the axis, from the positive end towards the origin.
              -I assume a coordinate system:
                          ^ Y
                          |
                          |
                          |
                          |
                          +---------> X
                         /
                        /
                       / Z  (towards the viewer).

              -The problem is posed as:  
                 given a rotation matrix ROT, determine the rotations
                 rx,ry,rz applied in order to give ROT
               solution:
                 assume the rot matrix is equivalent to a normalized
                 orthogonal local coord sys.  i.e.  row 1 of ROT is
                 the local x direction, row 2 is the local y and row 3
                 is the local z.
             
                 (note local is lower case, world is UPPER)

                 1- find RZ that brings local x into the XZ plane
                 2- find RY that brings local x*RZ onto X
                 3- find RX that brings local y*RZ*RY onto Y

                 the required rotations are -RX,-RY and -RZ!

@GLOBALS    : 
@CALLS      : mfmult(), rotx(),roty(),rotz()
@COPYRIGHT  :
              Copyright 1993 Louis Collins, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.

@CREATED    : Feb 9, 1992 lc
@MODIFIED   :  $Log: rotmat_to_ang.c,v $
@MODIFIED   :  Revision 96.6  2006-11-30 09:07:32  rotor
@MODIFIED   :   * many more changes for clean minc 2.0 build
@MODIFIED   :
@MODIFIED   :  Revision 96.5  2006/11/29 09:09:33  rotor
@MODIFIED   :   * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   :  Revision 96.4  2005/07/20 20:45:49  rotor
@MODIFIED   :      * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :      * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :      * Removed old VOLUME_IO cruft #defines
@MODIFIED   :      * Fixed up all Makefile.am's in subdirs
@MODIFIED   :      * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :      * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   :  Revision 96.3  2004/02/13 00:17:15  rotor
@MODIFIED   :   * removed /static defs
@MODIFIED   :
@MODIFIED   :  Revision 96.2  2002/03/26 14:15:41  stever
@MODIFIED   :  Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   :  Revision 96.1  1997/11/03 19:59:49  louis
@MODIFIED   :  - now include volume_io/internal_volume_io.h instead of volume_io.h
@MODIFIED   :
 * Revision 96.0  1996/08/21  18:21:58  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:21:55  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:15:50  louis
 * Never released version 0.95
 *
 * Revision 1.11  1996/08/12  14:15:47  louis
 * Pre-release
 *
 * Revision 1.10  1995/09/11  12:37:16  collins
 * All refs to numerical recipes routines have been replaced.
 * this is an updated working version - corresponds to mni_reg-0.1g
 *
 * Revision 1.9  1995/02/22  08:56:06  collins
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.8  94/04/06  11:48:49  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.7  94/02/21  16:36:53  louis
 * version before feb 22 changes
 * 
 * Revision 1.6  93/11/15  16:27:10  louis
 * working version, with new library, with RCS revision stuff,
 * before deformations included
 * 

Tue Jun  8 08:44:59 EST 1993 LC
   changes all vec*matrix to matrix*vec.  Std is premultiplication by matrix!

---------------------------------------------------------------------------- */


#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/minctracc/Numerical/rotmat_to_ang.c,v 96.6 2006-11-30 09:07:32 rotor Exp $";
#endif

#include <config.h>
#include <volume_io.h>
#include "local_macros.h"
#include "matrix_basics.h"

extern char *prog_name;


#define EPS  0.00000000001        /* epsilon, should be calculated */



VIO_BOOL rotmat_to_ang(float **rot, float *ang)
{

   float 
      rx,ry,rz,
      **t,**s,
      **R,
      **Rx,
      **Ry,
      **Rz,
      len,
      i,j,k;

   int
      m,n;

   VIO_ALLOC2D(t  ,5,5);        /* make two column vectors */
   VIO_ALLOC2D(s  ,5,5);

   VIO_ALLOC2D(R  ,5,5);        /* working space matrices */
   VIO_ALLOC2D(Rx ,5,5);
   VIO_ALLOC2D(Ry ,5,5);
   VIO_ALLOC2D(Rz ,5,5);

   nr_identf(R,1,4,1,4);        /* init R homogeneous matrix */

   for (m=1; m<=3; ++m)                /* copy rot matrix into R */
      for (n=1; n<=3; ++n)
         R[m][n] = rot[m][n];
   
/* ---------------------------------------------------------------
   step one,  find the RZ rotation reqd to bring 
              the local x into the world XZ plane
*/

   for (m=1; m<=3; ++m)                /* extract local x vector, ie the first column */
      t[m][1] = R[m][1];
   t[4][1] = 1.0;

   i = t[1][1];                        /* make local vector componants */
   j = t[2][1]; 
   k = t[3][1];

   if (i<EPS) {                        /* if i is not already in the positive X range, */
      print("WARNING: (%s:%d) %s\n",__FILE__, __LINE__,"step one: rz not in the range -pi/2..pi/2");
      return(FALSE);
   }

   len = sqrt(i*i + j*j);        /* length of vect x on XY plane */
   if (fabs(len)<EPS) {
      print("WARNING: (%s:%d) %s\n",__FILE__, __LINE__,"step one: length of vect x null.");
      return(FALSE);
   }

   if (fabs(i)>fabs(j)) {
      rz = fabs(asin((double)(j/len)));
   }
   else {
      rz = fabs(acos((double)(i/len)));
   }

   if (j>0)                        /* what is the counter clockwise angle */
      rz = -rz;                 /* necessary to bring vect x ont XY plane? */
      
  
/*---------------------------------------------------------------
   step two:  find the RY rotation reqd to align 
              the local x on the world X axis 

  (since i was positive above, RY should already by in range -pi/2..pi/2 
  but we'll check it  anyway)                                             */

   for (m=1; m<=3; ++m)                /* extract local x vector */
      t[m][1] = R[m][1];
   t[4][1] = 1.0;

   nr_rotzf(Rz,rz);             /* create the rotate Z matrix */
 
   nr_multf(Rz,1,4,1,4,  t,1,4,1,1,   s);   /* apply RZ, to get x in XZ plane */

   i = s[1][1];                        /* make local vector componants */
   j = s[2][1]; 
   k = s[3][1];

   if (i<EPS) {
      print("WARNING: (%s:%d) %s\n",__FILE__, __LINE__,"step two: ry not in the range -pi/2..pi/2");
      return(FALSE);
   }

   len = sqrt(i*i + k*k);                /* length of vect x in XZ plane, after RZ */

   if (fabs(len)<EPS) {
      print("WARNING: (%s:%d) %s\n",__FILE__, __LINE__,"step two: length of vect z null.");
      return(FALSE);
   }

   if (fabs(i)>fabs(k)) {
      ry = fabs(asin((double)(k/len)));
   }
   else {
      ry = fabs(acos((double)(i/len)));
   }

   /*    what is the counter clockwise angle necessary to bring  */
   /*    vect x onto X? */
   if (k < 0) { 
      ry = -ry;
   }

   /*--------------------------------------------------------------- */
   /*   step three,rotate around RX to */
   /*              align the local y with Y and z with Z */

   for (m=1; m<=3; ++m)                /* extract local z vector */
      t[m][1] = R[m][3];
   t[4][1] = 1.0;

   nr_rotyf(Ry,ry);             /* create the rotate Y matrix */

                                /* t =  roty(ry*180/pi) *(rotz(rz*180/pi) *r(3,:)); */
   nr_multf(Rz,1,4,1,4,  t,1,4,1,1,  s); /* apply RZ, to get x in XZ plane */
   nr_multf(Ry,1,4,1,4,  s,1,4,1,1,  t); /* apply RY, to get x onto X      */

   i = t[1][1];                        /* make local vector componants */
   j = t[2][1]; 
   k = t[3][1];

   len = sqrt(j*j + k*k);        /* length of vect x in Y,Z plane */

   if (fabs(len)<EPS) {
      print("WARNING: (%s:%d) %s\n",__FILE__, __LINE__,"step three: length of vect z null.");
      return(FALSE);
   }

   if (fabs(k)>fabs(j)) {
      rx = fabs(asin((double)(j/len)));
   }
   else {
      rx = fabs(acos((double)(k/len)));
   }

   if (j< 0) { 
      rx = -rx;
   }
        
   rx = -rx;  /* these are the required rotations */
   ry = -ry;
   rz = -rz;

   ang[1] = rx;
   ang[2] = ry;
   ang[3] = rz;

   VIO_FREE2D(t);
   VIO_FREE2D(s);
   VIO_FREE2D(R);
   VIO_FREE2D(Rx);
   VIO_FREE2D(Ry);
   VIO_FREE2D(Rz);

   return(TRUE);
}
