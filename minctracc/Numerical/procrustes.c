/* ----------------------------- MNI Header -----------------------------------
@NAME       : procrustes.c
@DESCRIPTION: File containing routines for doing procrustes calculations.
@METHOD     : Contains routines :
                 procrustes
@CALLS      : 
@COPYRIGHT  :
              Copyright 1993 Peter Neelin, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.

@CREATED    : January 29, 1992 (Peter Neelin)
@MODIFIED   :  $Log: procrustes.c,v $
@MODIFIED   :  Revision 1.10  2004-02-12 16:22:43  louis
@MODIFIED   :   * removed public/private defs
@MODIFIED   :
@MODIFIED   :  Revision 1.9  1995/09/11 12:37:16  louis
@MODIFIED   :  this file is not used in mni_reg-0.1g, and some numerical recipes
@MODIFIED   :  routines remain.
@MODIFIED   :
 * Revision 1.8  1995/02/22  08:56:06  louis
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.7  94/04/06  11:48:47  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.6  94/02/21  16:36:20  louis
 * version before feb 22 changes
 * 
 * Revision 1.5  93/11/15  16:27:10  louis
 * working version, with new library, with RCS revision stuff,
 * before deformations included
 * 

February 7, 1992 (Peter Neelin)
                 - added routine transformations_to_homogeneous
Fri Jun  4 14:10:34 EST 1993 LC

removed
    transformations_to_homogeneous
    translation_to_homogeneous
    rotation_to_homogeneous
and moved them to matrix_basics!
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Numerical/procrustes.c,v 1.10 2004-02-12 16:22:43 louis Exp $";
#endif

#include <volume_io.h>
#include "matrix_basics.h"

/* Routines called in this file */
void svdcmp(float **, int, int, float *, float **);

/* Routines defined in this file */
void procrustes(int npoints, int ndim, 
                       float **Apoints, float **Bpoints,
                       float *translation, float *centre_of_rotation,
                       float **rotation, float *scale);


/* ----------------------------- MNI Header -----------------------------------
@NAME       : procrustes
@INPUT      : npoints - number of input point pairs
              ndim    - number of dimensions for each point
              Apoints - Matrix of point set 1 (in zero offset
                 form). The dimensions of this matrix should be defined
                 to be 1 to npoints and 1 to ndim (when calling the numerical
                 recipes routine matrix).
              Bpoints - Matrix of point set 2 (in zero offset
                 form). The dimensions of this matrix should be defined
                 to be 1 to npoints and 1 to ndim (when calling the numerical
                 recipes routine matrix).
@OUTPUT     : translation - zero offset vector (1 to ndim) that 
                 specifies the translation to be applied to Bpoints to line
                 up the centroid with that of Apoints. Calling routine must
                 allocate space for this vector.
              centre_of_rotation - zero offset vector (1 to ndim) that
                 specifies the centre of rotation and scaling (this is 
                 in fact only the centroid of Apoints). Calling routine must
                 allocate space for this vector.
              rotation - zero offset matrix (1 to ndim by 1 to ndim) 
                 to rotate translated Bpoints so that they line up with 
                 Apoints. Calling routine must allocate space for this 
                 matrix.
              scale - Scalar value giving global scaling to be applied to
                 translated and rotated Bpoints to match Apoints.
@RETURNS    : (nothing)
@DESCRIPTION: Calculates n-dimensional linear transformation from one set 
              of points to another, minimizing distance between equivalent
              points. Transformation from Bpoints to Apoints is calculated.
@METHOD     : See Matrix Computations, Golub and Van Loan, pp. 425-426 and
              paper by Sibson, Robin, J.R.Statist.Soc. B(1978), Vol. 40,
              No. 2, pp 234-238.
              Steps of calculations are as follows :
                 1) Calculate translation that aligns the centroids of the
                    two point sets.
                 2) Calculate rotation/reflexion that minimizes residual.
                 3) Calculate scaling of points to minimize residual.
              The process can be broken into independent steps because the
              best translation aligns centroids independently of the choice
              of rotation/reflexion and scaling and the best rotation/reflexion
              can be found independently of scale (after the best translation
              has been found). (See Sibson for more).
@GLOBALS    : (none)
@CALLS      : calc_centroid
              translate
              transpose
              matrix_multiply
              svdcmp (zero offset)
              trace
@CREATED    : Long time ago (Sean Marrett)
@MODIFIED   : Some time later (Shyan Ku)
              Feb. 26, 1990 (Weiqian Dai)
              January 30, 1992 (Peter Neelin)
                 - complete rewrite for roughly NIL-abiding code. Modified
                 name and calling parameters.
---------------------------------------------------------------------------- */
void procrustes(int npoints, int ndim, 
                       float **Apoints, float **Bpoints,
                       float *translation, float *centre_of_rotation,
                       float **rotation, float *scale)
{
   int i;
   float *Atranslation, *Btranslation, *svd_W;
   float **Ashift, **Bshift, **Atranspose, **Btranspose;
   float **svd_U, **svd_V, **svd_VT;
   float **Brotated, **product;
   float trace1, trace2;
                                   
   /* Get the vectors for centroids */
   Atranslation=vector(1,ndim);
   Btranslation=vector(1,ndim);
   svd_W=vector(1,ndim);

   /* Get various matrices */
   Ashift=matrix(1,npoints,1,ndim);
   Bshift=matrix(1,npoints,1,ndim);
   Atranspose=matrix(1,ndim,1,npoints);
   Btranspose=matrix(1,ndim,1,npoints);
   svd_U=matrix(1,ndim,1,ndim);
   svd_V=matrix(1,ndim,1,ndim);
   svd_VT=matrix(1,ndim,1,ndim);
   Brotated=matrix(1,npoints,1,ndim);
   product=matrix(1,npoints,1,npoints);

   /* Calculate the centroids, remove them from A and B points and
    save the translation */

   calc_centroid(npoints, ndim, Apoints, centre_of_rotation); 
   for (i=1; i<=ndim; i++) Atranslation[i] = -centre_of_rotation[i];
   translate(npoints, ndim, Apoints, Atranslation, Ashift);
   calc_centroid(npoints, ndim, Bpoints, Btranslation); 
   for (i=1; i<=ndim; i++) Btranslation[i] *= -1;
   translate(npoints, ndim, Bpoints, Btranslation, Bshift);

   for (i=1; i<=ndim; i++) translation[i] = Btranslation[i] - Atranslation[i];


   /* Calculate the rotation/reflexion matrix */

   transpose(npoints, ndim, Bshift, Btranspose);
   matrix_multiply(ndim, npoints, ndim, Btranspose, Ashift, svd_U);
   svdcmp(svd_U, ndim, ndim, svd_W, svd_V);
   transpose(ndim, ndim, svd_V, svd_VT);
   matrix_multiply(ndim, ndim, ndim, svd_U, svd_VT, rotation);


   /* Calculate the scale */

   matrix_multiply(npoints, ndim, ndim, Bshift, rotation, Brotated);
   transpose(npoints, ndim, Ashift, Atranspose);
   matrix_multiply(npoints, ndim, npoints, Brotated, Atranspose, product);
   trace1 = trace(npoints, product);
   matrix_multiply(npoints, ndim, npoints, Bshift, Btranspose, product);
   trace2 = trace(npoints, product);
   if (trace2 != 0.0) {
      *scale = trace1 / trace2;
   }
   else {
      *scale = 0.0;
   }


   /* transpose back the rotation matrix */

   transpose(ndim, ndim, rotation, rotation);

   /* Free vectors */
   free_vector(Atranslation,1,ndim);
   free_vector(Btranslation,1,ndim);
   free_vector(svd_W,1,ndim);

   /* Free matrices */
   free_matrix(Ashift,1,npoints,1,ndim);
   free_matrix(Bshift,1,npoints,1,ndim);
   free_matrix(Atranspose,1,ndim,1,npoints);
   free_matrix(Btranspose,1,ndim,1,npoints);
   free_matrix(svd_U,1,ndim,1,ndim);
   free_matrix(svd_V,1,ndim,1,ndim);
   free_matrix(svd_VT,1,ndim,1,ndim);
   free_matrix(Brotated,1,npoints,1,ndim);
   free_matrix(product,1,npoints,1,npoints);
}

