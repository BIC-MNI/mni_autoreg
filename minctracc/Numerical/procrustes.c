/* ----------------------------- MNI Header -----------------------------------
@NAME       : procrustes.c
@DESCRIPTION: File containing routines for doing procrustes calculations.
@METHOD     : Contains routines :
                 procrustes
                 transformations_to_homogeneous
                 translation_to_homogeneous
                 rotation_to_homogeneous
@CALLS      : 
@CREATED    : January 29, 1992 (Peter Neelin)
@MODIFIED   : February 7, 1992 (Peter Neelin)
                 - added routine transformations_to_homogeneous
---------------------------------------------------------------------------- */
#include <def_mni.h>
#include <recipes.h>
#include <matrix_basics.h>

/* Routines called in this file */
void svdcmp(float **, int, int, float *, float **);

/* Routines defined in this file */
public void procrustes(int npoints, int ndim, 
                       float **Apoints, float **Bpoints,
                       float *translation, float *centre_of_rotation,
                       float **rotation, float *scale);
public void transformations_to_homogeneous(int ndim, 
                  float *translation, float *centre_of_rotation,
                  float **rotation, float scale,
                  float **transformation);
public void translation_to_homogeneous(int ndim, float *translation,
                                       float **transformation);
public void rotation_to_homogeneous(int ndim, float **rotation,
                                       float **transformation);



/* ----------------------------- MNI Header -----------------------------------
@NAME       : procrustes
@INPUT      : npoints - number of input point pairs
              ndim    - number of dimensions for each point
              Apoints - Matrix of point set 1 (in numerical recipes
                 form). The dimensions of this matrix should be defined
                 to be 1 to npoints and 1 to ndim (when calling the numerical
                 recipes routine matrix).
              Bpoints - Matrix of point set 2 (in numerical recipes
                 form). The dimensions of this matrix should be defined
                 to be 1 to npoints and 1 to ndim (when calling the numerical
                 recipes routine matrix).
@OUTPUT     : translation - Numerical recipes vector (1 to ndim) that 
                 specifies the translation to be applied to Bpoints to line
                 up the centroid with that of Apoints. Calling routine must
                 allocate space for this vector.
              centre_of_rotation - Numerical recipes vector (1 to ndim) that
                 specifies the centre of rotation and scaling (this is 
                 in fact only the centroid of Apoints). Calling routine must
                 allocate space for this vector.
              rotation - Numerical recipes matrix (1 to ndim by 1 to ndim) 
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
@CALLS      : numerical recipes stuff
              calc_centroid
              translate
              transpose
              matrix_multiply
              svdcmp (numerical recipes)
              trace
@CREATED    : Long time ago (Sean Marrett)
@MODIFIED   : Some time later (Shyan Ku)
              Feb. 26, 1990 (Weiqian Dai)
              January 30, 1992 (Peter Neelin)
                 - complete rewrite for roughly NIL-abiding code. Modified
                 name and calling parameters.
---------------------------------------------------------------------------- */
public void procrustes(int npoints, int ndim, 
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



/* ----------------------------- MNI Header -----------------------------------
@NAME       : transformations_to_homogeneous
@INPUT      : ndim    - number of dimensions
              translation - Numerical recipes vector (1 to ndim) that 
                 specifies the translation to be applied first.
              centre_of_rotation - Numerical recipes vector (1 to ndim) that
                 specifies the centre of rotation and scaling.
              rotation - Numerical recipes matrix (1 to ndim by 1 to ndim) 
                 for rotation about centre_of_rotation (applied after 
                 translation). Note that this matrix need not only specify
                 rotation/reflexion - any ndim x ndim matrix will work.
              scale - Scalar value giving global scaling to be applied after
                 translation and rotation.
@OUTPUT     : transformation - Numerical recipes matrix (1 to ndim+1 by
                 1 to ndim+1) specifying the transformation for homogeneous 
                 coordinates. To apply this transformation, a point
                 vector should be post-multiplied by this matrix, with the
                 last coordinate of the ndim+1 point vector having value
                 one. The calling routine must allocate space for this
                 matrix.
@RETURNS    : (nothing)
@DESCRIPTION: Computes a transformation matrix in homogeneous coordinates
              given a translation, a rotation matrix (or other 
              non-homogeneous matrix) and a global scaling factor.
              Transformations are applied in that order.
@METHOD     : Apply the following operations (multiply from left to right):
                 1) Translate by translation
                 2) Translate by -centre_of_rotation
                 3) Rotate
                 4) Scale
                 5) Translate by centre_of_rotation
@GLOBALS    : (none)
@CALLS      : numerical recipes stuff
              translation_to_homogeneous
              matrix_multiply
              matrix_scalar_multiply
@CREATED    : February 7, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void transformations_to_homogeneous(int ndim, 
                  float *translation, float *centre_of_rotation,
                  float **rotation, float scale,
                  float **transformation)
{
   int i;
   int size;
   float *centre_translate;
   float **trans1, **trans2;
   float **trans_temp, **rotation_and_scale;

   size=ndim+1;

   /* Allocate matrices and vectors */
   centre_translate = vector(1,ndim);
   trans1 = matrix(1,size,1,size);
   trans2 = matrix(1,size,1,size);
   trans_temp = matrix(1,size,1,size);
   rotation_and_scale = matrix(1,ndim,1,ndim);


   /* Construct translation matrix */
   translation_to_homogeneous(ndim, translation, trans1);


   /* Construct translation matrix for centre of rotation and
      apply it */
   for (i=1; i<=ndim; i++) centre_translate[i] = -centre_of_rotation[i];
   translation_to_homogeneous(ndim, centre_translate, trans_temp);
   matrix_multiply(size, size, size, trans1, trans_temp, trans2);


   /* Scale rotation matrix, then convert it to homogeneous coordinates and
      apply it */
   matrix_scalar_multiply(ndim, ndim, scale, rotation, rotation_and_scale);
   rotation_to_homogeneous(ndim, rotation_and_scale, trans_temp);
   matrix_multiply(size, size, size, trans2, trans_temp, trans1);


   /* Return to centre of rotation */
   translation_to_homogeneous(ndim, centre_of_rotation, trans_temp);
   matrix_multiply(size, size, size, trans1, trans_temp, transformation);


   /* Free matrices */
   free_vector(centre_translate,1,ndim);
   free_matrix(trans1,1,size,1,size);
   free_matrix(trans2,1,size,1,size);
   free_matrix(trans_temp,1,size,1,size);
   free_matrix(rotation_and_scale,1,ndim,1,ndim);

}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : translation_to_homogeneous
@INPUT      : ndim    - number of dimensions
              translation - Numerical recipes vector (1 to ndim) that 
                 specifies the translation.
@OUTPUT     : transformation - Numerical recipes matrix (1 to ndim+1 by
                 1 to ndim+1) specifying the transformation for homogeneous 
                 coordinates. To apply this transformation, a point
                 vector should be post-multiplied by this matrix, with the
                 last coordinate of the ndim+1 point vector having value
                 one. The calling routine must allocate space for this
                 matrix.
@RETURNS    : (nothing)
@DESCRIPTION: Computes a transformation matrix in homogeneous coordinates
              given a translation.
@METHOD     : 
@GLOBALS    : (none)
@CALLS      : 
@CREATED    : February 7, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void translation_to_homogeneous(int ndim, float *translation,
                                       float **transformation)
{
   int i,j;
   int size;

   size=ndim+1;

   /* Construct translation matrix */
   for (i=1; i<=ndim; i++) {
      for (j=1; j<=size; j++) {
         if (i == j) {
            transformation[i][j] = 1.0;
         }
         else {
            transformation[i][j] = 0.0;
         }
      }
   }
   for (j=1; j<=ndim; j++) {
      transformation[size][j] = translation[j];
   }

   transformation[size][size] = 1.0;

}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : rotation_to_homogeneous
@INPUT      : ndim    - number of dimensions
              rotation - Numerical recipes matrix (1 to ndim by 1 to ndim) 
                 for rotation about origin. Note that this matrix need not 
                 only specify rotation/reflexion - any ndim x ndim matrix 
                 will work.
@OUTPUT     : transformation - Numerical recipes matrix (1 to ndim+1 by
                 1 to ndim+1) specifying the transformation for homogeneous 
                 coordinates. To apply this transformation, a point
                 vector should be post-multiplied by this matrix, with the
                 last coordinate of the ndim+1 point vector having value
                 one. The calling routine must allocate space for this
                 matrix.
@RETURNS    : (nothing)
@DESCRIPTION: Computes a transformation matrix in homogeneous coordinates
              given a rotation matrix.
@METHOD     : 
@GLOBALS    : (none)
@CALLS      : 
@CREATED    : February 7, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void rotation_to_homogeneous(int ndim, float **rotation,
                                       float **transformation)
{
   int i,j;
   int size;

   size=ndim+1;

   /* Construct  matrix */
   for (i=1; i<=size; i++) {
      for (j=1; j<=size; j++) {
         if ((i==size) || (j==size)) {
            transformation[i][j] = 0.0;
         }
         else {
            transformation[i][j] = rotation[i][j];
         }
      }
   }

   transformation[size][size] = 1.0;

}
