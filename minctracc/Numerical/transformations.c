/* ----------------------------- MNI Header -----------------------------------
@NAME       : transformations.c
@DESCRIPTION: File containing routines to transform points in world
              space using linear or non-linear transformations.
@METHOD     : 
@GLOBALS    : 
@CREATED    : Wed May 26 13:05:44 EST 1993 LC (from NEELIN)
@MODIFIED   : Tue Jun  1 09:16:46 EST 1993 
                 separated the interpolation stuff from the transformation
		 routines.


---------------------------------------------------------------------------- */


#include <def_mni.h>
#include <recipes.h>
#include "minctracc.h"

				/* Some external functions used in this file */

void lubksb(float **a, int n, int *indx, float *b);
void ludcmp(float **a, int n, int *indx, float *d);

				/* prototypes for this file: */

private float return_r(double *cor1, double *cor2, int dim);

private float FU(double r, int dim);

private void mappingf(double **bdefor, double **INVMLY, int num_marks, 
		      double *icor, double *rcor, int dim);


/* ----------------------------- MNI Header -----------------------------------
@NAME       : do_linear_transformation
@INPUT      : trans_data - pointer to transformation data
              coordinate - point to be transformed
@OUTPUT     : result - resulting coordinate
@RETURNS    : (nothing)
@DESCRIPTION: Routine to apply a linear transformation. 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : February 10, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void do_linear_transformation(Coord_Vector *result, void *trans_data, 
                                     Coord_Vector *coordinate)
{
   Linear_Transformation *matrx;
   int idim, jdim;
   double lcoord[WORLD_NDIMS], lresult[WORLD_NDIMS];

   /* Get linear transformation info */
   matrx = trans_data;

   /* Make our own coord vector */
   lcoord[X] = coordinate->x;
   lcoord[Y] = coordinate->y;
   lcoord[Z] = coordinate->z;

   /* Calculate transformation */
   for (idim=0; idim<WORLD_NDIMS; idim++) {
      lresult[idim] = matrx->mat[idim][MAT_NDIMS-1];
      for (jdim=0; jdim<WORLD_NDIMS; jdim++) {
         lresult[idim] += matrx->mat[idim][jdim] * lcoord[jdim];
      }
   }

   /* Save the result */
   result->x = lresult[X];
   result->y = lresult[Y];
   result->z = lresult[Z];

   return;
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : do_non_linear_transformation
@INPUT      : trans_data - pointer to transformation data
              coordinate - point to be transformed
@OUTPUT     : result - resulting coordinate
@RETURNS    : (nothing)
@DESCRIPTION: Routine to apply a linear transformation. 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : February 10, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void do_non_linear_transformation(Coord_Vector *result, void *trans_data, 
                                     Coord_Vector *coordinate)
{
  Thin_plate_spline *tps;
  double lcoord[WORLD_NDIMS], lresult[WORLD_NDIMS];

  /* Get non-linear transformation info */
  tps = trans_data;

  /* Make our own coord vector */
  lcoord[X] = coordinate->x;
  lcoord[Y] = coordinate->y;
  lcoord[Z] = coordinate->z;

  mappingf(tps->coords, tps->warps, tps->num_points, lcoord, lresult, tps->dim);
  
  /* Save the result */
  result->x = lresult[X];
  result->y = lresult[Y];
  result->z = lresult[Z];
  
  return;
}



/* ----------------------------- MNI Header -----------------------------------
@NAME       : invert_transformation
@INPUT      : transformation - transformation to invert
@OUTPUT     : result - resultant transformation
@RETURNS    : (nothing)
@DESCRIPTION: Routine to invert a transformation. Currently only works on
              linear transformations.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : February 9, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void invert_transformation(Transformation *result, 
                                  Transformation *transformation)
{
   Linear_Transformation *matrx, *result_matrix;
   float **nrmatrix, **nrresult, nrd, nrcol[MAT_NDIMS+1];
   int nrindex[MAT_NDIMS+1], idim, jdim;

   /* Check that transformation is linear */
   if (!IS_LINEAR(transformation)) {
      (void) fprintf(stderr, "Unable to invert non-linear transformations!\n");
      exit(EXIT_FAILURE);
   }
   matrx = transformation->trans_data;

   /* Set up numerical recipes matrices */
   nrmatrix = matrix(1, MAT_NDIMS, 1, MAT_NDIMS);
   nrresult = matrix(1, MAT_NDIMS, 1, MAT_NDIMS);

   for (idim=1; idim<=MAT_NDIMS; idim++) {
      for (jdim=1; jdim<=MAT_NDIMS; jdim++) {
         if (idim<=WORLD_NDIMS)
            nrmatrix[idim][jdim] = matrx->mat[idim-1][jdim-1];
         else if (jdim<=WORLD_NDIMS)
            nrmatrix[idim][jdim] = 0.0;
         else 
            nrmatrix[idim][jdim] = 1.0;
      }
   }

   /* Invert matrix */
   ludcmp( nrmatrix, MAT_NDIMS, nrindex, &nrd );

   for (jdim=1; jdim<=MAT_NDIMS; jdim++) {
      for (idim=1; idim<=MAT_NDIMS; idim++)
         nrcol[idim] = 0.0;
      nrcol[jdim] = 1.0;
      lubksb( nrmatrix, MAT_NDIMS, nrindex, nrcol );
      for (idim=1; idim<=MAT_NDIMS; idim++)
         nrresult[idim][jdim] = nrcol[idim];
   }

   /* Save the result */
   if (result->trans_data != NULL) FREE(result->trans_data);
   *result = *transformation;
   ALLOC ( matrx , 1);  
   result->trans_data = matrx;
   for (idim=1; idim<=WORLD_NDIMS; idim++)
      for (jdim=1; jdim<=MAT_NDIMS; jdim++)
         matrx->mat[idim-1][jdim-1] = nrresult[idim][jdim];

   /* free the nr matrices */
   free_matrix(nrmatrix, 1, MAT_NDIMS, 1, MAT_NDIMS);
   free_matrix(nrresult, 1, MAT_NDIMS, 1, MAT_NDIMS);


   return;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : mult_linear_transform
@INPUT      : transform1 - first transformation
              transform2 - second transformation
@OUTPUT     : result - resulting transformation
@RETURNS    : (nothing)
@DESCRIPTION: Routine to multiply two linear matrices.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : February 10, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void mult_linear_transform(Transformation *result, 
                                  Transformation *transform1, 
                                  Transformation *transform2)
{
   int idim, jdim, kdim;
   Linear_Transformation *matrix1, *matrix2, *result_matrix;

   /* Check for linear transformations */
   if (!IS_LINEAR(transform1) || !IS_LINEAR(transform2)) {
      (void) fprintf(stderr, 
                     "Unable to multiply two non-linear transformations.\n");
      exit(EXIT_FAILURE);
   }
   matrix1 = transform1->trans_data;
   matrix2 = transform2->trans_data;

   /* Get space for the result */
   ALLOC( result_matrix, 1 );

   /* Multiply the matrices */
   for (idim=0; idim < WORLD_NDIMS; idim++) {
      for (jdim=0; jdim < MAT_NDIMS; jdim++) {
         if (jdim < WORLD_NDIMS) 
            result_matrix->mat[idim][jdim] = 0.0;
         else
            result_matrix->mat[idim][jdim] = 
               matrix1->mat[idim][MAT_NDIMS-1];
         for (kdim=0; kdim < WORLD_NDIMS; kdim++) {
            result_matrix->mat[idim][jdim] += 
               matrix1->mat[idim][kdim] * matrix2->mat[kdim][jdim];
         }
      }
   }

   /* Save the result */
   if (result->trans_data != NULL) FREE(result->trans_data);
   *result = *transform1;
   result->trans_data = result_matrix;

   return;
}


 

/* ----------------------------- MNI Header -----------------------------------
@NAME       : mappingf - thin plate spline mapping.
@INPUT      : bdefor - numerical recipes array with 'num_marks' rows, and 'dim' cols,
                       contains the list of landmarks points in the 'source' volume
	      INVMLY - numerical recipes array with 'num_marks+dim+1' rows, and 'dim' cols,
	               contains the deformation vectors that define the thin plate spline
	      num_marks - number of landmark points
              icor   - array of float [1..dim] for the input coordinate in the
                       'source' volume space.
	      dim    - number of dimensions (either 2 or 3).
@OUTPUT     : rcor   - array of float [1..dim] of the output coordinate in the 
                       'target' volume.
@RETURNS    : nothing
@DESCRIPTION: 
              INVMLY to 'in' to get 'out'
@METHOD     : 
@GLOBALS    : none
@CALLS      : return_r, vector, free_vector
@CREATED    : Mon Apr  5 09:00:54 EST 1993
@MODIFIED   : 
---------------------------------------------------------------------------- */
private void mappingf(double **bdefor, double **INVMLY, int num_marks, 
	      double *icor, double *rcor, int dim)
{

  /*  note that :
      INVMLY is defined from [0..num_marks+dim][0..dim-1]
      bdefor is defined from [0..num_marks-1][0..dim-1]
      icor[0..dim-1]
      rcor[0..dim-1]
      */


  int i,j,markpoint;
  double dist,weight,tempcor[WORLD_NDIMS];
  
  
  /* f(x,y) =a_{1} + a_{x}x + a_{y}y + sum_{1}^{n}
   *          w_{i}U(|P_{i} - (x,y)|) 
   */
  
  for (i=0;i<WORLD_NDIMS;i++){
    tempcor[i] = 0;
  }

   for (i=0; i<num_marks; i++){

      markpoint = 1;		/* set the point as the landmark point 
				 check the point where is landmark */

      for (j=0;(j<dim)&&(markpoint==1);j++){
	if ((bdefor[i][j] != icor[j])){
	  markpoint = 0;
	}
      }

      /* because the point at the landmark is not deformable, we use 
       * the close point to the landmark.
       */

      if (markpoint == 1){
         icor[1] = icor[1] + 0.00001;
      }

      dist = return_r(bdefor[i],icor,dim); 
      weight = FU(dist,dim);

      for (j=0;j<dim;j++){
         tempcor[j] = INVMLY[i][j]*weight + tempcor[j];
      }
   } 

   for (j=0;j<dim;j++){
      rcor[j] = INVMLY[num_marks][j]+tempcor[j];
   }
   

   for (j = 0; j<dim; j++){   
      for (i = 0; i<dim; i++){
         rcor[j] = INVMLY[num_marks+i+1][j]*icor[i] + rcor[j];
      }
   }
}

private float return_r(double *cor1, double *cor2, int dim)
{
   double r1,r2,r3;


   if (dim == 1){
      r1 = cor1[0] - cor2[0]; 
      r1 = fabs(r1);
      return(r1);
   } else
   if (dim == 2){
      r1 = cor1[0] - cor2[0];
      r2 = cor1[1] - cor2[1];
      return(r1*r1+r2*r2);
   } else
   if (dim == 3){
      r1 = cor1[0] - cor2[0];
      r2 = cor1[1] - cor2[1];
      r3 = cor1[2] - cor2[2];
      return(fsqrt(r1*r1 + r2*r2 + r3*r3));
   } 
   else { printf(" impossible error in mapping.c (return_r)\n"); exit(-1); }
}

private float FU(double r, int dim)
/* This function will calculate the U(r) function.
 * if dim = 1, the function returns |r|^3 
 * if dim = 2, the function retruns r^2*log r^2
 * if dim = 3, the function returns |r|
 */ 
{
   double z;

   if (dim==1){
      z = r*r*r;
      return(fabs(z));
   } else
   if (dim==2){/* r is stored as r^2 */
      z = r * log(r);
      return(z);
   } else
   if (dim==3){
      return(fabs(r));
   } else { printf(" impossible error in mapping.c (FU)\n"); exit(-1); }
}
