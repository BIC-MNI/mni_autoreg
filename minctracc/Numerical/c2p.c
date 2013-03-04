#include <def_mni.h>
#include <recipes.h>

/* Routine defined in this file */
void cov_to_praxes(int ndim, float **covar, float **pr_axes);


/* ----------------------------- MNI Header -----------------------------------
@NAME       : cov_to_praxes
@INPUT      : ndim    - number of dimensions (rank of covariance matrix)
              covar   - covariance matrix (in numerical recipes form).
                        This matrix is partly destroyed by the routine.
                        The dimensions of this matrix (and of pr_axes)
                        should be defined to be 1 to ndim (when calling
                        the numerical recipes routine matrix).
@OUTPUT     : pr_axes - principal axes (matrix in numerical recipes form).
                        Note that indices of vectors vary fastest, ie.
                        first vector is pr_axes[1][i], where i is the vector 
                        index ranging from 1 to ndim. Vector length is equal
                        to the standard deviation along the axis. The axes
                        are sorted into descending order of length.
@RETURNS    : nothing
@DESCRIPTION: Calculates the principal axes of an n-dimensional covariance
              matrix.
@METHOD     : Calculates the matrix of eigenvectors and the eigenvalues.
              See Gonzalez and Wintz, Digital Image Processing, 2nd ed. 
              Section 3.6 "The Hotelling Transform", pp. 122 - 125.
@GLOBALS    : none
@CALLS      : jacobi, eigsrt, matrix, vector, free_matrix, free_vector,
              fabs, sqrt
@CREATED    : July 16, 1991 (Peter Neelin)
@MODIFIED   : September 30, 1991 (P.N.)
                 - copied from calc_princ_axes and modified so that calling
                   parameters are closer to numerical recipes calls, and so
                   that code adheres to the NIL coding standard.
              October 4, 1991 (P.N.)
                 - modified so that covar matrix is not copied (so it is
                   changed by jacobi).
---------------------------------------------------------------------------- */
void cov_to_praxes(int ndim, float **covar, float **pr_axes)
{
   float *eigval;
   int nrot,i,j;

   /* Get vector for eigenvalues */
   eigval=vector(1,ndim);

   /* Get the eigenvector matrix */
   jacobi(covar, ndim, eigval, pr_axes, &nrot);
   eigsrt(eigval, pr_axes, ndim);

   /* Get the principal axes */
   for (i=1; i<=ndim; i++) {
      for (j=1; j<=ndim; j++) {
         pr_axes[i][j] *= sqrt(fabs((double)eigval[i]));
      }
   }

   /* Free up the vector */
   (void) free_vector(eigval,1,ndim);

}
