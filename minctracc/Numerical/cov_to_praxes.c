#include <stdio.h>
#include <math.h>
#include <recipes.h>

/* ----------------------------- MNI Header -----------------------------------
@NAME       : cov_to_praxes
@INPUT      : ndim    - number of dimensions (rank of covariance matrix)
              covar   - covariance matrix (in numerical recipes form)
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
              See ..., 
              Section 3.6 "The Hotelling Transform", pp. 122 - 125.
@GLOBALS    : none
@CALLS      : jacobi, eigsrt, matrix, vector, free_matrix, free_vector,
              fabs, sqrt
@CREATED    : July 16, 1991 (Peter Neelin)
@MODIFIED   : September 30, 1991 (P.N.)
                 - copied from calc_princ_axes and modified so that calling
                   parameters are closer to numerical recipes calls, and so
                   that code adheres to the NIL coding standard.
---------------------------------------------------------------------------- */
void cov_to_praxes(int ndim, float **covar, float **pr_axes)
{
   float **amat,*eigval,**eigvec;
   int nrot,i,j;

   /* Set up the matrices, copying covar since jacobi modifies its
      input matrix */
   amat=matrix(1,ndim,1,ndim);
   eigvec=matrix(1,ndim,1,ndim);
   eigval=vector(1,ndim);
   for (i=1; i<=ndim; i++) {
      for (j=1; j<=ndim; j++) {
         amat[i][j]=covar[i][j];
      }
   }

   /* Get the eigenvector matrix */
   jacobi(amat, ndim, eigval, eigvec, &nrot);
   eigsrt(eigval, eigvec, ndim);

   /* Get the principal axes */
   for (i=1; i<=ndim; i++) {
      for (j=1; j<=ndim; j++) {
         pr_axes[i][j]=sqrt(fabs((double)eigval[i]))*eigvec[i][j];
      }
   }

   /* Free up the matrices */
   (void) free_matrix(amat,1,ndim,1,ndim);
   (void) free_matrix(eigvec,1,ndim,1,ndim);
   (void) free_vector(eigval,1,ndim);

}
