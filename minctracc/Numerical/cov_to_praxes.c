/* ----------------------------- MNI Header -----------------------------------
@NAME       : cov_to_praxes
@INPUT      : ndim    - number of dimensions (rank of covariance matrix)
              covar   - covariance matrix (in zero offset form)
@OUTPUT     : pr_axes - principal axes (matrix in zero offset form).
                        Note that indices of vectors vary fastest, ie.
                        first vector is pr_axes[1][i], where i is the vector 
                        index ranging from 1 to ndim. Vector length is equal
                        to the standard deviation along the axis. The axes
                        are sorted into descending order of length.
@RETURNS    : nothing
@DESCRIPTION: Calculates the principal axes of an n-dimensional covariance
              matrix.
@METHOD     : Calculates the matrix of eigenvectors and the eigenvalues.
@GLOBALS    : none
@CALLS      : 
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

@CREATED    : July 16, 1991 (Peter Neelin)
                      using numerical recipes routines jacobi() and eigsrt().  
		      See Hotelling Transform
@MODIFIED   : $Log: cov_to_praxes.c,v $
@MODIFIED   : Revision 1.8  1995-09-11 12:32:29  louis
@MODIFIED   : all refs to numerical recipes routines have been removed. The procedure
@MODIFIED   : eigen() replaces calls to jacobi() and eigsrt(). eigen calculates the
@MODIFIED   : eigen vectors and eigen values of a symmetric matrix, and returns the
@MODIFIED   : eigen vectors (sorted by eigen value) in columns of the eigen_vec matrix.
@MODIFIED   :

   the current bit of code for eigen originally came from Pascal @
   Rennes for the analysis of 37 element meg data.  It has been
   modified by LC to replace jacobi() and eigsrt() originally used in
   cov_to_praxes.
 
---------------------------------------------------------------------------- */
#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Numerical/cov_to_praxes.c,v 1.8 1995-09-11 12:32:29 louis Exp $";
#endif

#include <volume_io.h>

public void eigen(double **inputMat, 
		  int    ndim, 
		  double *eigen_val, 
		  double **eigen_vec, 
		  int    *iters)

{

  double
    vec[37], test_vec[37],
    x, 
    vec_sum, previous_sum, 
    val,     previous_val,
    max_val;
  int 
    i,j,k,iterations;

  
  *iters = 0;

  for_less(k,0,ndim) {

    for_less(i,0,ndim)		/* init vec to ones */
      vec[i]=1;

    iterations=0;
    val=0.0;
    vec_sum=(double)ndim;	/* since  sum += vec[i]  */


    do {			/* interate  */

      iterations++;
      previous_val=val;
      previous_sum=vec_sum;
      max_val=0.0;
      vec_sum=0.0;

      for_less(i,0,ndim) {	/* for each row in the inputMat */

	test_vec[i]=0.0;	/* test_vec = inputMat * vec'   */
	for_less(j,0,ndim)	  
	  test_vec[i] += inputMat[i][j]*vec[j];

	vec_sum += test_vec[i];

	x = ABS(test_vec[i]);	/* find max element in test_vec */
	if(x>max_val) max_val=x;

      }

      val     = vec_sum/previous_sum;
      vec_sum = vec_sum/max_val;

      for_less(i,0,ndim)	/* nomalize */
	vec[i]=test_vec[i]/max_val;

    } while((iterations<3000)&&(fabs(1-previous_val/val)>1e-15));
    
    *iters += iterations;

    if((fabs(val)>1e-6)&&(iterations>2900))	
      printf(" Possible numerical error(%d iter, %12.9f)\n",iterations,fabs(val));
    
    vec_sum=0.0;
    for_less(i,0,ndim) 
      vec_sum+=vec[i]*vec[i];

    for_less(i,0,ndim)
      for_less(j,0,ndim)
	inputMat[i][j]-=vec[i]*vec[j]/vec_sum*val;

    for_less(i,0,ndim)
      eigen_vec[k][i]=vec[i]/sqrt(vec_sum); /* store eig vec in row 'k' */

    eigen_val[k]=val;
  }
				/* sort eig vectors in order of eigval */
  for_less(i,0,ndim) {
    max_val=eigen_val[k=i];
    for_less (j,i+1,ndim)
      if (eigen_val[j] >= max_val) max_val=eigen_val[k=j];

    if (k != i) {
      eigen_val[k]=eigen_val[i]; /* swap eigen values */
      eigen_val[i]=max_val;

      for_less (j,0,ndim) {	/* swap eigen vectors (rows of eigen_vec)*/
	val=            eigen_vec[i][j];
	eigen_vec[i][j]=eigen_vec[k][j];
	eigen_vec[k][j]=val;
      }

    }

  }
				/* transpose eigen vector (so that
                                   each eigen-vector is in a column,
                                   so that the result is similar to
                                   the old call to jacobi() and
                                   eigsrt() */
  for_less(i,0,ndim) {
    for_less (j,i+1,ndim) {
      	val=eigen_vec[i][j];
	eigen_vec[i][j]=eigen_vec[j][i];
	eigen_vec[j][i]=val;
    }
  }

}


void cov_to_praxes(int ndim, float **covar, float **pr_axes)
{
  double **amat,*eigval,**eigvec;
  int nrot,i,j;
  
  /* Set up the matrices, copying covar since jacobi modifies its
     input matrix */
  ALLOC2D(amat,ndim+1,ndim+1);
  ALLOC2D(eigvec,ndim+1,ndim+1);
  ALLOC(eigval,ndim+1);
  
				/* copy the input matrix, since
				   eigen() modifies it */
  for (i=1; i<=ndim; i++) {
    for (j=1; j<=ndim; j++) {
      amat[i-1][j-1]=covar[i][j];
    }
  }
  
  /* Get the sorted eigen values and eigen vectors  */
  eigen(amat, ndim, eigval, eigvec, &nrot);
  
  /* Get the principal axes /* each eigen vector 'j' is in column 'j' */
  for (i=1; i<=ndim; i++) {
    for (j=1; j<=ndim; j++) {
      pr_axes[i][j]=sqrt(fabs((double)eigval[i-1]))*eigvec[i-1][j-1];
    }
  }
  
  /* Free up the matrices */
  FREE2D(amat);
  FREE2D(eigvec);
  FREE(eigval);
  
}






