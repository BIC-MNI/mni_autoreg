/* ----------------------------- MNI Header -----------------------------------
@NAME       : cov_to_praxes
@INPUT      : ndim    - number of dimensions (rank of covariance matrix)
              covar   - covariance matrix (in zero offset form)
@OUTPUT     : pr_axes - principal axes (matrix in zero offset form).
                        Note that indices of vectors vary fastest, ie.
                        first vector is pr_axes[1][i], where i is the vector 
                        index ranging from 1 to ndim. VIO_Vector length is equal
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
                      See Hotelling VIO_Transform
@MODIFIED   : $Log: cov_to_praxes.c,v $
@MODIFIED   : Revision 96.7  2006-11-30 09:07:32  rotor
@MODIFIED   :  * many more changes for clean minc 2.0 build
@MODIFIED   :
@MODIFIED   : Revision 96.6  2006/11/29 09:09:33  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 96.5  2005/07/20 20:45:48  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 96.4  2004/02/12 05:54:27  rotor
@MODIFIED   :  * removed /static defs
@MODIFIED   :
@MODIFIED   : Revision 96.3  2002/03/26 14:15:39  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 96.2  2000/03/15 08:42:43  stever
@MODIFIED   : Code cleanup: all functions prototyped (except ParseArgs.c), no useless declarations, etc
@MODIFIED   :
@MODIFIED   : Revision 96.1  1997/11/03 19:59:49  louis
@MODIFIED   : now include volume_io/internal_volume_io.h instead of volume_io.h
@MODIFIED   :
 * Revision 96.0  1996/08/21  18:21:58  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:21:52  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:15:50  louis
 * Never released version 0.95
 *
 * Revision 1.9  1996/08/12  14:15:43  louis
 * Pre-release
 *
 * Revision 1.8  1995/09/11  12:32:29  collins
 * all refs to numerical recipes routines have been removed. The procedure
 * eigen() replaces calls to jacobi() and eigsrt(). eigen calculates the
 * eigen vectors and eigen values of a symmetric matrix, and returns the
 * eigen vectors (sorted by eigen value) in columns of the eigen_vec matrix.
 *
 * Revision 1.8  1995/09/11  12:32:29  collins
 * all refs to numerical recipes routines have been removed. The procedure
 * eigen() replaces calls to jacobi() and eigsrt(). eigen calculates the
 * eigen vectors and eigen values of a symmetric matrix, and returns the
 * eigen vectors (sorted by eigen value) in columns of the eigen_vec matrix.
 *

---------------------------------------------------------------------------- */
#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/minctracc/Numerical/cov_to_praxes.c,v 96.7 2006-11-30 09:07:32 rotor Exp $";
#endif

#include <volume_io.h>

VIO_BOOL eigen(double **inputMat, int ndim, 
                     double *eigen_val, double **eigen_vec, 
                     int    *iters);


void cov_to_praxes(int ndim, float **covar, float **pr_axes)
{
  double **amat,*eigval,**eigvec;
  int nrot,i,j;
  
  VIO_ALLOC2D(amat,ndim+1,ndim+1);
  VIO_ALLOC2D(eigvec,ndim+1,ndim+1);
  ALLOC(eigval,ndim+1);

  /* copy the input matrix and transpose it, since eigen() modifies
     it, and eigen returns 'right-eigenvectors' ie 
          amat * eigvec = eigvec * eigval */

  for (i=1; i<=ndim; i++) {
    for (j=1; j<=ndim; j++) {
      amat[i-1][j-1]=covar[i][j];
    }
  }
  
  /* Get the sorted eigen values and eigen vectors, while ignoring the
     return value from eigen */

  (void)eigen(amat, ndim, eigval, eigvec, &nrot);
  
  /* Get the principal axes: 
     each eigen vector 'j' (which is in  column 'j') */

  for (i=1; i<=ndim; i++) {
    for (j=1; j<=ndim; j++) {
      pr_axes[j][i]=sqrt(fabs((double)eigval[i-1]))*eigvec[j-1][i-1];
    }
  }
  
  /* Free up the matrices */
  VIO_FREE2D(amat);
  VIO_FREE2D(eigvec);
  FREE(eigval);
  
}



#include <math.h>


#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
        a[k][l]=h+s*(g-h*tau);

static VIO_BOOL jacobi(double **a,
                       int n,
                       double *d,
                       double **v,
                       int *nrot)
{
  int j,iq,ip,i;
  double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
  
  ALLOC(b,n+1);
  ALLOC(z,n+1);
  for (ip=1;ip<=n;ip++) {
    for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
  for (ip=1;ip<=n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=1;i<=150;i++) {
    sm=0.0;
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++)
        sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0) {
      FREE(z);
      FREE(b);
      return(TRUE);                /* we have a result! */
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++) {
        g=100.0*fabs(a[ip][iq]);
        if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
            && (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
          a[ip][iq]=0.0;
        else {
          if (fabs(a[ip][iq]) > tresh) {
            h=d[iq]-d[ip];
            if ((double)(fabs(h)+g) == (double)fabs(h))
              t=(a[ip][iq])/h;
            else {
              theta=0.5*h/(a[ip][iq]);
              t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
              if (theta < 0.0) t = -t;
            }
            c=1.0/sqrt(1+t*t);
            s=t*c;
            tau=s/(1.0+c);
            h=t*a[ip][iq];
            z[ip] -= h;
            z[iq] += h;
            d[ip] -= h;
            d[iq] += h;
            a[ip][iq]=0.0;
            for (j=1;j<=ip-1;j++) {
              ROTATE(a,j,ip,j,iq)
            }
            for (j=ip+1;j<=iq-1;j++) {
              ROTATE(a,ip,j,j,iq)
            }
            for (j=iq+1;j<=n;j++) {
              ROTATE(a,ip,j,iq,j)
            }
            for (j=1;j<=n;j++) {
              ROTATE(v,j,ip,j,iq)
            }
            ++(*nrot);
          }
        }
      }
    }
    for (ip=1;ip<=n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  return(FALSE);                /* if we fall through here, then the 
                                   maximum number of iterations has been
                                   exceeded, and there is no result possible */
}

#undef ROTATE

void eigsrt(double d[], double** v, int n)
{
  int max_index,j,i;
  double max_val;
  
  for (i=1;i<n;i++) {
      
      max_index=i;
      max_val=d[max_index];

      for (j=i+1;j<=n;j++)
          if (d[j] >= max_val) {
              max_index=j;
              max_val=d[max_index];
          }
      
      if (max_index != i) {
          d[max_index]=d[i];
          d[i]=max_val;
          for (j=1;j<=n;j++) {        /* swap columns */
              max_val=v[j][i];
              v[j][i]=v[j][max_index];
              v[j][max_index]=max_val;
          }
      }
  }
  
}


VIO_BOOL eigen(double **inputMat, 
                     int    ndim, 
                     double *eigen_val, 
                     double **eigen_vec, 
                     int    *iters)
{

  double sum,**copy_of_input,**eigvec,*eigval;
  int eig_flag, i,j;

  VIO_ALLOC2D(copy_of_input,ndim+1,ndim+1);
  VIO_ALLOC2D(eigvec,ndim+1,ndim+1);
  ALLOC(eigval,ndim+1);

  for(i=0; i<ndim; i++) {
    for(j=0; j<ndim; j++) {
      copy_of_input[i+1][j+1] = inputMat[i][j];
    }
  }

  eig_flag = jacobi(copy_of_input,ndim,eigval,eigvec,iters);

  if (eig_flag) {
    eigsrt(eigval, eigvec, ndim);
    
    for(i=0; i<ndim; i++)                /* copy to calling parameters */
      for(j=0; j<ndim; j++)
        eigen_vec[i][j] = eigvec[i+1][j+1];
    
    for(i=0; i<ndim; i++)
      eigen_val[i] = eigval[i+1];
    
    for(i=0; i<ndim; i++) {                /* normalize the eigen_Vectors */
      
      sum = 0.0;
      for(j=0; j<ndim; j++)
        sum += eigen_vec[j][i]*eigen_vec[j][i];
      
      sum = sqrt(sum);
      if (sum > 0.0) {
        for(j=0; j<ndim; j++)
          eigen_vec[j][i] /= sum;
      }
      else {
        print ("Can't norm %d (%f) %f %f %f\n",i,sum,eigen_vec[0][i],eigen_vec[1][i],eigen_vec[2][i]);
      }
      
    }
  }
  else {

    for(i=0; i<ndim; i++)                
      for(j=0; j<ndim; j++)
        eigen_vec[i][j] = 0.0;

    for(i=0; i<ndim; i++) {
      eigen_val[i] = 0.0;
      eigen_vec[i][i] = 1.0;
    }

  }
  
  VIO_FREE2D(copy_of_input);
  VIO_FREE2D(eigvec);
  FREE(eigval);

  return(eig_flag);
}




/* ----------------------------------------------------------------------------
@NAME       : eigen
@INPUT      : **inputMat - symetric square input matrix 
              ndim       - dimension of input matrix
@OUTPUT     : *eigen_val - array of sorted (decending order) eigen values
             **eigen_vec - matrix of eigen vectors (stored in columns)
              *iters     - the number of iterations required to find the 
                           result
@RETURNS    : FALSE if possible numerical error, TRUE otherwise
@DESCRIPTION: 

   the goal is to find eigen_val & eigen_vec such that

   inputmat * eigen_vec = eigen_vec * eigen_val

   note: inputMat is modified by the routine.

   the current bit of code for eigen originally came from Pascal @
   Rennes for the analysis of 37 element meg data.  It has been
   modified by LC to replace jacobi() and eigsrt() originally used in
   cov_to_praxes.

@GLOBALS    : none
@CALLS      : none
@CREATED    : sept 95 Louis Collins
@MODIFIED   :
---------------------------------------------------------------------------- */

VIO_BOOL eigen2(double **inputMat, 
                     int    ndim, 
                     double *eigen_val, 
                     double **eigen_vec, 
                     int    *iters)

{

  double
    *vec, *test_vec,
    x, 
    vec_sum, previous_sum, 
    val,     previous_val,
    max_val;
  int 
    possible_error, i,j,k,iterations;

  
  *iters = 0;
  possible_error = FALSE;

  ALLOC(vec, ndim);
  ALLOC(test_vec, ndim);


/*
  print ("[\n");
  for(i=0; i<ndim; i++) {
    for(j=0; j<ndim; j++) {
      print ("%f ",inputMat[i][j] );
    }
    print ("\n");
  }
*/

  /* NOTE: the code here actually treats and stores the eigen vectors
     in row format.  The matrix will be transposed below */

  for(k=0; k<ndim; k++) {

    for(i=0; i<ndim; i++)                /* init vec to ones */
      vec[i]=1;

    iterations=0;
    val=0.0;
    vec_sum=(double)ndim;        /* since  sum += vec[i]  */


    do {                        /* iterate  */

      iterations++;
      previous_val=val;
      previous_sum=vec_sum;
      max_val=0.0;
      vec_sum=0.0;

      for(i=0; i<ndim; i++) {        /* for each row in the inputMat */

        test_vec[i]=0.0;        /* test_vec = inputMat * vec'   */
        for(j=0; j<ndim; j++)          
          test_vec[i] += inputMat[i][j]*vec[j];

        vec_sum += test_vec[i];

        x = fabs(test_vec[i]);        /* find max element in test_vec */
        if(x>max_val) max_val=x;

      }

      val     = vec_sum/previous_sum;
      vec_sum = vec_sum/max_val;

      for(i=0; i<ndim; i++)        /* normalize to largest element*/
        vec[i]=test_vec[i]/max_val;

    } while( (iterations<3000) && (fabs(1.0-previous_val/val)>1e-15) );

print ("%3d (%5d) -> (%15.10f) %15.10f %15.10f %15.10f  \n", k, iterations, val, vec[0], vec[1],vec[2]);

    *iters += iterations;

    if((fabs(val)>1e-6)&&(iterations>2900)) {
      possible_error = TRUE;
    }

    vec_sum=0.0;
    for(i=0; i<ndim; i++) 
      vec_sum += vec[i]*vec[i];

    for(i=0; i<ndim; i++)
      for(j=0; j<ndim; j++)
        inputMat[i][j] -= (vec[i]*vec[j]/vec_sum)*val;

    vec_sum = sqrt(vec_sum);
    for(i=0; i<ndim; i++)
      eigen_vec[k][i]=vec[i]/vec_sum; /* store eig vec in row 'k' */

    eigen_val[k]=val;
  }



                                /* sort eig vectors in order of eigval */
  for(i=0; i<ndim; i++) {
    max_val=eigen_val[k=i];
    for(j=i+1; j<ndim; j++)
      if (eigen_val[j] >= max_val) max_val=eigen_val[k=j];

    if (k != i) {
      eigen_val[k]=eigen_val[i]; /* swap eigen values */
      eigen_val[i]=max_val;

      for(j=0; j<ndim; j++) {        /* swap eigen vectors (rows of eigen_vec)*/
        val=            eigen_vec[i][j];
        eigen_vec[i][j]=eigen_vec[k][j];
        eigen_vec[k][j]=val;
      }

    }

  }
                                /* transpose eigen vector matrix, to return
                                   eigen vectors in the columns */

  for(i=0; i<ndim; i++) {
    for(j=i+1; j<ndim; j++) {
      val = eigen_vec[i][j];
      eigen_vec[i][j] = eigen_vec[j][i];
      eigen_vec[j][i] = val;
    }
  }


  FREE(vec);
  FREE(test_vec);

  return(!possible_error);

}






