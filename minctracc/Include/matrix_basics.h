/* ----------------------------- MNI Header -----------------------------------
@NAME       : matrix_basics.h
@DESCRIPTION: File containing declarations of routines for doing basic 
              matrix calculations
@CREATED    : January 31, 1992 (Peter Neelin)
@MODIFIED   : Mon Jun  7 11:41:02 EST 1993 (LC) update with new routines
---------------------------------------------------------------------------- */

public void printmatrix(int rows, int cols, float **the_matrix);

public void calc_centroid(int npoints, int ndim, float **points, 
                          float *centroid);

public void translate(int npoints, int ndim, float **points, 
                      float *translation, float **newpoints);

public void transpose(int rows, int cols, float **mat, float **mat_transpose);

public void invertmatrix(int n, float **mat, float **mat_invert);

public void raw_matrix_multiply(int ldim, int mdim, int ndim, 
                                float **Amat, float **Bmat, float **Cmat);

public void matrix_multiply(int ldim, int mdim, int ndim, 
                            float **Amat, float **Bmat, float **Cmat);

public float trace(int size, float **the_matrix);

public void matrix_scalar_multiply(int rows, int cols, float scalar, 
                            float **the_matrix, float **product);

public void nr_identd(double **A, int m1, int m2, int n1, int n2 );
public void nr_identf(float **A, int m1, int m2, int n1, int n2 );

public void nr_copyd(double **A, int m1, int m2, int n1, int n2, double **B );
public void nr_copyf(float  **A, int m1, int m2, int n1, int n2, float **B );

public void nr_rotxd(double **M, double a);
public void nr_rotxf(float **M, float a);

public void nr_rotyd(double **M,double a);
public void nr_rotyf(float **M, float a);

public void nr_rotzd(double **M,double a);
public void nr_rotzf(float **M, float a);

public void nr_multd(double **A, int mA1, int mA2, int nA1, int nA2,
		     double **B, int mB1, int mB2, int nB1, int nB2, 
		     double **C);

public void nr_multf(float **A, int mA1, int mA2, int nA1, int nA2,
		     float **B, int mB1, int mB2, int nB1, int nB2, 
		     float **C);


public void transformations_to_homogeneous(int ndim, 
                  float *translation, float *centre_of_rotation,
                  float **rotation, float scale,
                  float **transformation);

public void translation_to_homogeneous(int ndim, float *translation,
                                       float **transformation);

public void rotation_to_homogeneous(int ndim, float **rotation,
                                       float **transformation);

