/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_rots.h
@DESCRIPTION: File containing declarations of routines for 
              prototypes for make_rots.c
@CREATED    : Mon Jun  7 11:41:02 EST 1993 (LC) 
@MODIFIED   : 
---------------------------------------------------------------------------- */

public void   
  make_resampling_matrix(float **xmat,
			 float trans_x,float trans_y,float trans_z,
			 float center_x,float center_y,float center_z,
			 float rot_x,float rot_y,float rot_z,
			 float scale_x,float scale_y,float scale_z);
public void
  make_inverted_resampling_matrix(float **xmat,
				  float *trans_x,float *trans_y,float *trans_z,
				  float *center_x,float *center_y,float *center_z,
				  float *rot_x,float *rot_y,float *rot_z,
				  float *scale_x,float *scale_y,float *scale_z);

public void build_transformation_matrix(double lt[3][4], 
					double *center,
					double *translations,
					double *scales,
					double *rotations);

public void build__inverse_transformation_matrix(double lt[3][4], 
						 double *center,
						 double *translations,
						 double *scales,
						 double *rotations);

public Boolean extract_parameters_from_matrix(double lt[3][4], 
					   double *center,
					   double *translations,
					   double *scales,
					   double *rotations);
