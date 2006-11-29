/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_rots.h
@DESCRIPTION: File containing declarations of routines for 
              prototypes for make_rots.c
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

@CREATED    : Mon Jun  7 11:41:02 EST 1993 (LC) 
@MODIFIED   : 
---------------------------------------------------------------------------- */


void   
  make_resampling_matrix(float **xmat,
                         float trans_x,float trans_y,float trans_z,
                         float center_x,float center_y,float center_z,
                         float rot_x,float rot_y,float rot_z,
                         float scale_x,float scale_y,float scale_z);
void
  make_inverted_resampling_matrix(float **xmat,
                                  float *trans_x,float *trans_y,float *trans_z,
                                  float *center_x,float *center_y,float *center_z,
                                  float *rot_x,float *rot_y,float *rot_z,
                                  float *scale_x,float *scale_y,float *scale_z);



void build_transformation_matrix_quater(VIO_Transform *trans,
                                               double *center,
                                               double *translations,
                                               double *scales,
                                               double *shears,
                                               double *quaternions);


void build_inverse_transformation_matrix_quater(VIO_Transform *trans,
                                                       double *center,
                                                       double *translations,
                                                       double *scales,
                                                       double *shears,
                                                       double *quaternions);


VIO_BOOL extract2_parameters_from_matrix_quater(VIO_Transform *trans,
                                                      double *center,
                                                      double *translations,
                                                      double *scales,
                                                      double *shears,
                                                      double *quaternions);
