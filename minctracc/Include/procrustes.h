/* ----------------------------- MNI Header -----------------------------------
@NAME       : procrustes.h
@DESCRIPTION: Header file containing declarations of routines for doing 
              procrustes calculations.
@METHOD     : 
@CALLS      : 
@CREATED    : February 7, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */

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

