/* ----------------------------- MNI Header -----------------------------------
@NAME       : procrustes.h
@DESCRIPTION: Header file containing declarations of routines for doing 
              procrustes calculations.
@METHOD     : 
@CALLS      : 
@COPYRIGHT  :
              Copyright 1993 Peter Neelin, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.

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

