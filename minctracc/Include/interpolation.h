#ifndef MINCTRACC_INTERPOLATION_H
#define MINCTRACC_INTERPOLATION_H

#include "point_vector.h"


public int trilinear_interpolant(Volume volume, 
                                 PointR *coord, double *result);

public int tricubic_interpolant(Volume volume, 
                                PointR *coord, double *result);

public void do_Ncubic_interpolation(Volume volume, 
                                    long index[], int cur_dim, 
                                    double frac[], double *result);

public int nearest_neighbour_interpolant(Volume volume, 
                                         PointR *coord, double *result);

/* A point is not masked if it is a point we should consider.
   If the mask volume is NULL, we consider all points.
   Otherwise, consider a point if the mask volume value is > 0.
*/
public int point_not_masked(Volume volume, 
			    Real wx, Real wy, Real wz);


public int voxel_point_not_masked(Volume volume, 
                                  Real vx, Real vy, Real vz);


#endif
