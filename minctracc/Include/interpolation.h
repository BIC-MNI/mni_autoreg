#ifndef MINCTRACC_INTERPOLATION_H
#define MINCTRACC_INTERPOLATION_H

#include "point_vector.h"


int trilinear_interpolant(Volume volume, 
                                 PointR *coord, double *result);

int tricubic_interpolant(Volume volume, 
                                PointR *coord, double *result);

void do_Ncubic_interpolation(Volume volume, 
                                    long index[], int cur_dim, 
                                    double frac[], double *result);

int nearest_neighbour_interpolant(Volume volume, 
                                         PointR *coord, double *result);

/* A point is not masked if it is a point we should consider.
   If the mask volume is NULL, we consider all points.
   Otherwise, consider a point if the mask volume value is > 0.
*/
int point_not_masked(Volume volume, 
			    Real wx, Real wy, Real wz);


int voxel_point_not_masked(Volume volume, 
                                  Real vx, Real vy, Real vz);


#endif
