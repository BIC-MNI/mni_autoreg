/*------------------------------ MNI Header ----------------------------------
@NAME       : sub_lattice.h
@DESCRIPTION: prototypes for Optimize/sub_lattice.c
@CREATED    : Mon Nov  3, 1997 , Louis Collins
@MODIFIED   : not yet!
@VERSION    : $Id: sub_lattice.h,v 1.4 2006-11-29 09:09:32 rotor Exp $
-----------------------------------------------------------------------------*/


void    
build_source_lattice(VIO_Real x, VIO_Real y, VIO_Real z,
                     float PX[], float PY[], float PZ[],
                     VIO_Real width_x, VIO_Real width_y, VIO_Real width_z, 
                     int nx, int ny, int nz,
                     int ndim, int *length);

void 
go_get_samples_in_source(VIO_Volume data, VIO_Volume mask,
                         float x[], float y[], float z[],
                         float samples[],
                         VIO_BOOL masked_samples_in_source[],
                         int len,
                         int inter_type);

float 
go_get_samples_with_offset(VIO_Volume data, VIO_Volume mask,
                           float *x, float *y, float *z,
                           VIO_Real  dx, VIO_Real  dy, VIO_Real dz,
                           int obj_func,
                           int len,  int *sample_target_count, 
                           float sqrt_s1, float *a1, VIO_BOOL *m1,
                           VIO_BOOL use_nearest_neighbour);

void    
build_target_lattice(float px[], float py[], float pz[],
                     float tx[], float ty[], float tz[],
                     int len, int dim);

void    
build_target_lattice_using_super_sampled_def(float px[], float py[], float pz[],
                                             float tx[], float ty[], float tz[],
                                             int len, int dim);


