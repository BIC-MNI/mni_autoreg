/*------------------------------ MNI Header ----------------------------------
@NAME       : sub_lattice.h
@DESCRIPTION: prototypes for Optimize/sub_lattice.c
@CREATED    : Mon Nov  3, 1997 , Louis Collins
@MODIFIED   : not yet!
@VERSION    : $Id: sub_lattice.h,v 1.2 2004-02-12 05:54:16 rotor Exp $
-----------------------------------------------------------------------------*/


void    
build_source_lattice(Real x, Real y, Real z,
                     float PX[], float PY[], float PZ[],
                     Real width_x, Real width_y, Real width_z, 
                     int nx, int ny, int nz,
                     int ndim, int *length);

void 
go_get_samples_in_source(Volume data,
                         float x[], float y[], float z[],
                         float samples[],
                         int len,
                         int inter_type);

float 
go_get_samples_with_offset(Volume data,
                           float *x, float *y, float *z,
                           Real  dx, Real  dy, Real dz,
                           int obj_func,
                           int len, 
                           float sqrt_s1, float *a1,
                           BOOLEAN use_nearest_neighbour);

void    
build_target_lattice(float px[], float py[], float pz[],
                     float tx[], float ty[], float tz[],
                     int len, int dim);

void    
build_target_lattice_using_super_sampled_def(float px[], float py[], float pz[],
                                             float tx[], float ty[], float tz[],
                                             int len, int dim);
