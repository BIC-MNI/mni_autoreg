/*
#------------------------------ MNI Header ----------------------------------
#@NAME       : vox_space.h
#@DESCRIPTION: struct and function prototypes for world to voxel
               space manipulations with the lattice 
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : Wed Jun 25, 1997, Louis Collins
#@MODIFIED   : not yet!
#@VERSION    : $Id: vox_space.h,v 1.4 2006-11-29 09:09:32 rotor Exp $
#-----------------------------------------------------------------------------
*/

typedef struct {
   VIO_Real              start[3];
   VectorR           directions[3];
   VIO_General_transform *voxel_to_voxel_space;
} Voxel_space_struct;

Voxel_space_struct* new_voxel_space_struct(void);

void delete_voxel_space_struct( Voxel_space_struct *vox_space);

void build_reorder_matrix_vox2xyz(VIO_General_transform *trans, VIO_Volume volume);

void build_reorder_matrix_xyz2vox(VIO_General_transform *trans, VIO_Volume volume);

void get_into_voxel_space(Arg_Data *globals,
                                 Voxel_space_struct *vox,
                                 VIO_Volume v1, VIO_Volume v2);

void  my_homogenous_transform_point(VIO_Transform  *transform,
                                           VIO_Real       x,
                                           VIO_Real       y,
                                           VIO_Real       z,
                                           VIO_Real       w,
                                           VIO_Real       *x_trans,
                                           VIO_Real       *y_trans,
                                           VIO_Real       *z_trans );

