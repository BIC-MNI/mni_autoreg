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
#@VERSION    : $Id: vox_space.h,v 1.3 2004-02-12 05:54:16 rotor Exp $
#-----------------------------------------------------------------------------
*/

typedef struct {
   Real              start[3];
   VectorR           directions[3];
   General_transform *voxel_to_voxel_space;
} Voxel_space_struct;

Voxel_space_struct* new_voxel_space_struct(void);

void delete_voxel_space_struct( Voxel_space_struct *vox_space);

void build_reorder_matrix_vox2xyz(General_transform *trans, Volume volume);

void build_reorder_matrix_xyz2vox(General_transform *trans, Volume volume);

void get_into_voxel_space(Arg_Data *globals,
                                 Voxel_space_struct *vox,
                                 Volume v1, Volume v2);

void  my_homogenous_transform_point(Transform  *transform,
                                           Real       x,
                                           Real       y,
                                           Real       z,
                                           Real       w,
                                           Real       *x_trans,
                                           Real       *y_trans,
                                           Real       *z_trans );

