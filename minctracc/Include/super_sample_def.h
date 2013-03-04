#ifndef MINCTRACC_SUPER_SAMPLE_DEF_H
#define MINCTRACC_SUPER_SAMPLE_DEF_H

/* build the volume structure and allocate the data space to store
   a super-sampled GRID_TRANSFORM.

   *super_sampled must be ALLOCed before this call
   this routine will alloc the volume data space.

   super_step specifies the number of times to super-sample the data.
 */

void 
create_super_sampled_data_volumes(VIO_General_transform *orig_deformation,
                                  VIO_General_transform *super_sampled,
                                  int super_step);


/*
   This procedure will use cubic interpolation to resample the deformation
   field stored in *orig_deformation->displacement_volume onto the voxel
   lattice defined in *super_sampled.

   The procedure essentially steps though all voxels of 
   *super_sampled->displacement_volume, figures out the world coordinate for
   the voxel position, and uses this info to interpolate into *orig_deformation.

   see interpolate_super_sampled_data_by2() (below) for speedy code for the
   special case when the super sampling rate is euqal to 2.
*/
void interpolate_super_sampled_data(VIO_General_transform *orig_deformation,
                                           VIO_General_transform *super_sampled);


/* build the volume structure and allocate the data space to store
   a super-sampled GRID_TRANSFORM.

   *super_sampled must be ALLOCed before this call
   this routine will alloc the volume data space.

   this is a specialized version of create_super_sampled_data_volumes(),
   specifically coded to speed the interpolation of the super_sampled
   data when the super-sampling rate is equal to 2.

   we will have:
   orig samples:   X   X   X   X   X   X   X
   new samples:    x x x x x x x x x x x x x

   new_steps = old_steps/2;
   new_count = old_count*2 - 1;

 */
void create_super_sampled_data_volumes_by2(VIO_General_transform *orig_deformation,
                                                  VIO_General_transform *super_sampled);



/*
   This procedure will use cubic interpolation to resample the deformation
   field stored in *orig_deformation->displacement_volume onto the voxel
   lattice defined in *super_sampled.

   Since the super_sampled data is exactly aligned with the original data,
   (but has twice as many samples in each direction) it is possible to do
   all interpolation in the voxel space - no calls to the voxel-to-world
   transformations are necessary - thus resulting in a serious speed 
   increase.

   Description:

   suppose we have three successive planes to be interpolated...

            X e X   e f e   X e X       X - original sample
            e f e   f c f   e f e       e - 1st level edge interpolation
            X e X   e f e   X e X       f - 2nd level face interpolation
                                        c - 3rd level center interpolation

   Since the super-sampled data is aligned with the original data, the 
   samples at X need only to be copied from the original data volume.

   The 'e' voxels need be interpolated only using cubic interpolation 
   along the direction of their corresponding edge.  (the spline 
   parameters essentially nul other data).  Only the X voxels are
   used to feed the interpolation.

   The 'f' voxels are interpolated in 2D, using data only from the 
   plane containing the correspoinding face.  Only the 'e' voxels
   need to be used to do the interpolation.

   The 'c' voxels are interpolation in 3D, using data only from the
   'f' voxels.
   
*/
void interpolate_super_sampled_data_by2( VIO_General_transform *orig_deformation,
                                                VIO_General_transform *super_sampled);

#endif
