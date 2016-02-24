/* ----------------------------- MNI Header -----------------------------------
@NAME       : super_sample_def.c
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: procedures to super sample the deformation field.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@COPYRIGHT  :
              Copyright 1995 Louis Collins, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.

@CREATED    : 
@MODIFIED   : $Log: super_sample_def.c,v $
@MODIFIED   : Revision 96.13  2006-11-30 09:07:33  rotor
@MODIFIED   :  * many more changes for clean minc 2.0 build
@MODIFIED   :
@MODIFIED   : Revision 96.12  2006/11/29 09:09:34  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 96.11  2005/07/20 20:45:51  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 96.10  2004/02/12 06:08:21  rotor
@MODIFIED   :  * removed /static defs
@MODIFIED   :
@MODIFIED   : Revision 96.9  2004/02/04 20:44:13  lenezet
@MODIFIED   : *** empty log message ***
@MODIFIED   :
@MODIFIED   : Revision 96.8  2003/02/26 00:56:38  lenezet
@MODIFIED   : for 2D : now computes all 3 coordinates for the "start" (to take into account the slice position).
@MODIFIED   : simplification of build_lattices.
@MODIFIED   : bug correction in amoeba_NL_obj_function.
@MODIFIED   :
@MODIFIED   : Revision 96.7  2002/12/13 21:18:20  lenezet
@MODIFIED   :
@MODIFIED   : A memory leak has been repaired
@MODIFIED   :
@MODIFIED   : Revision 96.6  2002/11/20 21:39:16  lenezet
@MODIFIED   :
@MODIFIED   : Fix the code to take in consideration the direction cosines especially in the grid transform.
@MODIFIED   : Add an option to choose the maximum expected deformation magnitude.
@MODIFIED   :
@MODIFIED   : Revision 96.5  2002/03/26 14:15:46  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 96.4  2002/03/07 19:08:56  louis
@MODIFIED   : Added -lattice_diameter as an optionto minctracc to account for a
@MODIFIED   : problem with the automated calculation of the sub-lattice diameter.
@MODIFIED   : It used to be step*3*2 - which was pretty big, when step = 8mm.
@MODIFIED   :
@MODIFIED   : Now, the sub lattice diameter can be input on the command line, and I
@MODIFIED   : suggest a lattice size 3 times greater than the step size.
@MODIFIED   :
@MODIFIED   : If not on the command line, the default is = 24mm.
@MODIFIED   :
@MODIFIED   : Revision 96.3  2000/04/05 04:38:25  stever
@MODIFIED   : * Reordered code in super_sample_def.c so that it compiles with
@MODIFIED   :   GCC under IRIX 5.3.
@MODIFIED   : * set VERSION to 0.98i.
@MODIFIED   :
@MODIFIED   : Revision 96.2  1997/11/12 21:07:43  louis
@MODIFIED   : no changes, other than rcsid...
@MODIFIED   :
 * Revision 96.1  1997/11/03  15:06:29  louis
 * working version, before creation of mni_animal package, and before inserting
 * distance transforms
 *
 * Revision 96.1  1997/11/03  15:06:29  louis
 * working version, before creation of mni_animal package, and before inserting
 * distance transforms
 *
 * Revision 96.0  1996/08/21  18:22:10  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:22:06  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:16:03  louis
 * Never released version 0.95
 *
 * Revision 1.2  1996/08/12  14:15:58  louis
 * Pre-release
 *
 * Revision 1.1  1996/04/01  09:16:29  collins
 * Initial revision
 *
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/minctracc/Optimize/super_sample_def.c,v 96.13 2006-11-30 09:07:33 rotor Exp $";
#endif

#include <config.h>
#include <volume_io.h>
#include <Proglib.h>
#include "constants.h"
#include "minctracc_point_vector.h"

                                /* prototypes called: */

void get_volume_XYZV_indices(VIO_Volume data, int xyzv[]);
void init_the_volume_to_zero(VIO_Volume volume);
void interpolate_deformation_slice(VIO_Volume volume, 
                                          VIO_Real wx, VIO_Real wy, VIO_Real wz,
                                          VIO_Real def[]);
void set_up_lattice(VIO_Volume data, 
                           double *user_step, 
                           double *start,     
                           double *wstart,     
                           int    *count,     
                           double *step,      
                           VectorR directions[]);

static void interpolate_super_sampled_data_by2_dim2(VIO_General_transform *orig_deformation,
                                                     VIO_General_transform *super_sampled,
                                                     int i);
static void interpolate_super_sampled_data_by2_dim2_X(VIO_General_transform *orig_deformation,
                                                       VIO_General_transform *super_sampled);
static void interpolate_super_sampled_data_by2_dim2_Y(VIO_General_transform *orig_deformation,
                                                       VIO_General_transform *super_sampled);
static void interpolate_super_sampled_data_by2_dim2_Z(VIO_General_transform *orig_deformation,
                                                       VIO_General_transform *super_sampled);
/* build the volume structure and allocate the data space to store
   a super-sampled GRID_TRANSFORM.

   *super_sampled must be ALLOCed before this call
   this routine will alloc the volume data space.

   super_step specifies the number of times to super-sample the data.
 */
void create_super_sampled_data_volumes(VIO_General_transform *orig_deformation,
                                              VIO_General_transform *super_sampled,
                                              int super_step)

{

 int
    i,
    xyzv[VIO_MAX_DIMENSIONS],
    orig_count[VIO_MAX_DIMENSIONS], 
    xyz_count[VIO_MAX_DIMENSIONS], 
    new_count[VIO_MAX_DIMENSIONS];
  VIO_Real
    voxel[VIO_MAX_DIMENSIONS],
    start[VIO_MAX_DIMENSIONS],
    wstart[VIO_MAX_DIMENSIONS],
    orig_steps[VIO_MAX_DIMENSIONS],
    new_steps[VIO_MAX_DIMENSIONS],
    xyz_steps[VIO_MAX_DIMENSIONS];

  VectorR 
    directions[3];



  if (orig_deformation->type != GRID_TRANSFORM) {
    print_error_and_line_num("create_super_sampled_data_volumes not called with GRID_TRANSFORM",
                             __FILE__, __LINE__);
  }

                                /* copy the transform definition */
  *super_sampled = *orig_deformation; 
	super_sampled->displacement_volume_file=NULL;

                                /* copy the GRID_TRANSFORM definition */
  super_sampled->displacement_volume = 
    copy_volume_definition_no_alloc(orig_deformation->displacement_volume,
                                    NC_UNSPECIFIED, FALSE, 0.0, 0.0);

                                /* prepare to modify the GRID_TRANSFORM */

  get_volume_XYZV_indices(orig_deformation->displacement_volume, xyzv);
  get_volume_sizes(       orig_deformation->displacement_volume, orig_count);
  get_volume_separations( orig_deformation->displacement_volume, orig_steps);
  
  for(i=0; i<3; i++) {                /* set_up_lattice needs XYZ order.... */
    new_steps[i] = orig_steps[xyzv[i]]/super_step;
  }
  
  set_up_lattice(orig_deformation->displacement_volume, 
                 new_steps, start, wstart, xyz_count, xyz_steps, directions);
      /* use the sign of the step returned to set the true step size */
  for(i=0; i<3; i++)        
    if (xyz_steps[i]<0) 
      xyz_steps[i] = -1.0 * fabs(new_steps [i]); 
    else 
      xyz_steps[i] = fabs(new_steps[i]);
  
  for(i=0; i<VIO_MAX_DIMENSIONS; i++) {
    voxel[i] = 0.0;
    new_count[i] = orig_count[i];
    new_steps[i] = orig_steps[i];
  }
  for(i=0; i<3; i++) {
    new_count[ xyzv[i] ] = xyz_count[i];
    new_steps[ xyzv[i] ] = xyz_steps[i];
  }

  


  set_volume_sizes(       super_sampled->displacement_volume, new_count);
  set_volume_separations( super_sampled->displacement_volume, new_steps);
  set_volume_starts( super_sampled->displacement_volume,  start);
  
  alloc_volume_data(      super_sampled->displacement_volume );


}
 
/*
   This procedure will use cubic interpolation to resample the deformation
   field stored in *orig_deformation->displacement_volume onto the voxel
   lattice defined in *super_sampled.

   The procedure essentially steps though all voxels of 
   *super_sampled->displacement_volume, figures out the world coordinate for
   the voxel position, and uses this info to interpolate into *orig_deformation.

   see interpolate_super_sampled_data_by2() (below) for speedy code for the
   special case when the super sampling rate is euqal to 2.



         W A R N I N G : there appears to be a bug in volume_io's evaluate_volume
                         call, so this routine should not be called.


*/
void interpolate_super_sampled_data(VIO_General_transform *orig_deformation,
                                           VIO_General_transform *super_sampled)
{
  VIO_Volume
    orig_vol,
    super_vol;

  int 
    i,
    index[VIO_MAX_DIMENSIONS],
    orig_xyzv[VIO_MAX_DIMENSIONS],
    xyzv[VIO_MAX_DIMENSIONS],
    count[VIO_MAX_DIMENSIONS];
  VIO_Real
    def_vector[VIO_N_DIMENSIONS],
    voxel[VIO_MAX_DIMENSIONS],
    wx,wy,wz;
  VIO_progress_struct
    progress;

  print ("W A R N I N G : there appears to be a bug in volume_io's evaluate_volume\ncall, so interpolate_super_sampled_data should not be called.\n");




  orig_vol = orig_deformation->displacement_volume;
  super_vol= super_sampled->displacement_volume;

  get_volume_sizes(super_vol, count);
  get_volume_XYZV_indices(super_vol, xyzv);
  get_volume_XYZV_indices(orig_vol, orig_xyzv);
  
  initialize_progress_report( &progress, FALSE, count[ xyzv[VIO_X] ],
                             "Super-sampling defs:" );
  
  for(i=0; i<VIO_MAX_DIMENSIONS; i++) index[i]=0;
  

  /* check here to see if special case of super sampling by 2, and if so, call the correct routine dec 02 */

  for(index[xyzv[VIO_X]]=0; index[xyzv[VIO_X]]<count[xyzv[VIO_X]]; index[xyzv[VIO_X]]++) {
    for(index[xyzv[VIO_Y]]=0; index[xyzv[VIO_Y]]<count[xyzv[VIO_Y]]; index[xyzv[VIO_Y]]++) {
      for(index[xyzv[VIO_Z]]=0; index[xyzv[VIO_Z]]<count[xyzv[VIO_Z]]; index[xyzv[VIO_Z]]++) {
        
        for(i=0; i<VIO_MAX_DIMENSIONS; i++) 
          voxel[i]=(VIO_Real)index[i];
        convert_voxel_to_world(super_vol, 
                               voxel,
                               &wx, &wy, &wz);

        evaluate_volume_in_world(orig_vol, wx,wy,wz, 2, TRUE, 0.0,
                                 def_vector,
                                 NULL, NULL, NULL,
                                 NULL, NULL, NULL, NULL, NULL, NULL);

        for(index[xyzv[VIO_Z+1]]=0; index[xyzv[VIO_Z+1]]<count[xyzv[VIO_Z+1]]; index[xyzv[VIO_Z+1]]++) 
          set_volume_real_value(super_vol,
                                index[0],index[1],index[2],index[3],index[4],
                                def_vector[ index[ xyzv[VIO_Z+1] ] ]);
        
      }
    }
    update_progress_report( &progress, index[ xyzv[VIO_X] ]+1 );
  }
  terminate_progress_report( &progress );


}

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
                                                  VIO_General_transform *super_sampled)

{

  int
    i,
    xyzv[VIO_MAX_DIMENSIONS],
    new_xyzv[VIO_MAX_DIMENSIONS], 
    orig_count[VIO_MAX_DIMENSIONS], 
    new_count[VIO_MAX_DIMENSIONS];
  VIO_Real
    dirs[VIO_MAX_DIMENSIONS],
    voxel[VIO_MAX_DIMENSIONS],
    wstart[VIO_MAX_DIMENSIONS],
    start[VIO_MAX_DIMENSIONS],
    orig_steps[VIO_MAX_DIMENSIONS],
    new_steps[VIO_MAX_DIMENSIONS];

 

  if (orig_deformation->type != GRID_TRANSFORM) {
    print_error_and_line_num("create_super_sampled_data_volumes not called with GRID_TRANSFORM",
                             __FILE__, __LINE__);
  }

                                /* copy the transform definition */
  *super_sampled = *orig_deformation; 

                                /* copy the GRID_TRANSFORM definition */
  super_sampled->displacement_volume = 
    copy_volume_definition_no_alloc(orig_deformation->displacement_volume,
                                    NC_UNSPECIFIED, FALSE, 0.0, 0.0);

  super_sampled->displacement_volume_file = NULL;

                                /* prepare to modify the GRID_TRANSFORM */

  get_volume_XYZV_indices(orig_deformation->displacement_volume, xyzv);
  get_volume_XYZV_indices(super_sampled->displacement_volume,    new_xyzv);
  get_volume_sizes(       orig_deformation->displacement_volume, orig_count);
  get_volume_separations( orig_deformation->displacement_volume, orig_steps);

  for(i=0; i<get_volume_n_dimensions(orig_deformation->displacement_volume); i++) {                
    new_steps[new_xyzv[i]] = orig_steps[xyzv[i]];
    new_count[new_xyzv[i]] = orig_count[xyzv[i]];
  }
                                /* set up the new step size */
  for(i=0; i<3; i++) {                
    new_steps[new_xyzv[i]] = orig_steps[xyzv[i]]/2.0;
  }
   
                                /* set up the new number of elements size */
  for(i=0; i<3; i++) {                 
    if (orig_count[xyzv[i]] > 1) { 
      new_count[new_xyzv[i]] = orig_count[xyzv[i]]*2.0 - 1;
    }
    else {
      new_count[new_xyzv[i]] = orig_count[xyzv[i]];
    }
  }

  for(i=0; i<VIO_MAX_DIMENSIONS; i++) 
    voxel[i] = 0.0;
  convert_voxel_to_world(orig_deformation->displacement_volume,
                         voxel,
                         &start[VIO_X], &start[VIO_Y], &start[VIO_Z]);
 
  
                                /* write the new info to the volume struct
                                   and allocate the space for the volume data */
  set_volume_sizes(       super_sampled->displacement_volume, new_count);
  set_volume_separations( super_sampled->displacement_volume, new_steps);
  set_volume_translation( super_sampled->displacement_volume, voxel, start);
  alloc_volume_data(      super_sampled->displacement_volume );

}

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
#define MY_CUBIC_05(a1,a2,a3,a4)  \
   ( ( -(a1) + 9*(a2) + 9*(a3) -(a4) ) / 16.0 )

static void interpolate_super_sampled_data_by2_dim3( 
    VIO_General_transform *orig_deformation,
    VIO_General_transform *super_sampled )
{
  VIO_Volume
    orig_vol,
    super_vol;

  int 
    i, save_index, save_index1, save_index2,
    count,
    index[VIO_MAX_DIMENSIONS],
    sindex[VIO_MAX_DIMENSIONS],
    orig_xyzv[VIO_MAX_DIMENSIONS],
    xyzv[VIO_MAX_DIMENSIONS],
    orig_count[VIO_MAX_DIMENSIONS],
    super_count[VIO_MAX_DIMENSIONS];
  VIO_Real
    value1, value2, value3, v1,v2,v3,v4,v5,v6;
  VIO_progress_struct
    progress;



    orig_vol = orig_deformation->displacement_volume;
    super_vol= super_sampled->displacement_volume;
    
    get_volume_sizes(       super_vol, super_count);
    get_volume_sizes(       orig_vol,  orig_count);
    get_volume_XYZV_indices(super_vol, xyzv);
    get_volume_XYZV_indices(orig_vol,  orig_xyzv);
    
    init_the_volume_to_zero(super_sampled->displacement_volume);

    initialize_progress_report(&progress, FALSE, 
                               super_count[ xyzv[VIO_X] ]*super_count[ xyzv[VIO_Y] ]*
                               super_count[ xyzv[VIO_Z] ]*super_count[ xyzv[VIO_Z+1] ]+1,
                               "Super-sampling defs:" );
    count = 0;



    /* LEVEL 0: copy original 'corner' nodes, identified as 'X' in desc above */

    for(i=0; i<VIO_MAX_DIMENSIONS; i++) sindex[i]=index[i]=0;
    
    for(index[orig_xyzv[VIO_X]]=0; index[orig_xyzv[VIO_X]]<orig_count[orig_xyzv[VIO_X]]; index[orig_xyzv[VIO_X]]++) {
      for(index[orig_xyzv[VIO_Y]]=0; index[orig_xyzv[VIO_Y]]<orig_count[orig_xyzv[VIO_Y]]; index[orig_xyzv[VIO_Y]]++) {
        for(index[orig_xyzv[VIO_Z]]=0; index[orig_xyzv[VIO_Z]]<orig_count[orig_xyzv[VIO_Z]]; index[orig_xyzv[VIO_Z]]++) {

          for(i=0; i<VIO_N_DIMENSIONS; i++)
            sindex[ xyzv[i] ] = 2*index[ orig_xyzv[i] ];

          for(index[orig_xyzv[VIO_Z+1]]=0; index[orig_xyzv[VIO_Z+1]]<orig_count[xyzv[VIO_Z+1]]; index[orig_xyzv[VIO_Z+1]]++) {
            sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];
            GET_VOXEL_4D(value1, orig_vol, index[0],index[1],index[2],index[3]);
            SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          }          

          count += 3;
        }
      }
      update_progress_report( &progress, count+1 );
    }


    /* LEVEL 1: edge interpolation, identified as 'e' in desc above */

               /* do edges along the index[ orig_xyzv[VIO_X] ] dir */

    for(i=0; i<VIO_MAX_DIMENSIONS; i++) sindex[i]=index[i]=0;
    
                                /* loop over all x-dirs  */

    for(index[orig_xyzv[VIO_Y]]=0; index[orig_xyzv[VIO_Y]]<orig_count[orig_xyzv[VIO_Y]]; index[orig_xyzv[VIO_Y]]++) {


      for(index[orig_xyzv[VIO_Z]]=0; index[orig_xyzv[VIO_Z]]<orig_count[orig_xyzv[VIO_Z]]; index[orig_xyzv[VIO_Z]]++) {

        sindex[ xyzv[VIO_Y] ] = 2*index[ orig_xyzv[VIO_Y] ];
        sindex[ xyzv[VIO_Z] ] = 2*index[ orig_xyzv[VIO_Z] ];

        for(index[orig_xyzv[VIO_Z+1]]=0; index[orig_xyzv[VIO_Z+1]]<orig_count[xyzv[VIO_Z+1]]; index[orig_xyzv[VIO_Z+1]]++) {

                                /* do linear interp at ends */

          sindex[ xyzv[VIO_X] ] = 1;                             /* beginning end */
          index[ orig_xyzv[VIO_X] ]=0;             
          GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
          index[ orig_xyzv[VIO_X] ]=1;
          GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
          value1 = (v1+v2)/2;

          sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];

          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;

          sindex[ xyzv[VIO_X] ] = 2*orig_count[ orig_xyzv[VIO_X] ] -3;      /* ending end */
          index[ orig_xyzv[VIO_X] ]= orig_count[ orig_xyzv[VIO_X] ]-2;              
          GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
          index[ orig_xyzv[VIO_X] ]= orig_count[ orig_xyzv[VIO_X] ]-1;
          GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
          value1 = (v1+v2)/2;

          sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];

          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
                                /* now do the voxels between the two ends */


          for(index[orig_xyzv[VIO_X]]=1; index[orig_xyzv[VIO_X]]<=orig_count[orig_xyzv[VIO_X]]-3; index[orig_xyzv[VIO_X]]++) {

            save_index = index[ orig_xyzv[VIO_X] ];
            sindex[ xyzv[VIO_X] ] = index[ orig_xyzv[VIO_X] ]*2 + 1;

            index[ orig_xyzv[VIO_X] ]--;             
            GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_X] ]++;             
            GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_X] ]++;             
            GET_VOXEL_4D(v3, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_X] ]++;             
            GET_VOXEL_4D(v4, orig_vol, index[0],index[1],index[2],index[3]);

            value1 = MY_CUBIC_05(v1,v2,v3,v4);

            SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
            count++;

            index[ orig_xyzv[VIO_X] ] = save_index;
          }

        }          

      }
      update_progress_report( &progress, count+1 );
    }
               /* do edges along the index[ orig_xyzv[VIO_Y] ] dir */

    for(i=0; i<VIO_MAX_DIMENSIONS; i++) sindex[i]=index[i]=0;
    
                                /* loop over all x-dirs  */

    for(index[orig_xyzv[VIO_Z]]=0; index[orig_xyzv[VIO_Z]]<orig_count[orig_xyzv[VIO_Z]]; index[orig_xyzv[VIO_Z]]++) {


      for(index[orig_xyzv[VIO_X]]=0; index[orig_xyzv[VIO_X]]<orig_count[orig_xyzv[VIO_X]]; index[orig_xyzv[VIO_X]]++) {


        sindex[ xyzv[VIO_Z] ] = 2*index[ orig_xyzv[VIO_Z] ];
        sindex[ xyzv[VIO_X] ] = 2*index[ orig_xyzv[VIO_X] ];

        for(index[orig_xyzv[VIO_Z+1]]=0; index[orig_xyzv[VIO_Z+1]]<orig_count[xyzv[VIO_Z+1]]; index[orig_xyzv[VIO_Z+1]]++) {

                                /* do linear interp at ends */

          sindex[ xyzv[VIO_Y] ] = 1;                             /* beginning end */
          index[ orig_xyzv[VIO_Y] ]=0;             
          GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
          index[ orig_xyzv[VIO_Y] ]=1;
          GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
          value1 = (v1+v2)/2;

          sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];

          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;

          sindex[ xyzv[VIO_Y] ] = 2*orig_count[ xyzv[VIO_Y] ] -3;      /* ending end */
          index[ orig_xyzv[VIO_Y] ]= orig_count[ orig_xyzv[VIO_Y] ]-2;              
          GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
          index[ orig_xyzv[VIO_Y] ]= orig_count[ orig_xyzv[VIO_Y] ]-1;
          GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
          value1 = (v1+v2)/2;

          sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];

          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;

                                /* now do the voxels between the two ends */


          for(index[orig_xyzv[VIO_Y]]=1; index[orig_xyzv[VIO_Y]]<=orig_count[orig_xyzv[VIO_Y]]-3; index[orig_xyzv[VIO_Y]]++) {

            save_index = index[ orig_xyzv[VIO_Y] ];
            sindex[ xyzv[VIO_Y] ] = index[ orig_xyzv[VIO_Y] ]*2 + 1;

            index[ orig_xyzv[VIO_Y] ]--;             
            GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_Y] ]++;             
            GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
             index[ orig_xyzv[VIO_Y] ]++;             
            GET_VOXEL_4D(v3, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_Y] ]++;             
            GET_VOXEL_4D(v4, orig_vol, index[0],index[1],index[2],index[3]);

            value1 = MY_CUBIC_05(v1,v2,v3,v4);

            SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
            count++;

            index[ orig_xyzv[VIO_Y] ] = save_index;
          }

        }          

      }
      update_progress_report( &progress, count+1 );
    }

               /* do edges along the index[ orig_xyzv[VIO_Z] ] dir */

    for(i=0; i<VIO_MAX_DIMENSIONS; i++) sindex[i]=index[i]=0;
    
                                /* loop over all Z-dirs  */

    for(index[orig_xyzv[VIO_X]]=0; index[orig_xyzv[VIO_X]]<orig_count[orig_xyzv[VIO_X]]; index[orig_xyzv[VIO_X]]++) {


      for(index[orig_xyzv[VIO_Y]]=0; index[orig_xyzv[VIO_Y]]<orig_count[orig_xyzv[VIO_Y]]; index[orig_xyzv[VIO_Y]]++) {

        sindex[ xyzv[VIO_X] ] = 2*index[ orig_xyzv[VIO_X] ];
        sindex[ xyzv[VIO_Y] ] = 2*index[ orig_xyzv[VIO_Y] ];

        for(index[orig_xyzv[VIO_Z+1]]=0; index[orig_xyzv[VIO_Z+1]]<orig_count[xyzv[VIO_Z+1]]; index[orig_xyzv[VIO_Z+1]]++) {

                                /* do linear interp at ends */

          sindex[ xyzv[VIO_Z] ] = 1;                             /* beginning end */
          index[ orig_xyzv[VIO_Z] ]=0;             
          GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
          index[ orig_xyzv[VIO_Z] ]=1;
          GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
          value1 = (v1+v2)/2;

          sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];

          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;

          sindex[ xyzv[VIO_Z] ] = 2*orig_count[ xyzv[VIO_Z] ] -3;      /* ending end */
          index[ orig_xyzv[VIO_Z] ]= orig_count[ orig_xyzv[VIO_Z] ]-2;              
          GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
          index[ orig_xyzv[VIO_Z] ]= orig_count[ orig_xyzv[VIO_Z] ]-1;
          GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
          value1 = (v1+v2)/2;

          sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];

          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;

                                /* now do the voxels between the two ends */


          for(index[orig_xyzv[VIO_Z]]=1; index[orig_xyzv[VIO_Z]]<=orig_count[orig_xyzv[VIO_Z]]-3; index[orig_xyzv[VIO_Z]]++) {

            save_index = index[ orig_xyzv[VIO_Z] ];
            sindex[ xyzv[VIO_Z] ] = index[ orig_xyzv[VIO_Z] ]*2 + 1;

            index[ orig_xyzv[VIO_Z] ]--;             
            GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_Z] ]++;             
            GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_Z] ]++;             
            GET_VOXEL_4D(v3, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_Z] ]++;             
            GET_VOXEL_4D(v4, orig_vol, index[0],index[1],index[2],index[3]);

            value1 = MY_CUBIC_05(v1,v2,v3,v4);

            SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
            count++;

            index[ orig_xyzv[VIO_Z] ] = save_index;
          }

        }          

      }
      update_progress_report( &progress, count+1 );
    }


    /* LEVEL 2: face interpolation, identified as 'f' in desc above */

               /* do faces in the plane with index[ orig_xyzv[VIO_X] ]=CONST */

    for(i=0; i<VIO_MAX_DIMENSIONS; i++)  sindex[i]=index[i]=0;
    
                                /* loop over all X planes  */

    for(index[orig_xyzv[VIO_X]]=0; index[orig_xyzv[VIO_X]]<orig_count[orig_xyzv[VIO_X]]-1; index[orig_xyzv[VIO_X]]++) {

      sindex[ xyzv[VIO_X] ] = 2*index[ orig_xyzv[VIO_X] ];

      for(index[orig_xyzv[VIO_Z+1]]=0; index[orig_xyzv[VIO_Z+1]]<orig_count[xyzv[VIO_Z+1]]; index[orig_xyzv[VIO_Z+1]]++) {

        sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];


                                /* do faces near the edge first */

        for(index[orig_xyzv[VIO_Z]]=0; index[orig_xyzv[VIO_Z]]<orig_count[orig_xyzv[VIO_Z]]-1; index[orig_xyzv[VIO_Z]]++) {

          sindex[ xyzv[VIO_Z] ] = index[ orig_xyzv[VIO_Z] ]*2 + 1;
          
          sindex[ xyzv[VIO_Y] ] = 0;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] = 2;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          
          value1 = (v1+v2)/2;
          
          sindex[ xyzv[VIO_Y] ] = 1;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
          
          sindex[ xyzv[VIO_Y] ] = 2*orig_count[ orig_xyzv[VIO_Y] ] -3;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] = 2*orig_count[ orig_xyzv[VIO_Y] ] -1;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          
          value1 = (v1+v2)/2;
          
          sindex[ xyzv[VIO_Y] ] = 2*orig_count[ orig_xyzv[VIO_Y] ] -2;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
          
          
        }
        update_progress_report( &progress, count+1 );

        for(index[orig_xyzv[VIO_Y]]=0; index[orig_xyzv[VIO_Y]]<orig_count[orig_xyzv[VIO_Y]]-1; index[orig_xyzv[VIO_Y]]++) {

          sindex[ xyzv[VIO_Y] ] = index[ orig_xyzv[VIO_Y] ]*2 + 1;
          
          sindex[ xyzv[VIO_Z] ] = 0;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Z] ] = 2;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          
          value1 = (v1+v2)/2;
          
          sindex[ xyzv[VIO_Z] ] = 1;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
          
          sindex[ xyzv[VIO_Z] ] = 2*orig_count[ orig_xyzv[VIO_Z] ] -3;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Z] ] = 2*orig_count[ orig_xyzv[VIO_Z] ] -1;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          
          value1 = (v1+v2)/2;
          
          sindex[ xyzv[VIO_Z] ] = 2*orig_count[ orig_xyzv[VIO_Z] ] -2;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
          
          
        }
        update_progress_report( &progress, count+1 );
      
                                /* now do faces in the middle */


        for(index[orig_xyzv[VIO_Y]]=1; index[orig_xyzv[VIO_Y]]<orig_count[orig_xyzv[VIO_Y]]-2; index[orig_xyzv[VIO_Y]]++) {
          for(index[orig_xyzv[VIO_Z]]=1; index[orig_xyzv[VIO_Z]]<orig_count[orig_xyzv[VIO_Z]]-2; index[orig_xyzv[VIO_Z]]++) {

            sindex[ xyzv[VIO_Y] ] = index[ orig_xyzv[VIO_Y] ]*2 + 1;
            sindex[ xyzv[VIO_Z] ] = index[ orig_xyzv[VIO_Z] ]*2 + 1;

            save_index1 =  sindex[ xyzv[VIO_Y] ];
            save_index2 =  sindex[ xyzv[VIO_Z] ];


            sindex[ xyzv[VIO_Y] ] -= 3;
            GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Y] ] += 2;
            GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Y] ] += 2;
            GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Y] ] += 2;
            GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Y] ] = save_index1 ;
            value1 = MY_CUBIC_05(v1,v2,v3,v4);


            sindex[ xyzv[VIO_Z] ] -= 3;
            GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Z] ] += 2;
            GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Z] ] += 2;
            GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Z] ] += 2;
            GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Z] ] = save_index2 ;
            value2 = MY_CUBIC_05(v1,v2,v3,v4);


            value1 = (value1 + value2)/2.0;

            SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
            count++;
          }
          update_progress_report( &progress, count+1 );
        }

      } /* index[ orig_xyzv[VIO_Z+1] ] */
      
    }
                                /* loop over all Y planes  */

    for(index[orig_xyzv[VIO_Y]]=0; index[orig_xyzv[VIO_Y]]<orig_count[orig_xyzv[VIO_Y]]-1; index[orig_xyzv[VIO_Y]]++) {

      sindex[ xyzv[VIO_Y] ] = 2*index[ orig_xyzv[VIO_Y] ];

      for(index[orig_xyzv[VIO_Z+1]]=0; index[orig_xyzv[VIO_Z+1]]<orig_count[xyzv[VIO_Z+1]]; index[orig_xyzv[VIO_Z+1]]++) {

        sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];


                                /* do faces near the edge first */

        for(index[orig_xyzv[VIO_X]]=0; index[orig_xyzv[VIO_X]]<orig_count[orig_xyzv[VIO_X]]-1; index[orig_xyzv[VIO_X]]++) {

          sindex[ xyzv[VIO_X] ] = index[ orig_xyzv[VIO_X] ]*2 + 1;
          
          sindex[ xyzv[VIO_Z] ] = 0;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Z] ] = 2;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          
          value1 = (v1+v2)/2;
          
          sindex[ xyzv[VIO_Z] ] = 1;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
          
          sindex[ xyzv[VIO_Z] ] = 2*orig_count[ orig_xyzv[VIO_Z] ] -3;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Z] ] = 2*orig_count[ orig_xyzv[VIO_Z] ] -1;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          
          value1 = (v1+v2)/2;
          
          sindex[ xyzv[VIO_Z] ] = 2*orig_count[ orig_xyzv[VIO_Z] ] -2;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
          
          
        }
        update_progress_report( &progress, count+1 );

        for(index[orig_xyzv[VIO_Z]]=0; index[orig_xyzv[VIO_Z]]<orig_count[orig_xyzv[VIO_Z]]-1; index[orig_xyzv[VIO_Z]]++) {

          sindex[ xyzv[VIO_Z] ] = index[ orig_xyzv[VIO_Z] ]*2 + 1;
          
          sindex[ xyzv[VIO_X] ] = 0;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] = 2;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          
          value1 = (v1+v2)/2;
          
          sindex[ xyzv[VIO_X] ] = 1;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
          
          sindex[ xyzv[VIO_X] ] = 2*orig_count[ orig_xyzv[VIO_X] ] -3;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] = 2*orig_count[ orig_xyzv[VIO_X] ] -1;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          
          value1 = (v1+v2)/2;
          
          sindex[ xyzv[VIO_X] ] = 2*orig_count[ orig_xyzv[VIO_X] ] -2;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
          
          
        }
        update_progress_report( &progress, count+1 );

                                /* now do faces in the middle */


        for(index[orig_xyzv[VIO_Z]]=1; index[orig_xyzv[VIO_Z]]<orig_count[orig_xyzv[VIO_Z]]-2; index[orig_xyzv[VIO_Z]]++) {
          for(index[orig_xyzv[VIO_X]]=1; index[orig_xyzv[VIO_X]]<orig_count[orig_xyzv[VIO_X]]-2; index[orig_xyzv[VIO_X]]++) {

            sindex[ xyzv[VIO_Z] ] = index[ orig_xyzv[VIO_Z] ]*2 + 1;
            sindex[ xyzv[VIO_X] ] = index[ orig_xyzv[VIO_X] ]*2 + 1;

            save_index1 =  sindex[ xyzv[VIO_Z] ];
            save_index2 =  sindex[ xyzv[VIO_X] ];


            sindex[ xyzv[VIO_Z] ] -= 3;
            GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Z] ] += 2;
            GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Z] ] += 2;
            GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Z] ] += 2;
            GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Z] ] = save_index1 ;
            value1 = MY_CUBIC_05(v1,v2,v3,v4);


            sindex[ xyzv[VIO_X] ] -= 3;
            GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_X] ] += 2;
            GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_X] ] += 2;
            GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_X] ] += 2;
            GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_X] ] = save_index2 ;
            value2 = MY_CUBIC_05(v1,v2,v3,v4);


            value1 = (value1 + value2)/2.0;

            SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
            count++;
          }
          update_progress_report( &progress, count+1 );
        }

      } /* index[ orig_xyzv[VIO_Z+1] ] */
      
    }

                                /* loop over all Z planes  */
    
    for(index[orig_xyzv[VIO_Z]]=0; index[orig_xyzv[VIO_Z]]<orig_count[orig_xyzv[VIO_Z]]-1; index[orig_xyzv[VIO_Z]]++) {

      sindex[ xyzv[VIO_Z] ] = 2*index[ orig_xyzv[VIO_Z] ];

      for(index[orig_xyzv[VIO_Z+1]]=0; index[orig_xyzv[VIO_Z+1]]<orig_count[xyzv[VIO_Z+1]]; index[orig_xyzv[VIO_Z+1]]++) {

        sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];


                                /* do faces near the edge first */

        for(index[orig_xyzv[VIO_Y]]=0; index[orig_xyzv[VIO_Y]]<orig_count[orig_xyzv[VIO_Y]]-1; index[orig_xyzv[VIO_Y]]++) {

          sindex[ xyzv[VIO_Y] ] = index[ orig_xyzv[VIO_Y] ]*2 + 1;
          
          sindex[ xyzv[VIO_X] ] = 0;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] = 2;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          
          value1 = (v1+v2)/2;
          
          sindex[ xyzv[VIO_X] ] = 1;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
          
          sindex[ xyzv[VIO_X] ] = 2*orig_count[ orig_xyzv[VIO_X] ] -3;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] = 2*orig_count[ orig_xyzv[VIO_X] ] -1;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          
          value1 = (v1+v2)/2;
          
          sindex[ xyzv[VIO_X] ] = 2*orig_count[ orig_xyzv[VIO_X] ] -2;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
          
          
        }
        update_progress_report( &progress, count+1 );

        for(index[orig_xyzv[VIO_X]]=0; index[orig_xyzv[VIO_X]]<orig_count[orig_xyzv[VIO_X]]-1; index[orig_xyzv[VIO_X]]++) {

          sindex[ xyzv[VIO_X] ] = index[ orig_xyzv[VIO_X] ]*2 + 1;
          
          sindex[ xyzv[VIO_Y] ] = 0;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] = 2;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          
          value1 = (v1+v2)/2;
          
          sindex[ xyzv[VIO_Y] ] = 1;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
          
          sindex[ xyzv[VIO_Y] ] = 2*orig_count[ orig_xyzv[VIO_Y] ] -3;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] = 2*orig_count[ orig_xyzv[VIO_Y] ] -1;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          
          value1 = (v1+v2)/2;
          
          sindex[ xyzv[VIO_Y] ] = 2*orig_count[ orig_xyzv[VIO_Y] ] -2;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
          
          
        }
        update_progress_report( &progress, count+1 );

                                /* now do faces in the middle */


        for(index[orig_xyzv[VIO_X]]=1; index[orig_xyzv[VIO_X]]<orig_count[orig_xyzv[VIO_X]]-2; index[orig_xyzv[VIO_X]]++) {
          for(index[orig_xyzv[VIO_Y]]=1; index[orig_xyzv[VIO_Y]]<orig_count[orig_xyzv[VIO_Y]]-2; index[orig_xyzv[VIO_Y]]++) {

            sindex[ xyzv[VIO_X] ] = index[ orig_xyzv[VIO_X] ]*2 + 1;
            sindex[ xyzv[VIO_Y] ] = index[ orig_xyzv[VIO_Y] ]*2 + 1;

            save_index1 =  sindex[ xyzv[VIO_X] ];
            save_index2 =  sindex[ xyzv[VIO_Y] ];


            sindex[ xyzv[VIO_X] ] -= 3;
            GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_X] ] += 2;
            GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_X] ] += 2;
            GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_X] ] += 2;
            GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_X] ] = save_index1 ;
            value1 = MY_CUBIC_05(v1,v2,v3,v4);


            sindex[ xyzv[VIO_Y] ] -= 3;
            GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Y] ] += 2;
            GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Y] ] += 2;
            GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Y] ] += 2;
            GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Y] ] = save_index2 ;
            value2 = MY_CUBIC_05(v1,v2,v3,v4);


            value1 = (value1 + value2)/2.0;

            SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
            count++;
          }
          update_progress_report( &progress, count+1 );
        }

      } /* index[ orig_xyzv[VIO_Z+1] ] */
      
    }
    /* LEVEL 3: center interpolation, identified as 'c' in desc above */

    for(index[orig_xyzv[VIO_Z+1]]=0; index[orig_xyzv[VIO_Z+1]]<orig_count[xyzv[VIO_Z+1]]; index[orig_xyzv[VIO_Z+1]]++) {

      sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];

                         /* do all centers that are 1 sample from
                            the edge of the super-sampled volume */

                                /* do the two X-planes */

      for(index[orig_xyzv[VIO_Y]]=0; index[orig_xyzv[VIO_Y]]<orig_count[orig_xyzv[VIO_Y]]-2; index[orig_xyzv[VIO_Y]]++) {
        for(index[orig_xyzv[VIO_Z]]=0; index[orig_xyzv[VIO_Z]]<orig_count[orig_xyzv[VIO_Z]]-2; index[orig_xyzv[VIO_Z]]++) {

          sindex[ xyzv[VIO_Y] ] = index[ orig_xyzv[VIO_Y] ]*2 + 1;
          sindex[ xyzv[VIO_Z] ] = index[ orig_xyzv[VIO_Z] ]*2 + 1;

          sindex[ xyzv[VIO_X] ] = 1;

          sindex[ xyzv[VIO_X] ] -= 1;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] += 2;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] -= 1;
          
          sindex[ xyzv[VIO_Y] ] -= 1;
          GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] += 2;
          GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] -= 1;
        
          sindex[ xyzv[VIO_Z] ] -= 1;
          GET_VOXEL_4D(v5, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Z] ] += 2;
          GET_VOXEL_4D(v6, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Z] ] -= 1;
        
          value1 = (v1 + v2 + v3 + v4 + v5 + v6) / 6.0;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
        
          sindex[ xyzv[VIO_X] ] = 2*orig_count[ orig_xyzv[VIO_X] ] -2;
      
          sindex[ xyzv[VIO_X] ] -= 1;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] += 2;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] -= 1;
          
          sindex[ xyzv[VIO_Y] ] -= 1;
          GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] += 2;
          GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] -= 1;
          
          sindex[ xyzv[VIO_Z] ] -= 1;
          GET_VOXEL_4D(v5, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Z] ] += 2;
          GET_VOXEL_4D(v6, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Z] ] -= 1;
          
          value1 = (v1 + v2 + v3 + v4 + v5 + v6) / 6.0;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
        }
      }
      update_progress_report( &progress, count+1 );
                                /* do the two Y-planes */

      for(index[orig_xyzv[VIO_X]]=0; index[orig_xyzv[VIO_X]]<orig_count[orig_xyzv[VIO_X]]-2; index[orig_xyzv[VIO_X]]++) {
        for(index[orig_xyzv[VIO_Z]]=0; index[orig_xyzv[VIO_Z]]<orig_count[orig_xyzv[VIO_Z]]-2; index[orig_xyzv[VIO_Z]]++) {

          sindex[ xyzv[VIO_X] ] = index[ orig_xyzv[VIO_X] ]*2 + 1;
          sindex[ xyzv[VIO_Z] ] = index[ orig_xyzv[VIO_Z] ]*2 + 1;

          sindex[ xyzv[VIO_Y] ] = 1;

          sindex[ xyzv[VIO_X] ] -= 1;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] += 2;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] -= 1;
          
          sindex[ xyzv[VIO_Y] ] -= 1;
          GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] += 2;
          GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] -= 1;
        
          sindex[ xyzv[VIO_Z] ] -= 1;
          GET_VOXEL_4D(v5, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Z] ] += 2;
          GET_VOXEL_4D(v6, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Z] ] -= 1;
        
          value1 = (v1 + v2 + v3 + v4 + v5 + v6) / 6.0;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
        
           sindex[ xyzv[VIO_Y] ] = 2*orig_count[ orig_xyzv[VIO_Y] ] -2;
      
          sindex[ xyzv[VIO_X] ] -= 1;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] += 2;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] -= 1;
          
          sindex[ xyzv[VIO_Y] ] -= 1;
          GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] += 2;
          GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] -= 1;
          
          sindex[ xyzv[VIO_Z] ] -= 1;
          GET_VOXEL_4D(v5, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Z] ] += 2;
          GET_VOXEL_4D(v6, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Z] ] -= 1;
          
          value1 = (v1 + v2 + v3 + v4 + v5 + v6) / 6.0;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
        }
      }
      update_progress_report( &progress, count+1 );

                                /* do the two Z-planes */

      for(index[orig_xyzv[VIO_X]]=0; index[orig_xyzv[VIO_X]]<orig_count[orig_xyzv[VIO_X]]-2; index[orig_xyzv[VIO_X]]++) {
        for(index[orig_xyzv[VIO_Y]]=0; index[orig_xyzv[VIO_Y]]<orig_count[orig_xyzv[VIO_Y]]-2; index[orig_xyzv[VIO_Y]]++) {

          sindex[ xyzv[VIO_X] ] = index[ orig_xyzv[VIO_X] ]*2 + 1;
          sindex[ xyzv[VIO_Y] ] = index[ orig_xyzv[VIO_Y] ]*2 + 1;

          sindex[ xyzv[VIO_Z] ] = 1;

          sindex[ xyzv[VIO_X] ] -= 1;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] += 2;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] -= 1;
          
          sindex[ xyzv[VIO_Y] ] -= 1;
          GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] += 2;
          GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] -= 1;
        
          sindex[ xyzv[VIO_Z] ] -= 1;
          GET_VOXEL_4D(v5, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Z] ] += 2;
          GET_VOXEL_4D(v6, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Z] ] -= 1;
        
          value1 = (v1 + v2 + v3 + v4 + v5 + v6) / 6.0;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
        
          sindex[ xyzv[VIO_Z] ] = 2*orig_count[ orig_xyzv[VIO_Z] ] -2;
      
          sindex[ xyzv[VIO_X] ] -= 1;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] += 2;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] -= 1;
          
          sindex[ xyzv[VIO_Y] ] -= 1;
          GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] += 2;
          GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] -= 1;
          
          sindex[ xyzv[VIO_Z] ] -= 1;
          GET_VOXEL_4D(v5, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Z] ] += 2;
          GET_VOXEL_4D(v6, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Z] ] -= 1;
          
          value1 = (v1 + v2 + v3 + v4 + v5 + v6) / 6.0;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
        }
      }
      update_progress_report( &progress, count+1 );

                               /* now do all central voxels */


      for(index[orig_xyzv[VIO_X]]=1; index[orig_xyzv[VIO_X]]<orig_count[orig_xyzv[VIO_X]]-2; index[orig_xyzv[VIO_X]]++) {
        for(index[orig_xyzv[VIO_Y]]=1; index[orig_xyzv[VIO_Y]]<orig_count[orig_xyzv[VIO_Y]]-2; index[orig_xyzv[VIO_Y]]++) {
          for(index[orig_xyzv[VIO_Z]]=1; index[orig_xyzv[VIO_Z]]<orig_count[orig_xyzv[VIO_Z]]-2; index[orig_xyzv[VIO_Z]]++) {
            
            sindex[ xyzv[VIO_X] ] = index[ orig_xyzv[VIO_X] ]*2 + 1;
            sindex[ xyzv[VIO_Y] ] = index[ orig_xyzv[VIO_Y] ]*2 + 1;
            sindex[ xyzv[VIO_Z] ] = index[ orig_xyzv[VIO_Z] ]*2 + 1;

            sindex[ xyzv[VIO_X] ] -= 3;
            GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_X] ] += 2;
            GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_X] ] += 2;
            GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_X] ] += 2;
            GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_X] ] -= 3;
            value1 = MY_CUBIC_05(v1,v2,v3,v4);
            sindex[ xyzv[VIO_Y] ] -= 3;
            GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Y] ] += 2;
            GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Y] ] += 2;
            GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Y] ] += 2;
            GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Y] ] -= 3;
             value2 = MY_CUBIC_05(v1,v2,v3,v4);
            
            sindex[ xyzv[VIO_Z] ] -= 3;
            GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Z] ] += 2;
            GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Z] ] += 2;
            GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Z] ] += 2;
            GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
              sindex[ xyzv[VIO_Z] ] -= 3;
            value3 = MY_CUBIC_05(v1,v2,v3,v4);
            
            value1 = (value1 + value2 + value3) / 3.0;
    
            SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
            count++;

            update_progress_report( &progress, count+1 );
          }
        }
      }
    }  /* index[ orig_xyzv[VIO_Z+1] ] */

    terminate_progress_report( &progress );

}


void interpolate_super_sampled_data_by2(
    VIO_General_transform *orig_deformation,
    VIO_General_transform *super_sampled)
{

  int
    i,num_dim,
    xyzv[VIO_MAX_DIMENSIONS],
    count[VIO_MAX_DIMENSIONS];
  VIO_Volume
    orig_vol;


  

  if (orig_deformation->type != GRID_TRANSFORM || super_sampled->type != GRID_TRANSFORM) {
    print_error_and_line_num("interpolate_super_sampled_data_by2 not called with GRID_TRANSFORM",
                             __FILE__, __LINE__);
  }

  orig_vol = orig_deformation->displacement_volume;
  get_volume_sizes(       orig_vol,  count);
  get_volume_XYZV_indices(orig_vol,  xyzv);

  num_dim = 0;

  for(i=0; i<VIO_N_DIMENSIONS; i++) {
    if ( count[xyzv[i]] > 1 ) 
      num_dim++;
  }
  if (num_dim == 3) {
    interpolate_super_sampled_data_by2_dim3( orig_deformation, super_sampled );
  } else {

    if (num_dim == 2){                /* then one dim == 1 */
      
      for(i=0; i<VIO_N_DIMENSIONS; i++) {

        if ( count[xyzv[i]] == 1 ) {
          interpolate_super_sampled_data_by2_dim2( orig_deformation, super_sampled,i );
        }

      }
    }
  }
}


static void interpolate_super_sampled_data_by2_dim2(VIO_General_transform *orig_deformation,
                                                     VIO_General_transform *super_sampled,
                                                     int i)
{

  if(i==0) interpolate_super_sampled_data_by2_dim2_X(orig_deformation,super_sampled );
  else
    { 
      if (i==1) interpolate_super_sampled_data_by2_dim2_Y(orig_deformation,super_sampled );
      else
        {
          if (i==2) interpolate_super_sampled_data_by2_dim2_Z(orig_deformation,super_sampled );
          else print("\n\n\n WARNING\n in interpolate_super_sampled_data_by2_dim2 \n i >2 \n\n");
        }
    }
}


static void interpolate_super_sampled_data_by2_dim2_X( 
    VIO_General_transform *orig_deformation,
    VIO_General_transform *super_sampled )
{
  VIO_Volume
    orig_vol,
    super_vol;

  int 
    i, save_index, save_index1, save_index2,
    count,
    index[VIO_MAX_DIMENSIONS],
    sindex[VIO_MAX_DIMENSIONS],
    orig_xyzv[VIO_MAX_DIMENSIONS],
    xyzv[VIO_MAX_DIMENSIONS],
    orig_count[VIO_MAX_DIMENSIONS],
    super_count[VIO_MAX_DIMENSIONS];
  VIO_Real
    value1, value2, value3, v1,v2,v3,v4,v5,v6;
  VIO_progress_struct
    progress;

    orig_vol = orig_deformation->displacement_volume;
    super_vol= super_sampled->displacement_volume;
    
    get_volume_sizes(       super_vol, super_count);
    get_volume_sizes(       orig_vol,  orig_count);
    get_volume_XYZV_indices(super_vol, xyzv);
    get_volume_XYZV_indices(orig_vol,  orig_xyzv);
    
    init_the_volume_to_zero(super_sampled->displacement_volume);

    initialize_progress_report(&progress, FALSE, 
                                      super_count[ xyzv[VIO_Y] ]*
                                super_count[ xyzv[VIO_Z] ]*super_count[ xyzv[VIO_Z+1] ]+1,
                               "Super-sampling defs:" );
    count = 0;


 

    /* LEVEL 0: copy original 'corner' nodes, identified as 'X' in desc above */

   
    for(i=0; i<VIO_MAX_DIMENSIONS; i++) sindex[i]=index[i]=0;
    
      for(index[orig_xyzv[VIO_Y]]=0; index[orig_xyzv[VIO_Y]]<orig_count[orig_xyzv[VIO_Y]]; index[orig_xyzv[VIO_Y]]++) {
        for(index[orig_xyzv[VIO_Z]]=0; index[orig_xyzv[VIO_Z]]<orig_count[orig_xyzv[VIO_Z]]; index[orig_xyzv[VIO_Z]]++) {

          for(i=1; i<VIO_N_DIMENSIONS; i++)
            sindex[ xyzv[i] ] = 2*index[ orig_xyzv[i] ];

          for(index[orig_xyzv[VIO_Z+1]]=0; index[orig_xyzv[VIO_Z+1]]<orig_count[xyzv[VIO_Z+1]]; index[orig_xyzv[VIO_Z+1]]++) {
            sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];
            GET_VOXEL_4D(value1, orig_vol, index[0],index[1],index[2],index[3]);
            SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          }          

          count += 3;
        }
      }
      update_progress_report( &progress, count+1 );



   /* LEVEL 1: edge interpolation, identified as 'e' in desc above */

 
               /* do edges along the index[ orig_xyzv[VIO_Y] ] dir */

    for(i=0; i<VIO_MAX_DIMENSIONS; i++) sindex[i]=index[i]=0;
    
                                /* loop over all x-dirs  */
    
    
    for(index[orig_xyzv[VIO_Z]]=0; index[orig_xyzv[VIO_Z]]<orig_count[orig_xyzv[VIO_Z]]; index[orig_xyzv[VIO_Z]]++) {

        sindex[ xyzv[VIO_Z] ] = 2*index[ orig_xyzv[VIO_Z] ];


        for(index[orig_xyzv[VIO_Z+1]]=0; index[orig_xyzv[VIO_Z+1]]<orig_count[xyzv[VIO_Z+1]]; index[orig_xyzv[VIO_Z+1]]++) {

                                /* do linear interp at ends */

          sindex[ xyzv[VIO_Y] ] = 1;                             /* beginning end */
          index[ orig_xyzv[VIO_Y] ]=0;             
          GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
          index[ orig_xyzv[VIO_Y] ]=1;
          GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
          value1 = (v1+v2)/2;

          sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];

          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;

          sindex[ xyzv[VIO_Y] ] = 2*orig_count[ orig_xyzv[VIO_Y] ] -3;      /* ending end */
          index[ orig_xyzv[VIO_Y] ]= orig_count[ orig_xyzv[VIO_Y] ]-2;              
          GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
          index[ orig_xyzv[VIO_Y] ]= orig_count[ orig_xyzv[VIO_Y] ]-1;
          GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
          value1 = (v1+v2)/2;

          sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];

          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
                                /* now do the voxels between the two ends */


          for(index[orig_xyzv[VIO_Y]]=1; index[orig_xyzv[VIO_Y]]<=orig_count[orig_xyzv[VIO_Y]]-3; index[orig_xyzv[VIO_Y]]++) {

            save_index = index[ orig_xyzv[VIO_Y] ];
            sindex[ xyzv[VIO_Y] ] = index[ orig_xyzv[VIO_Y] ]*2 + 1;

            index[ orig_xyzv[VIO_Y] ]--;             
            GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_Y] ]++;             
            GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_Y] ]++;             
            GET_VOXEL_4D(v3, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_Y] ]++;             
            GET_VOXEL_4D(v4, orig_vol, index[0],index[1],index[2],index[3]);

            value1 = MY_CUBIC_05(v1,v2,v3,v4);

            SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
            count++;

            index[ orig_xyzv[VIO_Y] ] = save_index;
          }
        }          
      update_progress_report( &progress, count+1 );
    }
               /* do edges along the index[ orig_xyzv[VIO_Z] ] dir */

    for(i=0; i<VIO_MAX_DIMENSIONS; i++) sindex[i]=index[i]=0;
    
                                /* loop over all Z-dirs  */


      for(index[orig_xyzv[VIO_Y]]=0; index[orig_xyzv[VIO_Y]]<orig_count[orig_xyzv[VIO_Y]]; index[orig_xyzv[VIO_Y]]++) {

        sindex[ xyzv[VIO_Y] ] = 2*index[ orig_xyzv[VIO_Y] ];

        for(index[orig_xyzv[VIO_Z+1]]=0; index[orig_xyzv[VIO_Z+1]]<orig_count[xyzv[VIO_Z+1]]; index[orig_xyzv[VIO_Z+1]]++) {

                                /* do linear interp at ends */

          sindex[ xyzv[VIO_Z] ] = 1;                             /* beginning end */
          index[ orig_xyzv[VIO_Z] ]=0;             
          GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
          index[ orig_xyzv[VIO_Z] ]=1;
          GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
          value1 = (v1+v2)/2;

          sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];

          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;

          sindex[ xyzv[VIO_Z] ] = 2*orig_count[ xyzv[VIO_Z] ] -3;      /* ending end */
          index[ orig_xyzv[VIO_Z] ]= orig_count[ orig_xyzv[VIO_Z] ]-2;              
          GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
          index[ orig_xyzv[VIO_Z] ]= orig_count[ orig_xyzv[VIO_Z] ]-1;
          GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
          value1 = (v1+v2)/2;

          sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];

          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
                                /* now do the voxels between the two ends */


          for(index[orig_xyzv[VIO_Z]]=1; index[orig_xyzv[VIO_Z]]<=orig_count[orig_xyzv[VIO_Z]]-3; index[orig_xyzv[VIO_Z]]++) {

            save_index = index[ orig_xyzv[VIO_Z] ];
            sindex[ xyzv[VIO_Z] ] = index[ orig_xyzv[VIO_Z] ]*2 + 1;

            index[ orig_xyzv[VIO_Z] ]--;             
            GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_Z] ]++;             
            GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_Z] ]++;             
            GET_VOXEL_4D(v3, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_Z] ]++;             
            GET_VOXEL_4D(v4, orig_vol, index[0],index[1],index[2],index[3]);

            value1 = MY_CUBIC_05(v1,v2,v3,v4);

            SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
            count++;

            index[ orig_xyzv[VIO_Z] ] = save_index;
          }

        }          

      update_progress_report( &progress, count+1 );
    }

  /* LEVEL 2: face interpolation, identified as 'f' in desc above */

               /* do faces in the plane with index[ orig_xyzv[VIO_X] ]=CONST */

    for(i=0; i<VIO_MAX_DIMENSIONS; i++)  sindex[i]=index[i]=0;
    
                                /* loop over the X plane  */

 
      for(index[orig_xyzv[VIO_Z+1]]=0; index[orig_xyzv[VIO_Z+1]]<orig_count[xyzv[VIO_Z+1]]; index[orig_xyzv[VIO_Z+1]]++) {

        sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];


                                /* do faces near the edge first */

        for(index[orig_xyzv[VIO_Z]]=0; index[orig_xyzv[VIO_Z]]<orig_count[orig_xyzv[VIO_Z]]-1; index[orig_xyzv[VIO_Z]]++) {

          sindex[ xyzv[VIO_Z] ] = index[ orig_xyzv[VIO_Z] ]*2 + 1;
          
          sindex[ xyzv[VIO_Y] ] = 0;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] = 2;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          
          value1 = (v1+v2)/2;
          
          sindex[ xyzv[VIO_Y] ] = 1;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
          
          sindex[ xyzv[VIO_Y] ] = 2*orig_count[ orig_xyzv[VIO_Y] ] -3;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] = 2*orig_count[ orig_xyzv[VIO_Y] ] -1;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          
          value1 = (v1+v2)/2;
          
          sindex[ xyzv[VIO_Y] ] = 2*orig_count[ orig_xyzv[VIO_Y] ] -2;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
          
          
        }
        update_progress_report( &progress, count+1 );

        for(index[orig_xyzv[VIO_Y]]=0; index[orig_xyzv[VIO_Y]]<=orig_count[orig_xyzv[VIO_Y]]-1; index[orig_xyzv[VIO_Y]]++) {

          sindex[ xyzv[VIO_Y] ] = index[ orig_xyzv[VIO_Y] ]*2 + 1;
          
          sindex[ xyzv[VIO_Z] ] = 0;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Z] ] = 2;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          
          value1 = (v1+v2)/2;
          
          sindex[ xyzv[VIO_Z] ] = 1;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
          
          sindex[ xyzv[VIO_Z] ] = 2*orig_count[ orig_xyzv[VIO_Z] ] -3;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Z] ] = 2*orig_count[ orig_xyzv[VIO_Z] ] -1;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          
          value1 = (v1+v2)/2;
          
          sindex[ xyzv[VIO_Z] ] = 2*orig_count[ orig_xyzv[VIO_Z] ] -2;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;

        }
        update_progress_report( &progress, count+1 );
      
                                /* now do faces in the middle */


        for(index[orig_xyzv[VIO_Y]]=1; index[orig_xyzv[VIO_Y]]<orig_count[orig_xyzv[VIO_Y]]-2; index[orig_xyzv[VIO_Y]]++) {
          for(index[orig_xyzv[VIO_Z]]=1; index[orig_xyzv[VIO_Z]]<orig_count[orig_xyzv[VIO_Z]]-2; index[orig_xyzv[VIO_Z]]++) {

            sindex[ xyzv[VIO_Y] ] = index[ orig_xyzv[VIO_Y] ]*2 + 1;
            sindex[ xyzv[VIO_Z] ] = index[ orig_xyzv[VIO_Z] ]*2 + 1;

            save_index1 =  sindex[ xyzv[VIO_Y] ];
            save_index2 =  sindex[ xyzv[VIO_Z] ];


            sindex[ xyzv[VIO_Y] ] -= 3;
            GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Y] ] += 2;
            GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Y] ] += 2;
            GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Y] ] += 2;
            GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Y] ] = save_index1 ;
            value1 = MY_CUBIC_05(v1,v2,v3,v4);


            sindex[ xyzv[VIO_Z] ] -= 3;
            GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Z] ] += 2;
            GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Z] ] += 2;
            GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Z] ] += 2;
            GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Z] ] = save_index2 ;
            value2 = MY_CUBIC_05(v1,v2,v3,v4);


            value1 = (value1 + value2)/2.0;

            SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
            count++;
          }
          update_progress_report( &progress, count+1 );
        }

      } /* index[ orig_xyzv[VIO_Z+1] ] */
      

    terminate_progress_report( &progress );

}


static void interpolate_super_sampled_data_by2_dim2_Y( 
    VIO_General_transform *orig_deformation,
    VIO_General_transform *super_sampled )
{
  VIO_Volume
    orig_vol,
    super_vol;

  int 
    i, save_index, save_index1, save_index2,
    count,
    index[VIO_MAX_DIMENSIONS],
    sindex[VIO_MAX_DIMENSIONS],
    orig_xyzv[VIO_MAX_DIMENSIONS],
    xyzv[VIO_MAX_DIMENSIONS],
    orig_count[VIO_MAX_DIMENSIONS],
    super_count[VIO_MAX_DIMENSIONS];
  VIO_Real
    value1, value2, value3, v1,v2,v3,v4,v5,v6;
  VIO_progress_struct
    progress;

    orig_vol = orig_deformation->displacement_volume;
    super_vol= super_sampled->displacement_volume;
    
    get_volume_sizes(       super_vol, super_count);
    get_volume_sizes(       orig_vol,  orig_count);
    get_volume_XYZV_indices(super_vol, xyzv);
    get_volume_XYZV_indices(orig_vol,  orig_xyzv);
    
    init_the_volume_to_zero(super_sampled->displacement_volume);

    initialize_progress_report(&progress, FALSE, 
                               super_count[ xyzv[VIO_X] ]*super_count[ xyzv[VIO_Y] ]*
                               super_count[ xyzv[VIO_Z] ]*super_count[ xyzv[VIO_Z+1] ]+1,
                               "Super-sampling defs:" );
    count = 0;



    /* LEVEL 0: copy original 'corner' nodes, identified as 'X' in desc above */

    for(i=0; i<VIO_MAX_DIMENSIONS; i++) sindex[i]=index[i]=0;
    
    for(index[orig_xyzv[VIO_X]]=0; index[orig_xyzv[VIO_X]]<orig_count[orig_xyzv[VIO_X]]; index[orig_xyzv[VIO_X]]++) {
        for(index[orig_xyzv[VIO_Z]]=0; index[orig_xyzv[VIO_Z]]<orig_count[orig_xyzv[VIO_Z]]; index[orig_xyzv[VIO_Z]]++) {

          for(i=0; i<VIO_N_DIMENSIONS; i++)
            sindex[ xyzv[i] ] = 2*index[ orig_xyzv[i] ];
          sindex[ xyzv[VIO_Y] ] = 0;

          for(index[orig_xyzv[VIO_Z+1]]=0; index[orig_xyzv[VIO_Z+1]]<orig_count[xyzv[VIO_Z+1]]; index[orig_xyzv[VIO_Z+1]]++) {
            sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];
            GET_VOXEL_4D(value1, orig_vol, index[0],index[1],index[2],index[3]);
            SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          }          

          count += 3;
        }
      update_progress_report( &progress, count+1 );
    }

    

    /* LEVEL 1: edge interpolation, identified as 'e' in desc above */

               /* do edges along the index[ orig_xyzv[VIO_X] ] dir */

    for(i=0; i<VIO_MAX_DIMENSIONS; i++) sindex[i]=index[i]=0;
    
                                /* loop over all x-dirs  */
    
    sindex[ xyzv[VIO_Y] ] = 0;

      for(index[orig_xyzv[VIO_Z]]=0; index[orig_xyzv[VIO_Z]]<orig_count[orig_xyzv[VIO_Z]]; index[orig_xyzv[VIO_Z]]++) {

        sindex[ xyzv[VIO_Z] ] = 2*index[ orig_xyzv[VIO_Z] ];

        for(index[orig_xyzv[VIO_Z+1]]=0; index[orig_xyzv[VIO_Z+1]]<orig_count[xyzv[VIO_Z+1]]; index[orig_xyzv[VIO_Z+1]]++) {

                                /* do linear interp at ends */

          sindex[ xyzv[VIO_X] ] = 1;                             /* beginning end */
          index[ orig_xyzv[VIO_X] ]=0;             
          GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
          index[ orig_xyzv[VIO_X] ]=1;
          GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
          value1 = (v1+v2)/2;

          sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];

          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;

          sindex[ xyzv[VIO_X] ] = 2*orig_count[ orig_xyzv[VIO_X] ] -3;      /* ending end */
          index[ orig_xyzv[VIO_X] ]= orig_count[ orig_xyzv[VIO_X] ]-2;              
          GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
          index[ orig_xyzv[VIO_X] ]= orig_count[ orig_xyzv[VIO_X] ]-1;
          GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
          value1 = (v1+v2)/2;

          sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];

          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;

                                /* now do the voxels between the two ends */


          for(index[orig_xyzv[VIO_X]]=1; index[orig_xyzv[VIO_X]]<=orig_count[orig_xyzv[VIO_X]]-3; index[orig_xyzv[VIO_X]]++) {

            save_index = index[ orig_xyzv[VIO_X] ];
            sindex[ xyzv[VIO_X] ] = index[ orig_xyzv[VIO_X] ]*2 + 1;

            index[ orig_xyzv[VIO_X] ]--;             
            GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_X] ]++;             
            GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_X] ]++;             
            GET_VOXEL_4D(v3, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_X] ]++;             
            GET_VOXEL_4D(v4, orig_vol, index[0],index[1],index[2],index[3]);

            value1 = MY_CUBIC_05(v1,v2,v3,v4);

            SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
            count++;

            index[ orig_xyzv[VIO_X] ] = save_index;
          }

        }          

      }
      update_progress_report( &progress, count+1 );


               /* do edges along the index[ orig_xyzv[VIO_Z] ] dir */

    for(i=0; i<VIO_MAX_DIMENSIONS; i++) sindex[i]=index[i]=0;
    
                                /* loop over all Z-dirs  */

    sindex[ xyzv[VIO_Y] ] = 0;

    for(index[orig_xyzv[VIO_X]]=0; index[orig_xyzv[VIO_X]]<orig_count[orig_xyzv[VIO_X]]; index[orig_xyzv[VIO_X]]++) {


        sindex[ xyzv[VIO_X] ] = 2*index[ orig_xyzv[VIO_X] ];


        for(index[orig_xyzv[VIO_Z+1]]=0; index[orig_xyzv[VIO_Z+1]]<orig_count[xyzv[VIO_Z+1]]; index[orig_xyzv[VIO_Z+1]]++) {

                                /* do linear interp at ends */

          sindex[ xyzv[VIO_Z] ] = 1;                             /* beginning end */
          index[ orig_xyzv[VIO_Z] ]=0;             
          GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
          index[ orig_xyzv[VIO_Z] ]=1;
          GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
          value1 = (v1+v2)/2;

          sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];

          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;

          sindex[ xyzv[VIO_Z] ] = 2*orig_count[ xyzv[VIO_Z] ] -3;      /* ending end */
          index[ orig_xyzv[VIO_Z] ]= orig_count[ orig_xyzv[VIO_Z] ]-2;              
          GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
          index[ orig_xyzv[VIO_Z] ]= orig_count[ orig_xyzv[VIO_Z] ]-1;
          GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
          value1 = (v1+v2)/2;

          sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];

          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;


                                /* now do the voxels between the two ends */


          for(index[orig_xyzv[VIO_Z]]=1; index[orig_xyzv[VIO_Z]]<=orig_count[orig_xyzv[VIO_Z]]-3; index[orig_xyzv[VIO_Z]]++) {

            save_index = index[ orig_xyzv[VIO_Z] ];
            sindex[ xyzv[VIO_Z] ] = index[ orig_xyzv[VIO_Z] ]*2 + 1;

            index[ orig_xyzv[VIO_Z] ]--;             
            GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_Z] ]++;             
            GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_Z] ]++;             
            GET_VOXEL_4D(v3, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_Z] ]++;             
            GET_VOXEL_4D(v4, orig_vol, index[0],index[1],index[2],index[3]);

            value1 = MY_CUBIC_05(v1,v2,v3,v4);

            SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
            count++;

            index[ orig_xyzv[VIO_Z] ] = save_index;
          }

        }          
      update_progress_report( &progress, count+1 );
    }


    /* LEVEL 2: face interpolation, identified as 'f' in desc above */

               /* do faces in the plane with index[ orig_xyzv[VIO_X] ]=CONST */

    for(i=0; i<VIO_MAX_DIMENSIONS; i++)  sindex[i]=index[i]=0;
    

                        /* loop over the Y plane  */


      sindex[ xyzv[VIO_Y] ] = 0;

      for(index[orig_xyzv[VIO_Z+1]]=0; index[orig_xyzv[VIO_Z+1]]<orig_count[xyzv[VIO_Z+1]]; index[orig_xyzv[VIO_Z+1]]++) {

        sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];


                                /* do faces near the edge first */

        for(index[orig_xyzv[VIO_X]]=0; index[orig_xyzv[VIO_X]]<orig_count[orig_xyzv[VIO_X]]-1; index[orig_xyzv[VIO_X]]++) {

          sindex[ xyzv[VIO_X] ] = index[ orig_xyzv[VIO_X] ]*2 + 1;
          
          sindex[ xyzv[VIO_Z] ] = 0;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Z] ] = 2;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          
          value1 = (v1+v2)/2;
          
          sindex[ xyzv[VIO_Z] ] = 1;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
          
          sindex[ xyzv[VIO_Z] ] = 2*orig_count[ orig_xyzv[VIO_Z] ] -3;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Z] ] = 2*orig_count[ orig_xyzv[VIO_Z] ] -1;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          
          value1 = (v1+v2)/2;
          
          sindex[ xyzv[VIO_Z] ] = 2*orig_count[ orig_xyzv[VIO_Z] ] -2;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
          
          
        }
        update_progress_report( &progress, count+1 );

        for(index[orig_xyzv[VIO_Z]]=0; index[orig_xyzv[VIO_Z]]<orig_count[orig_xyzv[VIO_Z]]-1; index[orig_xyzv[VIO_Z]]++) {

          sindex[ xyzv[VIO_Z] ] = index[ orig_xyzv[VIO_Z] ]*2 + 1;
          
          sindex[ xyzv[VIO_X] ] = 0;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] = 2;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          
          value1 = (v1+v2)/2;
          
          sindex[ xyzv[VIO_X] ] = 1;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
          
          sindex[ xyzv[VIO_X] ] = 2*orig_count[ orig_xyzv[VIO_X] ] -3;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] = 2*orig_count[ orig_xyzv[VIO_X] ] -1;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          
          value1 = (v1+v2)/2;
          
          sindex[ xyzv[VIO_X] ] = 2*orig_count[ orig_xyzv[VIO_X] ] -2;
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
          
          
        }
        update_progress_report( &progress, count+1 );

                                /* now do faces in the middle */


        for(index[orig_xyzv[VIO_Z]]=1; index[orig_xyzv[VIO_Z]]<orig_count[orig_xyzv[VIO_Z]]-2; index[orig_xyzv[VIO_Z]]++) {
          for(index[orig_xyzv[VIO_X]]=1; index[orig_xyzv[VIO_X]]<orig_count[orig_xyzv[VIO_X]]-2; index[orig_xyzv[VIO_X]]++) {

            sindex[ xyzv[VIO_Z] ] = index[ orig_xyzv[VIO_Z] ]*2 + 1;
            sindex[ xyzv[VIO_X] ] = index[ orig_xyzv[VIO_X] ]*2 + 1;

            save_index1 =  sindex[ xyzv[VIO_Z] ];
            save_index2 =  sindex[ xyzv[VIO_X] ];


            sindex[ xyzv[VIO_Z] ] -= 3;
            GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Z] ] += 2;
            GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Z] ] += 2;
            GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Z] ] += 2;
            GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_Z] ] = save_index1 ;
            value1 = MY_CUBIC_05(v1,v2,v3,v4);


            sindex[ xyzv[VIO_X] ] -= 3;
            GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_X] ] += 2;
            GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_X] ] += 2;
            GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_X] ] += 2;
            GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
            sindex[ xyzv[VIO_X] ] = save_index2 ;
            value2 = MY_CUBIC_05(v1,v2,v3,v4);


            value1 = (value1 + value2)/2.0;

            SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
            count++;
          }
          update_progress_report( &progress, count+1 );
        }

      } /* index[ orig_xyzv[VIO_Z+1] ] */
      

    terminate_progress_report( &progress );

}

static void interpolate_super_sampled_data_by2_dim2_Z( 
    VIO_General_transform *orig_deformation,
    VIO_General_transform *super_sampled )
{
  VIO_Volume
    orig_vol,
    super_vol;

  int 
    i, save_index, save_index1, save_index2,
    count,
    index[VIO_MAX_DIMENSIONS],
    sindex[VIO_MAX_DIMENSIONS],
    orig_xyzv[VIO_MAX_DIMENSIONS],
    xyzv[VIO_MAX_DIMENSIONS],
    orig_count[VIO_MAX_DIMENSIONS],
    super_count[VIO_MAX_DIMENSIONS];
  VIO_Real
    value1, value2, value3, v1,v2,v3,v4,v5,v6;
  VIO_progress_struct
    progress;

    orig_vol = orig_deformation->displacement_volume;
    super_vol= super_sampled->displacement_volume;
    
    get_volume_sizes(       super_vol, super_count);
    get_volume_sizes(       orig_vol,  orig_count);
    get_volume_XYZV_indices(super_vol, xyzv);
    get_volume_XYZV_indices(orig_vol,  orig_xyzv);
    
    init_the_volume_to_zero(super_sampled->displacement_volume);

    initialize_progress_report(&progress, FALSE, 
                               super_count[ xyzv[VIO_X] ]*super_count[ xyzv[VIO_Y] ]*
                               super_count[ xyzv[VIO_Z] ]*super_count[ xyzv[VIO_Z+1] ]+1,
                               "Super-sampling defs:" );
    count = 0;

 
    /* LEVEL 0: copy original 'corner' nodes, identified as 'X' in desc above */

    for(i=0; i<VIO_MAX_DIMENSIONS; i++) sindex[i]=index[i]=0;
    
    for(index[orig_xyzv[VIO_X]]=0; index[orig_xyzv[VIO_X]]<orig_count[orig_xyzv[VIO_X]]; index[orig_xyzv[VIO_X]]++) {
      for(index[orig_xyzv[VIO_Y]]=0; index[orig_xyzv[VIO_Y]]<orig_count[orig_xyzv[VIO_Y]]; index[orig_xyzv[VIO_Y]]++) {

          for(i=0; i<VIO_N_DIMENSIONS; i++)
            sindex[ xyzv[i] ] = 2*index[ orig_xyzv[i] ];
 
          for(index[orig_xyzv[VIO_Z+1]]=0; index[orig_xyzv[VIO_Z+1]]<orig_count[xyzv[VIO_Z+1]]; index[orig_xyzv[VIO_Z+1]]++) {
            sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];
            GET_VOXEL_4D(value1, orig_vol, index[0],index[1],index[2],index[3]);
            SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          }          


      }
      update_progress_report( &progress, count+1 );
    }


    /* LEVEL 1: edge interpolation, identified as 'e' in desc above */

               /* do edges along the index[ orig_xyzv[VIO_X] ] dir */

    for(i=0; i<VIO_MAX_DIMENSIONS; i++) sindex[i]=index[i]=0;
    
                                /* loop over all x-dirs  */

    sindex[ xyzv[VIO_Z] ] = 0;
    for(index[orig_xyzv[VIO_Y]]=0; index[orig_xyzv[VIO_Y]]<orig_count[orig_xyzv[VIO_Y]]; index[orig_xyzv[VIO_Y]]++) {

        sindex[ xyzv[VIO_Y] ] = 2*index[ orig_xyzv[VIO_Y] ];
        

        for(index[orig_xyzv[VIO_Z+1]]=0; index[orig_xyzv[VIO_Z+1]]<orig_count[xyzv[VIO_Z+1]]; index[orig_xyzv[VIO_Z+1]]++) {

                                /* do linear interp at ends */

          sindex[ xyzv[VIO_X] ] = 1;                             /* beginning end */
          index[ orig_xyzv[VIO_X] ]=0;             
          GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
          index[ orig_xyzv[VIO_X] ]=1;
          GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
          value1 = (v1+v2)/2;

          sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];

          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;

          sindex[ xyzv[VIO_X] ] = 2*orig_count[ orig_xyzv[VIO_X] ] -3;      /* ending end */
          index[ orig_xyzv[VIO_X] ]= orig_count[ orig_xyzv[VIO_X] ]-2;              
          GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
          index[ orig_xyzv[VIO_X] ]= orig_count[ orig_xyzv[VIO_X] ]-1;
          GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
          value1 = (v1+v2)/2;

          sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];

          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;

                                /* now do the voxels between the two ends */


          for(index[orig_xyzv[VIO_X]]=1; index[orig_xyzv[VIO_X]]<=orig_count[orig_xyzv[VIO_X]]-3; index[orig_xyzv[VIO_X]]++) {

            save_index = index[ orig_xyzv[VIO_X] ];
            sindex[ xyzv[VIO_X] ] = index[ orig_xyzv[VIO_X] ]*2 + 1;

            index[ orig_xyzv[VIO_X] ]--;             
            GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_X] ]++;             
            GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_X] ]++;             
            GET_VOXEL_4D(v3, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_X] ]++;             
            GET_VOXEL_4D(v4, orig_vol, index[0],index[1],index[2],index[3]);

            value1 = MY_CUBIC_05(v1,v2,v3,v4);

            SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
            count++;

            index[ orig_xyzv[VIO_X] ] = save_index;
          }

        }          

      update_progress_report( &progress, count+1 );
    }
               /* do edges along the index[ orig_xyzv[VIO_Y] ] dir */

    for(i=0; i<VIO_MAX_DIMENSIONS; i++) sindex[i]=index[i]=0;
    
                                /* loop over all x-dirs  */


    for(index[orig_xyzv[VIO_X]]=0; index[orig_xyzv[VIO_X]]<orig_count[orig_xyzv[VIO_X]]; index[orig_xyzv[VIO_X]]++) {



        sindex[ xyzv[VIO_X] ] = 2*index[ orig_xyzv[VIO_X] ];

        for(index[orig_xyzv[VIO_Z+1]]=0; index[orig_xyzv[VIO_Z+1]]<orig_count[xyzv[VIO_Z+1]]; index[orig_xyzv[VIO_Z+1]]++) {

                                /* do linear interp at ends */

          sindex[ xyzv[VIO_Y] ] = 1;                             /* beginning end */
          index[ orig_xyzv[VIO_Y] ]=0;             
          GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
          index[ orig_xyzv[VIO_Y] ]=1;
          GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
          value1 = (v1+v2)/2;

          sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];

          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;

          sindex[ xyzv[VIO_Y] ] = 2*orig_count[ orig_xyzv[VIO_Y] ] -3;      /* ending end */
          index[ orig_xyzv[VIO_Y] ]= orig_count[ orig_xyzv[VIO_Y] ]-2;              
          GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
          index[ orig_xyzv[VIO_Y] ]= orig_count[ orig_xyzv[VIO_Y] ]-1;
          GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
          value1 = (v1+v2)/2;

          sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];

          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;

                                /* now do the voxels between the two ends */


          for(index[orig_xyzv[VIO_Y]]=1; index[orig_xyzv[VIO_Y]]<=orig_count[orig_xyzv[VIO_Y]]-3; index[orig_xyzv[VIO_Y]]++) {

            save_index = index[ orig_xyzv[VIO_Y] ];
            sindex[ xyzv[VIO_Y] ] = index[ orig_xyzv[VIO_Y] ]*2 + 1;

            index[ orig_xyzv[VIO_Y] ]--;             
            GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_Y] ]++;             
            GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_Y] ]++;             
            GET_VOXEL_4D(v3, orig_vol, index[0],index[1],index[2],index[3]);
            index[ orig_xyzv[VIO_Y] ]++;             
            GET_VOXEL_4D(v4, orig_vol, index[0],index[1],index[2],index[3]);

            value1 = MY_CUBIC_05(v1,v2,v3,v4);

            SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
            count++;

            index[ orig_xyzv[VIO_Y] ] = save_index;
          }

        }          
      update_progress_report( &progress, count+1 );
} 


    /* LEVEL 2: face interpolation, identified as 'f' in desc above */

               /* do faces in the plane with index[ orig_xyzv[VIO_X] ]=CONST */

    for(i=0; i<VIO_MAX_DIMENSIONS; i++)  sindex[i]=index[i]=0;
    


                                /* loop over the Z plane  */
    
    sindex[ xyzv[VIO_Z] ] = 0;

    for(index[orig_xyzv[VIO_Z+1]]=0; index[orig_xyzv[VIO_Z+1]]<orig_count[xyzv[VIO_Z+1]]; index[orig_xyzv[VIO_Z+1]]++) {
      
      sindex[ xyzv[VIO_Z+1] ] = index[ orig_xyzv[VIO_Z+1] ];
      
      
                                /* do faces near the edge first */
      
      for(index[orig_xyzv[VIO_Y]]=0; index[orig_xyzv[VIO_Y]]<orig_count[orig_xyzv[VIO_Y]]-1; index[orig_xyzv[VIO_Y]]++) {
        
        sindex[ xyzv[VIO_Y] ] = index[ orig_xyzv[VIO_Y] ]*2 + 1;
        
        sindex[ xyzv[VIO_X] ] = 0;
        GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
        sindex[ xyzv[VIO_X] ] = 2;
        GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
        
        value1 = (v1+v2)/2;;
          
        sindex[ xyzv[VIO_X] ] = 1;
        SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
        count++;
        
        sindex[ xyzv[VIO_X] ] = 2*orig_count[ orig_xyzv[VIO_X] ] -3;
        GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
        sindex[ xyzv[VIO_X] ] = 2*orig_count[ orig_xyzv[VIO_X] ] -1;
        GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
        
        value1 = (v1+v2)/2;
        
        sindex[ xyzv[VIO_X] ] = 2*orig_count[ orig_xyzv[VIO_X] ] -2;
        SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
        count++;
          
          
      }
      update_progress_report( &progress, count+1 );
      
      for(index[orig_xyzv[VIO_X]]=0; index[orig_xyzv[VIO_X]]<orig_count[orig_xyzv[VIO_X]]-1; index[orig_xyzv[VIO_X]]++) {
        
        sindex[ xyzv[VIO_X] ] = index[ orig_xyzv[VIO_X] ]*2 + 1;
        
        sindex[ xyzv[VIO_Y] ] = 0;
        GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
        sindex[ xyzv[VIO_Y] ] = 2;
        GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
        
        value1 = (v1+v2)/2;
        
        sindex[ xyzv[VIO_Y] ] = 1;
        SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
        count++;
          
        sindex[ xyzv[VIO_Y] ] = 2*orig_count[ orig_xyzv[VIO_Y] ] -3;
        GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
        sindex[ xyzv[VIO_Y] ] = 2*orig_count[ orig_xyzv[VIO_Y] ] -1;
        GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
        
        value1 = (v1+v2)/2;
        
        sindex[ xyzv[VIO_Y] ] = 2*orig_count[ orig_xyzv[VIO_Y] ] -2;
        SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
        count++;
          
          
      }
      update_progress_report( &progress, count+1 );
      
                              /* now do faces in the middle */
      
     
      for(index[orig_xyzv[VIO_X]]=1; index[orig_xyzv[VIO_X]]<orig_count[orig_xyzv[VIO_X]]-2; index[orig_xyzv[VIO_X]]++) {
        for(index[orig_xyzv[VIO_Y]]=1; index[orig_xyzv[VIO_Y]]<orig_count[orig_xyzv[VIO_Y]]-2; index[orig_xyzv[VIO_Y]]++) {
          
          sindex[ xyzv[VIO_X] ] = index[ orig_xyzv[VIO_X] ]*2 + 1;
          sindex[ xyzv[VIO_Y] ] = index[ orig_xyzv[VIO_Y] ]*2 + 1;
          
          save_index1 =  sindex[ xyzv[VIO_X] ];
          save_index2 =  sindex[ xyzv[VIO_Y] ];
          
          
          sindex[ xyzv[VIO_X] ] -= 3;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] += 2;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] += 2;
          GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] += 2;
          GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_X] ] = save_index1 ;
          value1 = MY_CUBIC_05(v1,v2,v3,v4);
          
          
          sindex[ xyzv[VIO_Y] ] -= 3;
          GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] += 2;
          GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] += 2;
          GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] += 2;
          GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
          sindex[ xyzv[VIO_Y] ] = save_index2 ;
          value2 = MY_CUBIC_05(v1,v2,v3,v4);
          
          
          value1 = (value1 + value2)/2.0;
          
          SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
          count++;
        }
        update_progress_report( &progress, count+1 );
      }

    } /* index[ orig_xyzv[VIO_Z+1] ] */
 
    terminate_progress_report( &progress );

}
