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
@MODIFIED   : Revision 96.5  2002-03-26 14:15:46  stever
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
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Optimize/super_sample_def.c,v 96.5 2002-03-26 14:15:46 stever Exp $";
#endif

#include <config.h>
#include <volume_io/internal_volume_io.h>
#include <print_error.h>
#include "constants.h"
#include "point_vector.h"

				/* prototypes called: */

public void get_volume_XYZV_indices(Volume data, int xyzv[]);
public void init_the_volume_to_zero(Volume volume);
public void interpolate_deformation_slice(Volume volume, 
					  Real wx,Real wy,Real wz,
					  Real def[]);
public void set_up_lattice(Volume data, 
			   double *user_step, 
			   double *start,     
			   int    *count,     
			   double *step,      
			   VectorR directions[]);


/* build the volume structure and allocate the data space to store
   a super-sampled GRID_TRANSFORM.

   *super_sampled must be ALLOCed before this call
   this routine will alloc the volume data space.

   super_step specifies the number of times to super-sample the data.
 */
public void create_super_sampled_data_volumes(General_transform *orig_deformation,
					      General_transform *super_sampled,
					      int super_step)

{

  int 
    i,
    xyzv[MAX_DIMENSIONS],
    orig_count[MAX_DIMENSIONS], 
    xyz_count[MAX_DIMENSIONS], 
    new_count[MAX_DIMENSIONS];
  Real
    voxel[MAX_DIMENSIONS],
    start[MAX_DIMENSIONS],
    orig_steps[MAX_DIMENSIONS],
    new_steps[MAX_DIMENSIONS],
    xyz_steps[MAX_DIMENSIONS];

  VectorR 
    directions[3];

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

				/* prepare to modify the GRID_TRANSFORM */

  get_volume_XYZV_indices(orig_deformation->displacement_volume, xyzv);
  get_volume_sizes(       orig_deformation->displacement_volume, orig_count);
  get_volume_separations( orig_deformation->displacement_volume, orig_steps);
  
  for_less(i,0,3) {		/* set_up_lattice needs XYZ order.... */
    new_steps[i] = orig_steps[xyzv[i]]/super_step;
  }
  
  set_up_lattice(orig_deformation->displacement_volume, 
		 new_steps, start, xyz_count, xyz_steps, directions);

      /* use the sign of the step returned to set the true step size */
  for_less(i,0,3)	
    if (xyz_steps[i]<0) 
      xyz_steps[i] = -1.0 * ABS(new_steps [i]); 
    else 
      xyz_steps[i] = ABS(new_steps[i]);
  
  for_less(i,0,MAX_DIMENSIONS) {
    voxel[i] = 0.0;
    new_count[i] = orig_count[i];
    new_steps[i] = orig_steps[i];
  }
  for_less(i,0,3) {
    new_count[ xyzv[i] ] = xyz_count[i];
    new_steps[ xyzv[i] ] = xyz_steps[i];
  }

  set_volume_sizes(       super_sampled->displacement_volume, new_count);
  set_volume_separations( super_sampled->displacement_volume, new_steps);
  set_volume_translation( super_sampled->displacement_volume, voxel, start);
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
*/
public void interpolate_super_sampled_data(General_transform *orig_deformation,
					   General_transform *super_sampled,
					   int dim)
{
  Volume
    orig_vol,
    super_vol;

  int 
    i,
    index[MAX_DIMENSIONS],
    orig_xyzv[MAX_DIMENSIONS],
    xyzv[MAX_DIMENSIONS],
    count[MAX_DIMENSIONS];
  Real
    def_vector[N_DIMENSIONS],
    voxel[MAX_DIMENSIONS],
    wx,wy,wz;
  progress_struct
    progress;

  orig_vol = orig_deformation->displacement_volume;
  super_vol= super_sampled->displacement_volume;

  get_volume_sizes(super_vol, count);
  get_volume_XYZV_indices(super_vol, xyzv);
  get_volume_XYZV_indices(orig_vol, orig_xyzv);
  
  initialize_progress_report( &progress, FALSE, count[ xyzv[X] ],
			     "Super-sampling defs:" );
  
  for_less(i,0,MAX_DIMENSIONS) index[i]=0;
  
  for_less( index[ xyzv[X] ] , 0, count[ xyzv[X] ]) {
    for_less( index[ xyzv[Y] ] , 0, count[ xyzv[Y] ]) {
      for_less( index[ xyzv[Z] ] , 0, count[ xyzv[Z] ]) {
	
	for_less(i,0,MAX_DIMENSIONS) 
	  voxel[i]=(Real)index[i];
	convert_voxel_to_world(super_vol, 
			       voxel,
			       &wx, &wy, &wz);
	if (dim==2) 
	  interpolate_deformation_slice(orig_vol, wx, wy, wz, def_vector);
	else {
	  evaluate_volume_in_world(orig_vol, wx,wy,wz, 2, TRUE, 0.0,
				   def_vector,
				   NULL, NULL, NULL,
				   NULL, NULL, NULL, NULL, NULL, NULL);
	}

	for_less( index[ xyzv[Z+1] ], 0, count[ xyzv[Z+1] ]) 
	  set_volume_real_value(super_vol,
				index[0],index[1],index[2],index[3],index[4],
				def_vector[ index[ xyzv[Z+1] ] ]);
	
      }
    }
    update_progress_report( &progress, index[ xyzv[X] ]+1 );
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

public void create_super_sampled_data_volumes_by2(General_transform *orig_deformation,
						  General_transform *super_sampled)

{

  int 
    i,
    xyzv[MAX_DIMENSIONS],
    new_xyzv[MAX_DIMENSIONS], 
    orig_count[MAX_DIMENSIONS], 
    new_count[MAX_DIMENSIONS];
  Real
    voxel[MAX_DIMENSIONS],
    start[MAX_DIMENSIONS],
    orig_steps[MAX_DIMENSIONS],
    new_steps[MAX_DIMENSIONS];

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

				/* prepare to modify the GRID_TRANSFORM */

  get_volume_XYZV_indices(orig_deformation->displacement_volume, xyzv);
  get_volume_XYZV_indices(super_sampled->displacement_volume,    new_xyzv);
  get_volume_sizes(       orig_deformation->displacement_volume, orig_count);
  get_volume_separations( orig_deformation->displacement_volume, orig_steps);
  

  for_less(i,0, get_volume_n_dimensions(orig_deformation->displacement_volume)) {		
    new_steps[new_xyzv[i]] = orig_steps[xyzv[i]];
    new_count[new_xyzv[i]] = orig_count[xyzv[i]];
  }
				/* set up the new step size */
  for_less(i,0,3) {		
    new_steps[new_xyzv[i]] = orig_steps[xyzv[i]]/2.0;
  }
  
				/* set up the new number of elements size */
  for_less(i,0,3) {		
    new_count[new_xyzv[i]] = orig_count[xyzv[i]]*2.0 - 1;
  }
  
  for_less(i,0,MAX_DIMENSIONS) 
    voxel[i] = 0.0;

  convert_voxel_to_world(orig_deformation->displacement_volume,
			 voxel,
			 &start[X], &start[Y], &start[Z]);
  
  
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

private void interpolate_super_sampled_data_by2_dim3( 
    General_transform *orig_deformation,
    General_transform *super_sampled )
{
  Volume
    orig_vol,
    super_vol;

  int 
    i, save_index, save_index1, save_index2,
    count,
    index[MAX_DIMENSIONS],
    sindex[MAX_DIMENSIONS],
    orig_xyzv[MAX_DIMENSIONS],
    xyzv[MAX_DIMENSIONS],
    orig_count[MAX_DIMENSIONS],
    super_count[MAX_DIMENSIONS];
  Real
    value1, value2, value3, v1,v2,v3,v4,v5,v6;
  progress_struct
    progress;

    orig_vol = orig_deformation->displacement_volume;
    super_vol= super_sampled->displacement_volume;
    
    get_volume_sizes(       super_vol, super_count);
    get_volume_sizes(       orig_vol,  orig_count);
    get_volume_XYZV_indices(super_vol, xyzv);
    get_volume_XYZV_indices(orig_vol,  orig_xyzv);
    
    init_the_volume_to_zero(super_sampled->displacement_volume);

    initialize_progress_report(&progress, FALSE, 
			       super_count[ xyzv[X] ]*super_count[ xyzv[Y] ]*
			       super_count[ xyzv[Z] ]*super_count[ xyzv[Z+1] ]+1,
			       "Super-sampling defs:" );
    count = 0;



    /* LEVEL 0: copy original 'corner' nodes, identified as 'X' in desc above */

    for_less(i,0,MAX_DIMENSIONS) sindex[i]=index[i]=0;
    
    for_less( index[ orig_xyzv[X] ] , 0, orig_count[ orig_xyzv[X] ]) {
      for_less( index[ orig_xyzv[Y] ] , 0, orig_count[ orig_xyzv[Y] ]) {
	for_less( index[ orig_xyzv[Z] ] , 0, orig_count[ orig_xyzv[Z] ]) {

	  for_less(i,0,N_DIMENSIONS)
	    sindex[ xyzv[i] ] = 2*index[ orig_xyzv[i] ];

	  for_less( index[ orig_xyzv[Z+1] ], 0, orig_count[ xyzv[Z+1] ]) {
	    sindex[ xyzv[Z+1] ] = index[ orig_xyzv[Z+1] ];
	    GET_VOXEL_4D(value1, orig_vol, index[0],index[1],index[2],index[3]);
	    SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  }	  

	  count += 3;
	}
      }
      update_progress_report( &progress, count+1 );
    }


    /* LEVEL 1: edge interpolation, identified as 'e' in desc above */

               /* do edges along the index[ orig_xyzv[X] ] dir */

    for_less(i,0,MAX_DIMENSIONS) sindex[i]=index[i]=0;
    
				/* loop over all x-dirs  */

    for_less( index[ orig_xyzv[Y] ], 0, orig_count[ orig_xyzv[Y] ]) {


      for_less( index[ orig_xyzv[Z] ] , 0, orig_count[ orig_xyzv[Z] ]) {

	sindex[ xyzv[Y] ] = 2*index[ orig_xyzv[Y] ];
	sindex[ xyzv[Z] ] = 2*index[ orig_xyzv[Z] ];

	for_less( index[ orig_xyzv[Z+1] ], 0, orig_count[ xyzv[Z+1] ]) {

				/* do linear interp at ends */

	  sindex[ xyzv[X] ] = 1;                             /* beginning end */
	  index[ orig_xyzv[X] ]=0;             
	  GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
	  index[ orig_xyzv[X] ]=1;
	  GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
	  value1 = (v1+v2)/2;

	  sindex[ xyzv[Z+1] ] = index[ orig_xyzv[Z+1] ];

	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;

	  sindex[ xyzv[X] ] = super_count[ xyzv[X] ] -2;      /* ending end */
	  index[ orig_xyzv[X] ]= orig_count[ orig_xyzv[X] ]-2;              
	  GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
	  index[ orig_xyzv[X] ]= orig_count[ orig_xyzv[X] ]-1;
	  GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
	  value1 = (v1+v2)/2;

	  sindex[ xyzv[Z+1] ] = index[ orig_xyzv[Z+1] ];

	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;

				/* now do the voxels between the two ends */


	  for_inclusive( index[ orig_xyzv[X] ], 1, orig_count[ orig_xyzv[X] ]-3) {

	    save_index = index[ orig_xyzv[X] ];
	    sindex[ xyzv[X] ] = index[ orig_xyzv[X] ]*2 + 1;

	    index[ orig_xyzv[X] ]--;             
	    GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
	    index[ orig_xyzv[X] ]++;             
	    GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
	    index[ orig_xyzv[X] ]++;             
	    GET_VOXEL_4D(v3, orig_vol, index[0],index[1],index[2],index[3]);
	    index[ orig_xyzv[X] ]++;             
	    GET_VOXEL_4D(v4, orig_vol, index[0],index[1],index[2],index[3]);

	    value1 = MY_CUBIC_05(v1,v2,v3,v4);

	    SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	    count++;

	    index[ orig_xyzv[X] ] = save_index;
	  }

	}	  

      }
      update_progress_report( &progress, count+1 );
    }
               /* do edges along the index[ orig_xyzv[Y] ] dir */

    for_less(i,0,MAX_DIMENSIONS) sindex[i]=index[i]=0;
    
				/* loop over all x-dirs  */

    for_less( index[ orig_xyzv[Z] ], 0, orig_count[ orig_xyzv[Z] ]) {


      for_less( index[ orig_xyzv[X] ] , 0, orig_count[ orig_xyzv[X] ]) {


	sindex[ xyzv[Z] ] = 2*index[ orig_xyzv[Z] ];
	sindex[ xyzv[X] ] = 2*index[ orig_xyzv[X] ];

	for_less( index[ orig_xyzv[Z+1] ], 0, orig_count[ xyzv[Z+1] ]) {

				/* do linear interp at ends */

	  sindex[ xyzv[Y] ] = 1;                             /* beginning end */
	  index[ orig_xyzv[Y] ]=0;             
	  GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
	  index[ orig_xyzv[Y] ]=1;
	  GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
	  value1 = (v1+v2)/2;

	  sindex[ xyzv[Z+1] ] = index[ orig_xyzv[Z+1] ];

	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;

	  sindex[ xyzv[Y] ] = super_count[ xyzv[Y] ] -2;      /* ending end */
	  index[ orig_xyzv[Y] ]= orig_count[ orig_xyzv[Y] ]-2;              
	  GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
	  index[ orig_xyzv[Y] ]= orig_count[ orig_xyzv[Y] ]-1;
	  GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
	  value1 = (v1+v2)/2;

	  sindex[ xyzv[Z+1] ] = index[ orig_xyzv[Z+1] ];

	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;

				/* now do the voxels between the two ends */


	  for_inclusive( index[ orig_xyzv[Y] ], 1, orig_count[ orig_xyzv[Y] ]-3) {

	    save_index = index[ orig_xyzv[Y] ];
	    sindex[ xyzv[Y] ] = index[ orig_xyzv[Y] ]*2 + 1;

	    index[ orig_xyzv[Y] ]--;             
	    GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
	    index[ orig_xyzv[Y] ]++;             
	    GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
	    index[ orig_xyzv[Y] ]++;             
	    GET_VOXEL_4D(v3, orig_vol, index[0],index[1],index[2],index[3]);
	    index[ orig_xyzv[Y] ]++;             
	    GET_VOXEL_4D(v4, orig_vol, index[0],index[1],index[2],index[3]);

	    value1 = MY_CUBIC_05(v1,v2,v3,v4);

	    SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	    count++;

	    index[ orig_xyzv[Y] ] = save_index;
	  }

	}	  

      }
      update_progress_report( &progress, count+1 );
    }

               /* do edges along the index[ orig_xyzv[Z] ] dir */

    for_less(i,0,MAX_DIMENSIONS) sindex[i]=index[i]=0;
    
				/* loop over all Z-dirs  */

    for_less( index[ orig_xyzv[X] ], 0, orig_count[ orig_xyzv[X] ]) {


      for_less( index[ orig_xyzv[Y] ] , 0, orig_count[ orig_xyzv[Y] ]) {

	sindex[ xyzv[X] ] = 2*index[ orig_xyzv[X] ];
	sindex[ xyzv[Y] ] = 2*index[ orig_xyzv[Y] ];

	for_less( index[ orig_xyzv[Z+1] ], 0, orig_count[ xyzv[Z+1] ]) {

				/* do linear interp at ends */

	  sindex[ xyzv[Z] ] = 1;                             /* beginning end */
	  index[ orig_xyzv[Z] ]=0;             
	  GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
	  index[ orig_xyzv[Z] ]=1;
	  GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
	  value1 = (v1+v2)/2;

	  sindex[ xyzv[Z+1] ] = index[ orig_xyzv[Z+1] ];

	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;

	  sindex[ xyzv[Z] ] = super_count[ xyzv[Z] ] -2;      /* ending end */
	  index[ orig_xyzv[Z] ]= orig_count[ orig_xyzv[Z] ]-2;              
	  GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
	  index[ orig_xyzv[Z] ]= orig_count[ orig_xyzv[Z] ]-1;
	  GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
	  value1 = (v1+v2)/2;

	  sindex[ xyzv[Z+1] ] = index[ orig_xyzv[Z+1] ];

	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;

				/* now do the voxels between the two ends */


	  for_inclusive( index[ orig_xyzv[Z] ], 1, orig_count[ orig_xyzv[Z] ]-3) {

	    save_index = index[ orig_xyzv[Z] ];
	    sindex[ xyzv[Z] ] = index[ orig_xyzv[Z] ]*2 + 1;

	    index[ orig_xyzv[Z] ]--;             
	    GET_VOXEL_4D(v1, orig_vol, index[0],index[1],index[2],index[3]);
	    index[ orig_xyzv[Z] ]++;             
	    GET_VOXEL_4D(v2, orig_vol, index[0],index[1],index[2],index[3]);
	    index[ orig_xyzv[Z] ]++;             
	    GET_VOXEL_4D(v3, orig_vol, index[0],index[1],index[2],index[3]);
	    index[ orig_xyzv[Z] ]++;             
	    GET_VOXEL_4D(v4, orig_vol, index[0],index[1],index[2],index[3]);

	    value1 = MY_CUBIC_05(v1,v2,v3,v4);

	    SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	    count++;

	    index[ orig_xyzv[Z] ] = save_index;
	  }

	}	  

      }
      update_progress_report( &progress, count+1 );
    }


    /* LEVEL 2: face interpolation, identified as 'f' in desc above */

               /* do faces in the plane with index[ orig_xyzv[X] ]=CONST */

    for_less(i,0,MAX_DIMENSIONS)  sindex[i]=index[i]=0;
    
				/* loop over all X planes  */

    for_less( index[ orig_xyzv[X] ], 0, orig_count[ orig_xyzv[X] ]-1) {

      sindex[ xyzv[X] ] = 2*index[ orig_xyzv[X] ];

      for_less( index[ orig_xyzv[Z+1] ], 0, orig_count[ xyzv[Z+1] ]) {

	sindex[ xyzv[Z+1] ] = index[ orig_xyzv[Z+1] ];


				/* do faces near the edge first */

	for_less( index[ orig_xyzv[Z] ] , 0, orig_count[ orig_xyzv[Z] ]-1) {

	  sindex[ xyzv[Z] ] = index[ orig_xyzv[Z] ]*2 + 1;
	  
	  sindex[ xyzv[Y] ] = 0;
	  GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Y] ] = 2;
	  GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  
	  value1 = (v1+v2)/2;
	  
	  sindex[ xyzv[Y] ] = 1;
	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;
	  
	  sindex[ xyzv[Y] ] = super_count[ xyzv[Y] ] -3;
	  GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Y] ] = super_count[ xyzv[Y] ] -1;
	  GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  
	  value1 = (v1+v2)/2;
	  
	  sindex[ xyzv[Y] ] = super_count[ xyzv[Y] ] -2;
	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;
	  
	  
	}
	update_progress_report( &progress, count+1 );

	for_less( index[ orig_xyzv[Y] ] , 0, orig_count[ orig_xyzv[Y] ]-1) {

	  sindex[ xyzv[Y] ] = index[ orig_xyzv[Y] ]*2 + 1;
	  
	  sindex[ xyzv[Z] ] = 0;
	  GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Z] ] = 2;
	  GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  
	  value1 = (v1+v2)/2;
	  
	  sindex[ xyzv[Z] ] = 1;
	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;
	  
	  sindex[ xyzv[Z] ] = super_count[ xyzv[Z] ] -3;
	  GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Z] ] = super_count[ xyzv[Z] ] -1;
	  GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  
	  value1 = (v1+v2)/2;
	  
	  sindex[ xyzv[Z] ] = super_count[ xyzv[Z] ] -2;
	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;
	  
	  
	}
	update_progress_report( &progress, count+1 );
      
				/* now do faces in the middle */


	for_less( index[ orig_xyzv[Y] ] , 1, orig_count[ orig_xyzv[Y] ]-2) {
	  for_less( index[ orig_xyzv[Z] ] , 1, orig_count[ orig_xyzv[Z] ]-2) {

	    sindex[ xyzv[Y] ] = index[ orig_xyzv[Y] ]*2 + 1;
	    sindex[ xyzv[Z] ] = index[ orig_xyzv[Z] ]*2 + 1;

	    save_index1 =  sindex[ xyzv[Y] ];
	    save_index2 =  sindex[ xyzv[Z] ];


	    sindex[ xyzv[Y] ] -= 3;
	    GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Y] ] += 2;
	    GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Y] ] += 2;
	    GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Y] ] += 2;
	    GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Y] ] = save_index1 ;
	    value1 = MY_CUBIC_05(v1,v2,v3,v4);


	    sindex[ xyzv[Z] ] -= 3;
	    GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Z] ] += 2;
	    GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Z] ] += 2;
	    GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Z] ] += 2;
	    GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Z] ] = save_index2 ;
	    value2 = MY_CUBIC_05(v1,v2,v3,v4);


	    value1 = (value1 + value2)/2.0;

	    SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	    count++;
	  }
	  update_progress_report( &progress, count+1 );
	}

      } /* index[ orig_xyzv[Z+1] ] */
      
    }
				/* loop over all Y planes  */

    for_less( index[ orig_xyzv[Y] ], 0, orig_count[ orig_xyzv[Y] ]-1) {

      sindex[ xyzv[Y] ] = 2*index[ orig_xyzv[Y] ];

      for_less( index[ orig_xyzv[Z+1] ], 0, orig_count[ xyzv[Z+1] ]) {

	sindex[ xyzv[Z+1] ] = index[ orig_xyzv[Z+1] ];


				/* do faces near the edge first */

	for_less( index[ orig_xyzv[X] ] , 0, orig_count[ orig_xyzv[X] ]-1) {

	  sindex[ xyzv[X] ] = index[ orig_xyzv[X] ]*2 + 1;
	  
	  sindex[ xyzv[Z] ] = 0;
	  GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Z] ] = 2;
	  GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  
	  value1 = (v1+v2)/2;
	  
	  sindex[ xyzv[Z] ] = 1;
	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;
	  
	  sindex[ xyzv[Z] ] = super_count[ xyzv[Z] ] -3;
	  GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Z] ] = super_count[ xyzv[Z] ] -1;
	  GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  
	  value1 = (v1+v2)/2;
	  
	  sindex[ xyzv[Z] ] = super_count[ xyzv[Z] ] -2;
	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;
	  
	  
	}
	update_progress_report( &progress, count+1 );

	for_less( index[ orig_xyzv[Z] ] , 0, orig_count[ orig_xyzv[Z] ]-1) {

	  sindex[ xyzv[Z] ] = index[ orig_xyzv[Z] ]*2 + 1;
	  
	  sindex[ xyzv[X] ] = 0;
	  GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[X] ] = 2;
	  GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  
	  value1 = (v1+v2)/2;
	  
	  sindex[ xyzv[X] ] = 1;
	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;
	  
	  sindex[ xyzv[X] ] = super_count[ xyzv[X] ] -3;
	  GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[X] ] = super_count[ xyzv[X] ] -1;
	  GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  
	  value1 = (v1+v2)/2;
	  
	  sindex[ xyzv[X] ] = super_count[ xyzv[X] ] -2;
	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;
	  
	  
	}
	update_progress_report( &progress, count+1 );

				/* now do faces in the middle */


	for_less( index[ orig_xyzv[Z] ] , 1, orig_count[ orig_xyzv[Z] ]-2) {
	  for_less( index[ orig_xyzv[X] ] , 1, orig_count[ orig_xyzv[X] ]-2) {

	    sindex[ xyzv[Z] ] = index[ orig_xyzv[Z] ]*2 + 1;
	    sindex[ xyzv[X] ] = index[ orig_xyzv[X] ]*2 + 1;

	    save_index1 =  sindex[ xyzv[Z] ];
	    save_index2 =  sindex[ xyzv[X] ];


	    sindex[ xyzv[Z] ] -= 3;
	    GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Z] ] += 2;
	    GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Z] ] += 2;
	    GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Z] ] += 2;
	    GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Z] ] = save_index1 ;
	    value1 = MY_CUBIC_05(v1,v2,v3,v4);


	    sindex[ xyzv[X] ] -= 3;
	    GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[X] ] += 2;
	    GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[X] ] += 2;
	    GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[X] ] += 2;
	    GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[X] ] = save_index2 ;
	    value2 = MY_CUBIC_05(v1,v2,v3,v4);


	    value1 = (value1 + value2)/2.0;

	    SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	    count++;
	  }
	  update_progress_report( &progress, count+1 );
	}

      } /* index[ orig_xyzv[Z+1] ] */
      
    }

				/* loop over all Z planes  */
    
    for_less( index[ orig_xyzv[Z] ], 0, orig_count[ orig_xyzv[Z] ]-1) {

      sindex[ xyzv[Z] ] = 2*index[ orig_xyzv[Z] ];

      for_less( index[ orig_xyzv[Z+1] ], 0, orig_count[ xyzv[Z+1] ]) {

	sindex[ xyzv[Z+1] ] = index[ orig_xyzv[Z+1] ];


				/* do faces near the edge first */

	for_less( index[ orig_xyzv[Y] ] , 0, orig_count[ orig_xyzv[Y] ]-1) {

	  sindex[ xyzv[Y] ] = index[ orig_xyzv[Y] ]*2 + 1;
	  
	  sindex[ xyzv[X] ] = 0;
	  GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[X] ] = 2;
	  GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  
	  value1 = (v1+v2)/2;
	  
	  sindex[ xyzv[X] ] = 1;
	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;
	  
	  sindex[ xyzv[X] ] = super_count[ xyzv[X] ] -3;
	  GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[X] ] = super_count[ xyzv[X] ] -1;
	  GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  
	  value1 = (v1+v2)/2;
	  
	  sindex[ xyzv[X] ] = super_count[ xyzv[X] ] -2;
	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;
	  
	  
	}
	update_progress_report( &progress, count+1 );

	for_less( index[ orig_xyzv[X] ] , 0, orig_count[ orig_xyzv[X] ]-1) {

	  sindex[ xyzv[X] ] = index[ orig_xyzv[X] ]*2 + 1;
	  
	  sindex[ xyzv[Y] ] = 0;
	  GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Y] ] = 2;
	  GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  
	  value1 = (v1+v2)/2;
	  
	  sindex[ xyzv[Y] ] = 1;
	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;
	  
	  sindex[ xyzv[Y] ] = super_count[ xyzv[Y] ] -3;
	  GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Y] ] = super_count[ xyzv[Y] ] -1;
	  GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  
	  value1 = (v1+v2)/2;
	  
	  sindex[ xyzv[Y] ] = super_count[ xyzv[Y] ] -2;
	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;
	  
	  
	}
	update_progress_report( &progress, count+1 );

				/* now do faces in the middle */


	for_less( index[ orig_xyzv[X] ] , 1, orig_count[ orig_xyzv[X] ]-2) {
	  for_less( index[ orig_xyzv[Y] ] , 1, orig_count[ orig_xyzv[Y] ]-2) {

	    sindex[ xyzv[X] ] = index[ orig_xyzv[X] ]*2 + 1;
	    sindex[ xyzv[Y] ] = index[ orig_xyzv[Y] ]*2 + 1;

	    save_index1 =  sindex[ xyzv[X] ];
	    save_index2 =  sindex[ xyzv[Y] ];


	    sindex[ xyzv[X] ] -= 3;
	    GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[X] ] += 2;
	    GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[X] ] += 2;
	    GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[X] ] += 2;
	    GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[X] ] = save_index1 ;
	    value1 = MY_CUBIC_05(v1,v2,v3,v4);


	    sindex[ xyzv[Y] ] -= 3;
	    GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Y] ] += 2;
	    GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Y] ] += 2;
	    GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Y] ] += 2;
	    GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Y] ] = save_index2 ;
	    value2 = MY_CUBIC_05(v1,v2,v3,v4);


	    value1 = (value1 + value2)/2.0;

	    SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	    count++;
	  }
	  update_progress_report( &progress, count+1 );
	}

      } /* index[ orig_xyzv[Z+1] ] */
      
    }
    /* LEVEL 3: center interpolation, identified as 'c' in desc above */

    for_less( index[ orig_xyzv[Z+1] ], 0, orig_count[ xyzv[Z+1] ]) {

      sindex[ xyzv[Z+1] ] = index[ orig_xyzv[Z+1] ];

                         /* do all centers that are 1 sample from
			    the edge of the super-sampled volume */

				/* do the two X-planes */

      for_less( index[ orig_xyzv[Y] ] , 0, orig_count[ orig_xyzv[Y] ]-2) {
	for_less( index[ orig_xyzv[Z] ] , 0, orig_count[ orig_xyzv[Z] ]-2) {

	  sindex[ xyzv[Y] ] = index[ orig_xyzv[Y] ]*2 + 1;
	  sindex[ xyzv[Z] ] = index[ orig_xyzv[Z] ]*2 + 1;

	  sindex[ xyzv[X] ] = 1;

	  sindex[ xyzv[X] ] -= 1;
	  GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[X] ] += 2;
	  GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[X] ] -= 1;
	  
	  sindex[ xyzv[Y] ] -= 1;
	  GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Y] ] += 2;
	  GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Y] ] -= 1;
	
	  sindex[ xyzv[Z] ] -= 1;
	  GET_VOXEL_4D(v5, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Z] ] += 2;
	  GET_VOXEL_4D(v6, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Z] ] -= 1;
	
	  value1 = (v1 + v2 + v3 + v4 + v5 + v6) / 6.0;
	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;
	
	  sindex[ xyzv[X] ] = super_count[ xyzv[X] ] -2;
      
	  sindex[ xyzv[X] ] -= 1;
	  GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[X] ] += 2;
	  GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[X] ] -= 1;
	  
	  sindex[ xyzv[Y] ] -= 1;
	  GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Y] ] += 2;
	  GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Y] ] -= 1;
	  
	  sindex[ xyzv[Z] ] -= 1;
	  GET_VOXEL_4D(v5, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Z] ] += 2;
	  GET_VOXEL_4D(v6, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Z] ] -= 1;
	  
	  value1 = (v1 + v2 + v3 + v4 + v5 + v6) / 6.0;
	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;
	}
      }
      update_progress_report( &progress, count+1 );
				/* do the two Y-planes */

      for_less( index[ orig_xyzv[X] ] , 0, orig_count[ orig_xyzv[X] ]-2) {
	for_less( index[ orig_xyzv[Z] ] , 0, orig_count[ orig_xyzv[Z] ]-2) {

	  sindex[ xyzv[X] ] = index[ orig_xyzv[X] ]*2 + 1;
	  sindex[ xyzv[Z] ] = index[ orig_xyzv[Z] ]*2 + 1;

	  sindex[ xyzv[Y] ] = 1;

	  sindex[ xyzv[X] ] -= 1;
	  GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[X] ] += 2;
	  GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[X] ] -= 1;
	  
	  sindex[ xyzv[Y] ] -= 1;
	  GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Y] ] += 2;
	  GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Y] ] -= 1;
	
	  sindex[ xyzv[Z] ] -= 1;
	  GET_VOXEL_4D(v5, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Z] ] += 2;
	  GET_VOXEL_4D(v6, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Z] ] -= 1;
	
	  value1 = (v1 + v2 + v3 + v4 + v5 + v6) / 6.0;
	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;
	
	  sindex[ xyzv[Y] ] = super_count[ xyzv[Y] ] -2;
      
	  sindex[ xyzv[X] ] -= 1;
	  GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[X] ] += 2;
	  GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[X] ] -= 1;
	  
	  sindex[ xyzv[Y] ] -= 1;
	  GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Y] ] += 2;
	  GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Y] ] -= 1;
	  
	  sindex[ xyzv[Z] ] -= 1;
	  GET_VOXEL_4D(v5, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Z] ] += 2;
	  GET_VOXEL_4D(v6, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Z] ] -= 1;
	  
	  value1 = (v1 + v2 + v3 + v4 + v5 + v6) / 6.0;
	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;
	}
      }
      update_progress_report( &progress, count+1 );

				/* do the two Z-planes */

      for_less( index[ orig_xyzv[X] ] , 0, orig_count[ orig_xyzv[X] ]-2) {
	for_less( index[ orig_xyzv[Y] ] , 0, orig_count[ orig_xyzv[Y] ]-2) {

	  sindex[ xyzv[X] ] = index[ orig_xyzv[X] ]*2 + 1;
	  sindex[ xyzv[Y] ] = index[ orig_xyzv[Y] ]*2 + 1;

	  sindex[ xyzv[Z] ] = 1;

	  sindex[ xyzv[X] ] -= 1;
	  GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[X] ] += 2;
	  GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[X] ] -= 1;
	  
	  sindex[ xyzv[Y] ] -= 1;
	  GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Y] ] += 2;
	  GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Y] ] -= 1;
	
	  sindex[ xyzv[Z] ] -= 1;
	  GET_VOXEL_4D(v5, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Z] ] += 2;
	  GET_VOXEL_4D(v6, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Z] ] -= 1;
	
	  value1 = (v1 + v2 + v3 + v4 + v5 + v6) / 6.0;
	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;
	
	  sindex[ xyzv[Z] ] = super_count[ xyzv[Z] ] -2;
      
	  sindex[ xyzv[X] ] -= 1;
	  GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[X] ] += 2;
	  GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[X] ] -= 1;
	  
	  sindex[ xyzv[Y] ] -= 1;
	  GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Y] ] += 2;
	  GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Y] ] -= 1;
	  
	  sindex[ xyzv[Z] ] -= 1;
	  GET_VOXEL_4D(v5, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Z] ] += 2;
	  GET_VOXEL_4D(v6, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	  sindex[ xyzv[Z] ] -= 1;
	  
	  value1 = (v1 + v2 + v3 + v4 + v5 + v6) / 6.0;
	  SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	  count++;
	}
      }
      update_progress_report( &progress, count+1 );

                         /* now do all central voxels */

      for_less( index[ orig_xyzv[X] ] , 1, orig_count[ orig_xyzv[X] ]-2) {
	for_less( index[ orig_xyzv[Y] ] , 1, orig_count[ orig_xyzv[Y] ]-2) {
	  for_less( index[ orig_xyzv[Z] ] , 1, orig_count[ orig_xyzv[Z] ]-2) {

	    sindex[ xyzv[X] ] = index[ orig_xyzv[X] ]*2 + 1;
	    sindex[ xyzv[Y] ] = index[ orig_xyzv[Y] ]*2 + 1;
	    sindex[ xyzv[Z] ] = index[ orig_xyzv[Z] ]*2 + 1;

	    sindex[ xyzv[X] ] -= 3;
	    GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[X] ] += 2;
	    GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[X] ] += 2;
	    GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[X] ] += 2;
	    GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[X] ] -= 3;
	    value1 = MY_CUBIC_05(v1,v2,v3,v4);

	    sindex[ xyzv[Y] ] -= 3;
	    GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Y] ] += 2;
	    GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Y] ] += 2;
	    GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Y] ] += 2;
	    GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Y] ] -= 3;
	    value2 = MY_CUBIC_05(v1,v2,v3,v4);

	    sindex[ xyzv[Z] ] -= 3;
	    GET_VOXEL_4D(v1, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Z] ] += 2;
	    GET_VOXEL_4D(v2, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Z] ] += 2;
	    GET_VOXEL_4D(v3, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Z] ] += 2;
	    GET_VOXEL_4D(v4, super_vol, sindex[0],sindex[1],sindex[2],sindex[3]);
	    sindex[ xyzv[Z] ] -= 3;
	    value3 = MY_CUBIC_05(v1,v2,v3,v4);

	    
	    value1 = (value1 + value2 + value3) / 3.0;
	    SET_VOXEL_4D(super_vol, sindex[0],sindex[1],sindex[2],sindex[3], value1);
	    count++;

	  }
	}
	update_progress_report( &progress, count+1 );
      }
    } /* index[ orig_xyzv[Z+1] ] */

    terminate_progress_report( &progress );

}


public void interpolate_super_sampled_data_by2(
    General_transform *orig_deformation,
    General_transform *super_sampled,
    int dim)
{
    if (orig_deformation->type != GRID_TRANSFORM || super_sampled->type != GRID_TRANSFORM) {
	print_error_and_line_num("interpolate_super_sampled_data_by2 not called with GRID_TRANSFORM",
				 __FILE__, __LINE__);
    }

    if (dim == 2) {
	interpolate_super_sampled_data(orig_deformation,super_sampled,dim);
    } else {
	interpolate_super_sampled_data_by2_dim3( orig_deformation, super_sampled );
    }
}


