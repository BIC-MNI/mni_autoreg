/* ----------------------------- MNI Header -----------------------------------
@NAME       : default_def.c
@DESCRIPTION: 
   build_default_deformation_field: using the info in *globals, this
   routine will build the structure that will represent the deformation
   field that is to be optimized in the program.

   there are two possibilities:

1- If there is NO non-linear deformation (as defined by a GRID_TRANSFORM) 
   in the globals->trans_info.transformation, then the routine will:

     concatenate a zero-valued non-linear deformation field transformation 
     to the end of the globals->trans_info.transformation, using the lattice 
     information (in *globals) to build the deformation volumes.

2- If there is already a deformation field in the input transformation, then 
   this deformation will be replaced by its equivalent, i.e. the field
   must be replaced by a structure that corresponds to the data in *globals,
   and the values in this new structure must be interpolated from the
   existing deformation field. 

@COPYRIGHT  :
              Copyright 1993 Louis Collins, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.

@CREATED    : Thu Apr 27 10:00:52 MET DST 1995
      created by removing build_default_deformation_field from 
      transformations.c
@MODIFIED   : $Log: default_def.c,v $
@MODIFIED   : Revision 9.5  1996-08-12 14:15:50  louis
@MODIFIED   : Release of MNI_AutoReg version 1.0
@MODIFIED   :
 * Revision 1.4  1996/08/12  14:15:44  louis
 * Pre-release
 *
 * Revision 1.3  1995/09/11  12:49:12  collins
 * removed reference to get_linear_part_of_transformation(), which was no
 * longer needed or used.
 *
 * Revision 1.2  1995/09/11  12:37:16  collins
 * updated working version - corresponds to mni_reg-0.1g
 *
 * Revision 1.1  1995/05/02  11:31:53  collins
 * Initial revision
 *
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="";
#endif

#include <config.h>		/* to have HAVE_RECENT_VOLUME_IO macro */
#include <volume_io.h>
#include <print_error.h>
#include "arg_data.h"
#include "local_macros.h"
#include "constants.h"

#define MY_MAX_VOX 32766.0
#define MY_MAX_REAL 50.0

static  char *dim_name_vector_vol[] = 
                { MIvector_dimension, MIzspace,MIyspace, MIxspace };

public void set_up_lattice(Volume data,       /* in: volume  */
			   double *user_step, /* in: user requested spacing for lattice */
			   double *start,     /* out:world starting position of lattice */
			   int    *count,     /* out:number of steps in each direction */
			   double *step,      /* out:step size in each direction */
			   VectorR directions[]);/* out: vector directions for each index*/

				/* note that user_step, start, count, step and directions
				   are in x,y,z order*/


private void append_new_default_deformation_field(Arg_Data *globals);

private void resample_the_deformation_field(Arg_Data *globals);

public void get_volume_XYZV_indices(Volume data, int xyzv[]);

public void build_default_deformation_field(Arg_Data *globals)
{

  if (get_n_concated_transforms(globals->trans_info.transformation)==1) {
    append_new_default_deformation_field(globals);
  }
  else {	/* we are starting with concatenated transforms,
		   so go get a point to the linear part....      */
    
    resample_the_deformation_field(globals);
  }

}


/* 
   This procedure will replace the existing deformation by a resampled
   equivalent, where the volumetric definition of the new deformation
   field corresponds to the data in *globals, and the values in this 
   new structure must be interpolated from the existing deformation field. 
*/

private void resample_the_deformation_field(Arg_Data *globals)
{

  Volume
    existing_field,
    new_field;
  Real 
    vector_val[3],
    XYZstart[ MAX_DIMENSIONS ], 
    start[    MAX_DIMENSIONS ], 
    XYZstep[  MAX_DIMENSIONS ], 
    step[     MAX_DIMENSIONS ],  
    step2[    MAX_DIMENSIONS ],
    s1[       MAX_DIMENSIONS ],
    voxel[    MAX_DIMENSIONS ];
  int 
    i,
    siz[      MAX_DIMENSIONS ],
    index[    MAX_DIMENSIONS ],
    xyzv[     MAX_DIMENSIONS ],
    XYZcount[ MAX_DIMENSIONS ],
    count[    MAX_DIMENSIONS ];
  General_transform 
    *non_lin_part;
  VectorR 
    XYZdirections[ MAX_DIMENSIONS ];
  Real 
    del_x, del_y, del_z, wx, wy,wz;
  progress_struct
    progress;
  char 
    **data_dim_names;

  
				/* get the nonlinear part
				   of the transformation           */
  
  existing_field = (Volume)NULL;
  non_lin_part = get_nth_general_transform(globals->trans_info.transformation,
					   get_n_concated_transforms(
					       globals->trans_info.transformation)
					   -1); 

  if (get_transform_type( non_lin_part ) == GRID_TRANSFORM)
    existing_field = (Volume)(non_lin_part->displacement_volume);
  else {
    for_less(i,0,get_n_concated_transforms(globals->trans_info.transformation))
      print ("Transform %d is of type %d\n",i, 
	     get_transform_type(
		get_nth_general_transform(globals->trans_info.transformation,
				i) ));
    
    print_error_and_line_num("Cannot find the deformation field transform to resample",
			     __FILE__, __LINE__);
  }

  /* build a vector volume to store the Grid Transform */
  
  new_field = create_volume(4, dim_name_vector_vol, NC_SHORT, TRUE, 0.0, 0.0);

  get_volume_XYZV_indices(new_field, xyzv);

  for_less(i,0,N_DIMENSIONS) 
    step2[i] = globals->step[i];  
				/* get new start, count, step and directions,
				   all returned in X, Y, Z order.          */
  set_up_lattice(existing_field, step2, XYZstart, XYZcount, XYZstep, XYZdirections);

				/* reset count and step to be in volume order */
  for_less(i,0,N_DIMENSIONS) {
    start[ i ] = XYZstart[ i ];
    count[ xyzv[i] ] = XYZcount[ i ];
    step[  xyzv[i] ] = XYZstep[  i ];
  }
				/* add info for the vector dimension */
  count[xyzv[Z+1]] = 3;
  step[xyzv[Z+1]] = 0.0;

         /* use the sign of the step returned to set the true step size */
  for_less(i,0,N_DIMENSIONS) {
    if (step[xyzv[i]]<0) 
      step[xyzv[i]] = -1.0 * ABS(globals->step[i]); 
    else 
      step[xyzv[i]] = ABS(globals->step[i]);
  }

  for_less(i,0,MAX_DIMENSIONS)  /* set the voxel origin, used in the vol def */
    voxel[i] = 0.0; 

  set_volume_sizes(       new_field, count); 
  set_volume_separations( new_field, step);
  set_volume_voxel_range( new_field, -MY_MAX_VOX, MY_MAX_VOX);
  set_volume_real_range(  new_field, -MY_MAX_REAL, MY_MAX_REAL);
  set_volume_translation( new_field, voxel, start);
  
				/* make sure that the vector dimension 
				   is named! */
  data_dim_names = get_volume_dimension_names(new_field);

  if( strcmp( data_dim_names[ xyzv[Z+1] ] , MIvector_dimension ) != 0 ) {
    ALLOC((new_field)->dimension_names[xyzv[Z+1]], \
	  strlen(MIvector_dimension  ) + 1 );
    (void) strcpy( (new_field)->dimension_names[xyzv[Z+1]], MIvector_dimension );
  }
#ifdef HAVE_RECENT_VOLUME_IO
  delete_dimension_names(new_field, data_dim_names);
#else
  delete_dimension_names(data_dim_names);
#endif

  if (globals->flags.debug) {
    print ("in resample_deformation_field:\n");
    print ("xyzv[axes] = %d, %d, %d, %d\n",xyzv[X],xyzv[Y],xyzv[Z],xyzv[Z+1]);

    get_volume_sizes(new_field, siz);
    get_volume_separations(new_field, s1);
    print ("seps: %7.3f %7.3f %7.3f %7.3f %7.3f \n",s1[0],s1[1],s1[2],s1[3],s1[4]);
    print ("size: %7d %7d %7d %7d %7d \n",siz[0],siz[1],siz[2],siz[3],siz[4]);
  }
    
  alloc_volume_data(new_field);
  
  if (globals->flags.verbose>0)
    initialize_progress_report( &progress, FALSE, count[0],
			       "Interpolating new field" );
  
  /* now resample the values from the input deformation */

  for_less(i,0,MAX_DIMENSIONS) {
    voxel[i] = 0.0;
    index[i] = 0;
  }

  for_less(index[xyzv[X]],0,count[xyzv[X]]) {
    voxel[xyzv[X]] = (Real)index[xyzv[X]];

    for_less(index[xyzv[Y]],0,count[xyzv[Y]]) {
      voxel[xyzv[Y]] = (Real)index[xyzv[Y]];

      for_less(index[xyzv[Z]],0,count[xyzv[Z]]) {
	voxel[xyzv[Z]] = (Real)index[xyzv[Z]];

	convert_voxel_to_world(new_field, voxel, &wx,&wy,&wz);


	   grid_transform_point(non_lin_part, wx, wy, wz, 
			     &del_x, &del_y, &del_z);


	
				/* get just the deformation part */
	del_x = del_x - wx;	
	del_y = del_y - wy;
	del_z = del_z - wz;


/*	del_x = del_y = del_z = 0.0;
*/
	vector_val[0] = CONVERT_VALUE_TO_VOXEL(new_field, del_x); 
	vector_val[1] = CONVERT_VALUE_TO_VOXEL(new_field, del_y); 
	vector_val[2] = CONVERT_VALUE_TO_VOXEL(new_field, del_z); 

	for_less(index[ xyzv[ Z+1] ], 0, 3) {
	  SET_VOXEL(new_field, \
		    index[0], index[1], index[2], index[3], index[4], \
		    vector_val[ index[ xyzv[ Z+1] ] ]);
	}
	
	
      }
    }
    if (globals->flags.verbose>0) 
      update_progress_report( &progress, i+1 );
  }
  if (globals->flags.verbose>0) 
    terminate_progress_report( &progress );
  
	 	/* delete and free up old data */
  delete_volume(non_lin_part->displacement_volume);
               /* set new volumes into transform */
  non_lin_part->displacement_volume = new_field;

}



/* 
  this procedure will concatenate a zero-valued non-linear deformation
  field transformation to the end of the globals->trans_info.transformation, 
  using the lattice information (in *globals) to build the deformation 
  volumes.
*/

private void append_new_default_deformation_field(Arg_Data *globals)
{

  Volume
    new_field;  
  Real 
    zero, 
    step[MAX_DIMENSIONS],  
    voxel[MAX_DIMENSIONS], 
    point[N_DIMENSIONS];
  int 
    index[MAX_DIMENSIONS],
    xyzv[MAX_DIMENSIONS],
    i,
    count[MAX_DIMENSIONS],
    count_extended[MAX_DIMENSIONS];
  
  General_transform 
    *grid_trans;
  
  
  /* build a vector volume to store the Grid Transform */
  
  ALLOC(new_field,1);
  
  new_field = create_volume(4, dim_name_vector_vol, NC_SHORT, TRUE, 0.0, 0.0);

  get_volume_XYZV_indices(new_field, xyzv);

				/* get the global voxel count and voxel size */
  for_less(i,0,N_DIMENSIONS) {
    count[xyzv[i]] = globals->count[i];
    count_extended[xyzv[i]] = count[i];
    step[xyzv[i]]  = globals->step[i];
  } 
				/* add info for the vector dimension */
  count[xyzv[Z+1]] = 3;
  count_extended[xyzv[Z+1]] = 3;   
  step[xyzv[Z+1]] = 0.0;

  for_less(i,0,MAX_DIMENSIONS)  /* set the voxel origin, used in the vol def */
    voxel[i] = 0.0; 

  set_volume_sizes(       new_field, count); 
  set_volume_separations( new_field, step);
  set_volume_voxel_range( new_field, -MY_MAX_VOX, MY_MAX_VOX);
  set_volume_real_range(  new_field, -MY_MAX_REAL, MY_MAX_REAL);
  set_volume_translation( new_field, voxel, globals->start);

		   /* now pad the volume along the spatial axis
		      to ensure good coverage of the data space
		      with the deformation field */
 
  for_less(i,0,N_DIMENSIONS) {
    if (globals->count[i]>1) {
      voxel[xyzv[i]] = -2.5;
      count_extended[xyzv[i]] = globals->count[i]+5;
    }
    else {
      voxel[xyzv[i]] = 0.0;
      count_extended[xyzv[i]] = 1;
    }
  }

  if (globals->flags.debug) {
    print("in append_new_default_deformation_field:\n\tcount_extended= %d %d %d %d\n",
	   count_extended[0],count_extended[1],count_extended[2],count_extended[3]);
  }

  set_volume_sizes(new_field, count_extended);
  for_less(i,0,MAX_DIMENSIONS) count[i] = count_extended[i];

              /* reset the first voxel position with the new origin */
  convert_voxel_to_world(new_field, voxel,
			 &(point[X]), &(point[Y]), &(point[Z]));
  for_less(i,0,MAX_DIMENSIONS) voxel[i] = 0;
  set_volume_translation(new_field, voxel, point);
  
              /* allocate space for the deformation field data */
  alloc_volume_data(new_field);

	      /* Initilize the field to zero deformation */
  
  zero = CONVERT_VALUE_TO_VOXEL(new_field, 0.0);

  for_less(index[0],0,count[0])
    for_less(index[1],0,count[1])
      for_less(index[2],0,count[2])
	for_less(index[3],0,count[3])
	  {
	    SET_VOXEL(new_field, index[0],index[1],index[2],index[3],0, zero);
	  }
  
              /* build the new GRID_TRANSFORM */
  ALLOC(grid_trans, 1);
  create_grid_transform(grid_trans, new_field);

              /* append the deforamation to the current transformation */
  concat_general_transforms(globals->trans_info.transformation, grid_trans, 
			    globals->trans_info.transformation);
  
  delete_volume(new_field);

  
}
