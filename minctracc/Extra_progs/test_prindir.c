/* ----------------------------- MNI Header -----------------------------------
@NAME       : test_prindir
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : status
@DESCRIPTION: Program to test calculation of principal directions.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Mon Sep 25 08:45:43 MET 1995
@MODIFIED   : $Log: test_prindir.c,v $
@MODIFIED   : Revision 1.3  2004-02-12 05:54:06  rotor
@MODIFIED   :  * removed public/private defs
@MODIFIED   :
@MODIFIED   : Revision 1.2  2002/03/26 14:15:34  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 1.1  1999/10/25 19:52:11  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :

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
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Extra_progs/test_prindir.c,v 1.3 2004-02-12 05:54:06 rotor Exp $";
#endif

#include <stdio.h>
#include <volume_io/internal_volume_io.h>
#include <time_stamp.h>

/* Constants */
#ifndef TRUE
#  define TRUE 1
#  define FALSE 0
#endif

#define VERY_SMALL_EPS 0.0001	/* this is data dependent! */

static char *default_dim_names[N_DIMENSIONS] = { MIxspace, MIyspace, MIzspace };

BOOLEAN return_local_eigen(Real r[3][3][3],
				  Real dir_1[3],
				  Real dir_2[3],
				  Real dir_3[3],
				  Real val[3]);

BOOLEAN return_principal_directions(Real r[3][3][3],
					   Real dir_1[3],
					   Real dir_2[3],
					   Real *r_K,
					   Real *r_S,
					   Real *r_k1,
					   Real *r_k2,
					   Real *r_norm,
					   Real *r_Lvv,
					   Real eps);

BOOLEAN return_local_eigen_from_hessian(Real r[3][3][3],
				  Real dir_1[3],
				  Real dir_2[3],
				  Real dir_3[3],
				  Real val[3]);

void build_def_field(Volume vol,
		     General_transform **grid_trans);

void init_the_volume_to_zero(Volume volume);

void get_volume_XYZV_indices(Volume data, int xyzv[]);

char *prog_name;

int main(int argc, char *argv[])
{

  progress_struct
    progress;

  Status 
    stat;

  Real 
    tmp, max_val, min_val, intensity_threshold,
    max_val_x, max_val_y, max_val_z,
    min_val_x, min_val_y, min_val_z,
    K, S, k1, k2, Lvv,
    dir_max[3],
    dir_mid[3],
    dir_min[3],
    dir_vals[3],
    val[3][3][3];
  
  Volume 
    data, defs, max;

  General_transform
    *output_trans;

  int
    vectors,
    count,
    index[MAX_DIMENSIONS],
    defs_xyzv[MAX_DIMENSIONS],
    data_xyzv[MAX_DIMENSIONS],
    sizes[MAX_DIMENSIONS],
    m,n,o,i,j,k;

  char 
    *history,
    output_filename[1024];

  prog_name = argv[0];
  history = time_stamp(argc, argv);

  if (argc!=3) {
    print("usage: %s input.mnc output_basename\n", prog_name);
    exit(EXIT_FAILURE);
  }

  stat = input_volume(argv[1],3,default_dim_names, NC_UNSPECIFIED, FALSE, 
		      0.0,0.0,TRUE, &data, (minc_input_options *)NULL);

  if (stat != OK) {
    print ("Error: cannot read %s.\n",argv[1]);
    exit(EXIT_FAILURE);
  }

  build_def_field(data, &output_trans);
  defs = output_trans->displacement_volume;

  max= copy_volume_definition(data, NC_SHORT, FALSE, 0.0, 0.0);

  set_volume_real_range(max, -10.0, 10.0);
  init_the_volume_to_zero(max);

  get_volume_sizes(data,sizes);
  get_volume_XYZV_indices(data,data_xyzv);
  get_volume_XYZV_indices(defs,defs_xyzv);


  max_val = max_val_x = max_val_y = max_val_z = -1000000.0;
  min_val = min_val_x = min_val_y = min_val_z =  1000000.0;


  for_less(i,0,sizes[0])
    for_less(j,0,sizes[1])
      for_less(k,0,sizes[2]) {
	tmp = get_volume_real_value(data, i,j,k,0,0);
	if (tmp>max_val) max_val = tmp;
	if (tmp<min_val) min_val = tmp;
      }

  intensity_threshold = 0.01 * max_val;

  print ("For %s, max = %f, min = %f, thresh = %f\n",argv[1], max_val, min_val, intensity_threshold);

  initialize_progress_report(&progress, FALSE, sizes[0]*sizes[1]*sizes[2]+1,
			     "Calc principal directions");
  vectors = count = 0;

  for_less(index[ defs_xyzv[X] ],1,sizes[data_xyzv[X]]-1)
    for_less(index[ defs_xyzv[Y] ],1,sizes[data_xyzv[Y]]-1)
      for_less(index[ defs_xyzv[Z] ],1,sizes[data_xyzv[Z]]-1) {
	
	tmp = get_volume_real_value(data, 
				    index[defs_xyzv[X]], 
				    index[defs_xyzv[Y]], 
				    index[defs_xyzv[Z]], 0,0);


	dir_max[X] = dir_max[Y] = dir_max[Z] = 0.0;
	dir_min[X] = dir_min[Y] = dir_min[Z] = 0.0;
	K = S = k1 = k2 = 0.0;
	
	if (tmp > intensity_threshold) {

	  for_inclusive(m,-1,1)
	    for_inclusive(n,-1,1)
	      for_inclusive(o,-1,1)
		val[m+1][n+1][o+1] = 
		  get_volume_real_value(data, 
					index[defs_xyzv[X]]+m, 
					index[defs_xyzv[Y]]+n, 
					index[defs_xyzv[Z]]+o, 0,0);
	  
	  if (FALSE) {		/* test principal directions */
	    if (!return_principal_directions(val, 
					     dir_mid, dir_min,
					     &K, &S, &k1, &k2, dir_max, &Lvv,
					     (Real)VERY_SMALL_EPS)) {
	      
	      dir_max[X] = dir_max[Y] = dir_max[Z] = 0.0;
	      dir_min[X] = dir_min[Y] = dir_min[Z] = 0.0;
	      K = S = k1 = k2 = Lvv = 0.0;
	    }
	    
	    tmp = sqrt(dir_max[0] * dir_max[0] +
		       dir_max[1] * dir_max[1] +
		       dir_max[2] * dir_max[2]);
	    
	    if (tmp>0.00001) {
	      dir_max[0] *= 10.0/tmp;
	      dir_max[1] *= 10.0/tmp;
	      dir_max[2] *= 10.0/tmp;
	      vectors++;
	    }

	    
	    if (max_val_x < dir_max[X]) max_val_x=dir_max[X];
	    if (min_val_x > dir_max[X]) min_val_x=dir_max[X];
	    if (max_val_y < dir_max[Y]) max_val_y=dir_max[Y];
	    if (min_val_y > dir_max[Y]) min_val_y=dir_max[Y];
	    if (max_val_z < dir_max[Z]) max_val_z=dir_max[Z];
	    if (min_val_z > dir_max[Z]) min_val_z=dir_max[Z];
	    
	  }
	
	  else {			/* test eigen values */
	    
	    if (return_local_eigen_from_hessian(val, dir_max, dir_mid, dir_min, 
						dir_vals)) {
	      
	      dir_max[0] *= 2.0;
	      dir_max[1] *= 2.0;
	      dir_max[2] *= 2.0;
	      
	      if (max_val_x < dir_max[X]) max_val_x=dir_max[X];
	      if (min_val_x > dir_max[X]) min_val_x=dir_max[X];
	      if (max_val_y < dir_max[Y]) max_val_y=dir_max[Y];
	      if (min_val_y > dir_max[Y]) min_val_y=dir_max[Y];
	      if (max_val_z < dir_max[Z]) max_val_z=dir_max[Z];
	      if (min_val_z > dir_max[Z]) min_val_z=dir_max[Z];	
	      
	      vectors++;

	    }
	    
	  }

	  
	  for_inclusive(index[ defs_xyzv[Z+1]],X,Z) 
	    set_volume_real_value(defs, 
				  index[0], index[1], index[2],
				  index[3], index[5], 
				  dir_max[index[ defs_xyzv[Z+1]] ]);
	  
	  index[ defs_xyzv[Z+1]] = 0;
	  set_volume_real_value(max, 
				index[defs_xyzv[X]], 
				index[defs_xyzv[Y]], 
				index[defs_xyzv[Z]], 0,0, dir_max[0]);
	}
	count++;
	update_progress_report( &progress, count );
	
      }
  terminate_progress_report(&progress);

  print ("%d vectors estimated\n", vectors);

  print ("Saving data...\n");

  sprintf(output_filename,"%s_def.xfm",argv[2]);
  stat = output_transform_file(output_filename, history, output_trans);

  if (stat != OK) {
    print ("Error: cannot write %s.\n",output_filename);
    exit(EXIT_FAILURE);
  }

  sprintf(output_filename,"%s_k1.mnc",argv[2]);
  stat = output_modified_volume(output_filename, NC_UNSPECIFIED, FALSE, 
				0.0, 0.0,  max, argv[1], history, NULL);

  if (stat != OK) {
    print ("Error: cannot write %s.\n",output_filename);
    exit(EXIT_FAILURE);
  }

  print ("done.\n");

  exit(EXIT_SUCCESS);
}



void init_the_volume_to_zero(Volume volume)
{
    int             v0, v1, v2, v3, v4;
    Real            zero;
  
    zero = 0.0;

    BEGIN_ALL_VOXELS( volume, v0, v1, v2, v3, v4 )
      
      set_volume_voxel_value( volume, v0, v1, v2, v3, v4, zero );
    
    END_ALL_VOXELS

}

void build_def_field(Volume vol,
		     General_transform **grid_trans)
{

  Volume
    new_field;  
  Real 
    start[MAX_DIMENSIONS],
    orig_steps[MAX_DIMENSIONS],  
    steps[MAX_DIMENSIONS],  
    voxel[MAX_DIMENSIONS];
  int 
    i,
    orig_xyzv[MAX_DIMENSIONS],
    orig_count[MAX_DIMENSIONS],
    xyzv[MAX_DIMENSIONS],
    count[MAX_DIMENSIONS],
    index[MAX_DIMENSIONS];
  
  
  static  char *dim_name_vector_vol[] = 
    { MIvector_dimension, MIzspace,MIyspace, MIxspace };
  
  /* build a vector volume to store the Grid Transform */
  
  ALLOC(new_field,1);
  
  new_field = create_volume(4, dim_name_vector_vol, NC_SHORT, TRUE, 0.0, 0.0);
  get_volume_XYZV_indices(new_field, xyzv);

				/* get the global voxel count and voxel size */
  get_volume_sizes(       vol, orig_count);
  get_volume_XYZV_indices(vol, orig_xyzv);
  get_volume_separations( vol, orig_steps);

  for_inclusive(i,X,Z) {
    count[xyzv[i]] = orig_count[orig_xyzv[i]];
    steps[ xyzv[i]] = orig_steps[ orig_xyzv[i]];
  } 
				/* add info for the vector dimension */
  count[xyzv[Z+1]] = 3;
  steps[ xyzv[Z+1]] = 0.0;

  for_less(i,0,MAX_DIMENSIONS)  /* set the voxel origin, used in the vol def */
    voxel[i] = 0.0; 

  get_volume_translation( vol, voxel, start);

 
  set_volume_sizes(       new_field, count); 
  set_volume_separations( new_field, steps);
  set_volume_real_range(  new_field, -16.0, 16.0);
  set_volume_starts( new_field,  start);
  
              /* allocate space for the deformation field data */
  alloc_volume_data(new_field);

	      /* Initilize the field to zero deformation */
  init_the_volume_to_zero(new_field);

  ALLOC(*grid_trans, 1);
  create_grid_transform(*grid_trans, new_field);

  delete_volume(new_field);

}

void get_volume_XYZV_indices(Volume data, int xyzv[])
{
  
  int 
    axis, i, vol_dims;
  char 
    **data_dim_names;
  
  vol_dims       = get_volume_n_dimensions(data);
  data_dim_names = get_volume_dimension_names(data);
  
  for_less(i,0,N_DIMENSIONS+1) xyzv[i] = -1;
  for_less(i,0,vol_dims) {
    if (convert_dim_name_to_spatial_axis(data_dim_names[i], &axis )) {
      xyzv[axis] = i; 
    } 
    else {     /* not a spatial axis */
      xyzv[Z+1] = i;
    }
  }
  delete_dimension_names(data_dim_names);
  
}

