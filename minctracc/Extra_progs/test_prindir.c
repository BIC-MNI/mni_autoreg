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
@MODIFIED   : Revision 1.5  2006-11-29 09:09:31  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 1.4  2005/07/20 20:45:47  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 1.3  2004/02/12 05:54:06  rotor
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
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Extra_progs/test_prindir.c,v 1.5 2006-11-29 09:09:31 rotor Exp $";
#endif

#include <stdio.h>
#include <volume_io.h>
#include <time_stamp.h>

/* Constants */
#ifndef TRUE
#  define TRUE 1
#  define FALSE 0
#endif

#define VERY_SMALL_EPS 0.0001        /* this is data dependent! */

static char *default_dim_names[VIO_N_DIMENSIONS] = { MIxspace, MIyspace, MIzspace };

VIO_BOOL return_local_eigen(VIO_Real r[3][3][3],
                                  VIO_Real dir_1[3],
                                  VIO_Real dir_2[3],
                                  VIO_Real dir_3[3],
                                  VIO_Real val[3]);

VIO_BOOL return_principal_directions(VIO_Real r[3][3][3],
                                           VIO_Real dir_1[3],
                                           VIO_Real dir_2[3],
                                           VIO_Real *r_K,
                                           VIO_Real *r_S,
                                           VIO_Real *r_k1,
                                           VIO_Real *r_k2,
                                           VIO_Real *r_norm,
                                           VIO_Real *r_Lvv,
                                           VIO_Real eps);

VIO_BOOL return_local_eigen_from_hessian(VIO_Real r[3][3][3],
                                  VIO_Real dir_1[3],
                                  VIO_Real dir_2[3],
                                  VIO_Real dir_3[3],
                                  VIO_Real val[3]);

void build_def_field(VIO_Volume vol,
                     VIO_General_transform **grid_trans);

void init_the_volume_to_zero(VIO_Volume volume);

void get_volume_XYZV_indices(VIO_Volume data, int xyzv[]);

char *prog_name;

int main(int argc, char *argv[])
{

  VIO_progress_struct
    progress;

  VIO_Status 
    stat;

  VIO_Real 
    tmp, max_val, min_val, intensity_threshold,
    max_val_x, max_val_y, max_val_z,
    min_val_x, min_val_y, min_val_z,
    K, S, k1, k2, Lvv,
    dir_max[3],
    dir_mid[3],
    dir_min[3],
    dir_vals[3],
    val[3][3][3];
  
  VIO_Volume 
    data, defs, max;

  VIO_General_transform
    *output_trans;

  int
    vectors,
    count,
    index[VIO_MAX_DIMENSIONS],
    defs_xyzv[VIO_MAX_DIMENSIONS],
    data_xyzv[VIO_MAX_DIMENSIONS],
    sizes[VIO_MAX_DIMENSIONS],
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


  for(i=0; i<sizes[0]; i++)
    for(j=0; j<sizes[1]; j++)
      for(k=0; k<sizes[2]; k++) {
        tmp = get_volume_real_value(data, i,j,k,0,0);
        if (tmp>max_val) max_val = tmp;
        if (tmp<min_val) min_val = tmp;
      }

  intensity_threshold = 0.01 * max_val;

  print ("For %s, max = %f, min = %f, thresh = %f\n",argv[1], max_val, min_val, intensity_threshold);

  initialize_progress_report(&progress, FALSE, sizes[0]*sizes[1]*sizes[2]+1,
                             "Calc principal directions");
  vectors = count = 0;

  for(index[defs_xyzv[VIO_X]]=1; index[defs_xyzv[VIO_X]]<sizes[data_xyzv[VIO_X]]-1; index[defs_xyzv[VIO_X]]++)
    for(index[defs_xyzv[VIO_Y]]=1; index[defs_xyzv[VIO_Y]]<sizes[data_xyzv[VIO_Y]]-1; index[defs_xyzv[VIO_Y]]++)
      for(index[defs_xyzv[VIO_Z]]=1; index[defs_xyzv[VIO_Z]]<sizes[data_xyzv[VIO_Z]]-1; index[defs_xyzv[VIO_Z]]++) {
        
        tmp = get_volume_real_value(data, 
                                    index[defs_xyzv[VIO_X]], 
                                    index[defs_xyzv[VIO_Y]], 
                                    index[defs_xyzv[VIO_Z]], 0,0);


        dir_max[VIO_X] = dir_max[VIO_Y] = dir_max[VIO_Z] = 0.0;
        dir_min[VIO_X] = dir_min[VIO_Y] = dir_min[VIO_Z] = 0.0;
        K = S = k1 = k2 = 0.0;
        
        if (tmp > intensity_threshold) {

          for(m=-1; m<=1; m++)
            for(n=-1; n<=1; n++)
              for(o=-1; o<=1; o++)
                val[m+1][n+1][o+1] = 
                  get_volume_real_value(data, 
                                        index[defs_xyzv[VIO_X]]+m, 
                                        index[defs_xyzv[VIO_Y]]+n, 
                                        index[defs_xyzv[VIO_Z]]+o, 0,0);
          
          if (FALSE) {                /* test principal directions */
            if (!return_principal_directions(val, 
                                             dir_mid, dir_min,
                                             &K, &S, &k1, &k2, dir_max, &Lvv,
                                             (VIO_Real)VERY_SMALL_EPS)) {
              
              dir_max[VIO_X] = dir_max[VIO_Y] = dir_max[VIO_Z] = 0.0;
              dir_min[VIO_X] = dir_min[VIO_Y] = dir_min[VIO_Z] = 0.0;
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

            
            if (max_val_x < dir_max[VIO_X]) max_val_x=dir_max[VIO_X];
            if (min_val_x > dir_max[VIO_X]) min_val_x=dir_max[VIO_X];
            if (max_val_y < dir_max[VIO_Y]) max_val_y=dir_max[VIO_Y];
            if (min_val_y > dir_max[VIO_Y]) min_val_y=dir_max[VIO_Y];
            if (max_val_z < dir_max[VIO_Z]) max_val_z=dir_max[VIO_Z];
            if (min_val_z > dir_max[VIO_Z]) min_val_z=dir_max[VIO_Z];
            
          }
        
          else {                        /* test eigen values */
            
            if (return_local_eigen_from_hessian(val, dir_max, dir_mid, dir_min, 
                                                dir_vals)) {
              
              dir_max[0] *= 2.0;
              dir_max[1] *= 2.0;
              dir_max[2] *= 2.0;
              
              if (max_val_x < dir_max[VIO_X]) max_val_x=dir_max[VIO_X];
              if (min_val_x > dir_max[VIO_X]) min_val_x=dir_max[VIO_X];
              if (max_val_y < dir_max[VIO_Y]) max_val_y=dir_max[VIO_Y];
              if (min_val_y > dir_max[VIO_Y]) min_val_y=dir_max[VIO_Y];
              if (max_val_z < dir_max[VIO_Z]) max_val_z=dir_max[VIO_Z];
              if (min_val_z > dir_max[VIO_Z]) min_val_z=dir_max[VIO_Z];        
              
              vectors++;

            }
            
          }

          
          for(index[defs_xyzv[Z+1]]=X; index[defs_xyzv[Z+1]]<=Z; index[defs_xyzv[Z+1]]++) 
            set_volume_real_value(defs, 
                                  index[0], index[1], index[2],
                                  index[3], index[5], 
                                  dir_max[index[ defs_xyzv[Z+1]] ]);
          
          index[ defs_xyzv[Z+1]] = 0;
          set_volume_real_value(max, 
                                index[defs_xyzv[VIO_X]], 
                                index[defs_xyzv[VIO_Y]], 
                                index[defs_xyzv[VIO_Z]], 0,0, dir_max[0]);
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



void init_the_volume_to_zero(VIO_Volume volume)
{
    int             v0, v1, v2, v3, v4;
    VIO_Real            zero;
  
    zero = 0.0;

    BEGIN_ALL_VOXELS( volume, v0, v1, v2, v3, v4 )
      
      set_volume_voxel_value( volume, v0, v1, v2, v3, v4, zero );
    
    END_ALL_VOXELS

}

void build_def_field(VIO_Volume vol,
                     VIO_General_transform **grid_trans)
{

  VIO_Volume
    new_field;  
  VIO_Real 
    start[VIO_MAX_DIMENSIONS],
    orig_steps[VIO_MAX_DIMENSIONS],  
    steps[VIO_MAX_DIMENSIONS],  
    voxel[VIO_MAX_DIMENSIONS];
  int 
    i,
    orig_xyzv[VIO_MAX_DIMENSIONS],
    orig_count[VIO_MAX_DIMENSIONS],
    xyzv[VIO_MAX_DIMENSIONS],
    count[VIO_MAX_DIMENSIONS],
    index[VIO_MAX_DIMENSIONS];
  
  
  static  char *dim_name_vector_vol[] = 
    { MIvector_dimension, MIzspace,MIyspace, MIxspace };
  
  /* build a vector volume to store the Grid VIO_Transform */
  
  ALLOC(new_field,1);
  
  new_field = create_volume(4, dim_name_vector_vol, NC_SHORT, TRUE, 0.0, 0.0);
  get_volume_XYZV_indices(new_field, xyzv);

                                /* get the global voxel count and voxel size */
  get_volume_sizes(       vol, orig_count);
  get_volume_XYZV_indices(vol, orig_xyzv);
  get_volume_separations( vol, orig_steps);

  for(i=X; i<=Z; i++) {
    count[xyzv[i]] = orig_count[orig_xyzv[i]];
    steps[ xyzv[i]] = orig_steps[ orig_xyzv[i]];
  } 
                                /* add info for the vector dimension */
  count[xyzv[Z+1]] = 3;
  steps[ xyzv[Z+1]] = 0.0;

  for(i=0; i<MAX_DIMENSIONS; i++)  /* set the voxel origin, used in the vol def */
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

void get_volume_XYZV_indices(VIO_Volume data, int xyzv[])
{
  
  int 
    axis, i, vol_dims;
  char 
    **data_dim_names;
  
  vol_dims       = get_volume_n_dimensions(data);
  data_dim_names = get_volume_dimension_names(data);
  
  for(i=0; i<N_DIMENSIONS+1; i++) xyzv[i] = -1;
  for(i=0; i<vol_dims; i++) {
    if (convert_dim_name_to_spatial_axis(data_dim_names[i], &axis )) {
      xyzv[axis] = i; 
    } 
    else {     /* not a spatial axis */
      xyzv[Z+1] = i;
    }
  }
  delete_dimension_names(data_dim_names);
  
}

