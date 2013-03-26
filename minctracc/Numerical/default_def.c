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
@MODIFIED   : Revision 96.13  2009-04-03 18:36:59  louis
@MODIFIED   : made changes to use only DOUBLES for input source and model volumes, and for all estimation of deformation fields
@MODIFIED   :
@MODIFIED   : Revision 96.12  2009/03/13 19:51:31  claude
@MODIFIED   : fixed bug in offsets for minctracc and free memory upon exit
@MODIFIED   :
@MODIFIED   : Revision 96.11  2006/11/30 09:07:32  rotor
@MODIFIED   :  * many more changes for clean minc 2.0 build
@MODIFIED   :
@MODIFIED   : Revision 96.10  2006/11/29 09:09:33  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 96.9  2005/07/20 20:45:49  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 96.8  2004/02/12 05:54:27  rotor
@MODIFIED   :  * removed /static defs
@MODIFIED   :
@MODIFIED   : Revision 96.7  2004/02/04 20:43:46  lenezet
@MODIFIED   :  * change the start to fix the dir cos
@MODIFIED   :
@MODIFIED   : Revision 96.6  2004/01/27 00:28:03  lenezet
@MODIFIED   : change init_params to correct the COG bug when there is not input transform.
@MODIFIED   : add the cosines director to the resampled field
@MODIFIED   :
@MODIFIED   : Revision 96.5  2002/12/13 21:16:30  lenezet
@MODIFIED   : nonlinear in 2D has changed. The option -2D-non-lin is no more necessary. The grid transform has been adapted to feet on the target volume whatever is size. The Optimization is done on the dimensions for which "count" is greater than 1.
@MODIFIED   :
@MODIFIED   : Revision 96.4  2002/11/20 21:38:49  lenezet
@MODIFIED   :
@MODIFIED   : Fix the code to take in consideration the direction cosines especially in the grid transform.
@MODIFIED   : Add an option to choose the maximum expected deformation magnitude.
@MODIFIED   :
@MODIFIED   : Revision 96.3  2002/03/26 14:15:40  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 96.2  2000/02/07 19:33:05  stever
@MODIFIED   : replaced HAVE_RECENT_VOLUME_IO with more specific feature tests.
@MODIFIED   :
@MODIFIED   : Revision 96.1  1997/11/03 19:59:49  louis
@MODIFIED   : - now include volume_io/internal_volume_io.h instead of volume_io.h
@MODIFIED   : - used xyzv[VIO_X] instead of 0 for the index for the x variable
@MODIFIED   :
 * Revision 96.0  1996/08/21  18:21:58  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:21:52  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:15:50  louis
 * Never released version 0.95
 *
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

#include <config.h>
#include <volume_io.h>
#include <Proglib.h>
#include "arg_data.h"
#include "local_macros.h"
#include "constants.h"

#define MY_MAX_VOX 32766.0
#define MY_MAX_REAL 50.0

static  char *dim_name_vector_vol[] =
                { MIvector_dimension, MIzspace,MIyspace, MIxspace };

void set_up_lattice(VIO_Volume data,       /* in: volume  */
                           double *user_step, /* in: user requested spacing for lattice */
                           double *start,     /* out:world starting position of lattice  in volume dircos coords*/
                           double *wstart,     /* out:world starting position of lattice */
                           int    *count,     /* out:number of steps in each direction */
                           double *step,      /* out:step size in each direction */
                           VectorR directions[]);/* out: vector directions for each index*/

                                /* note that user_step, start, count, step and directions
                                   are in x,y,z order*/


static void append_new_default_deformation_field(Arg_Data *globals);

static void resample_the_deformation_field(Arg_Data *globals);

void get_volume_XYZV_indices(VIO_Volume data, int xyzv[]);

void build_default_deformation_field(Arg_Data *globals)
{

  if (get_n_concated_transforms(globals->trans_info.transformation)==1) {

    append_new_default_deformation_field(globals);

  }
  else {        /* we are starting with concatenated transforms,
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

static void resample_the_deformation_field(Arg_Data *globals)
{

  VIO_Volume
    existing_field,
    new_field;
  VIO_Real
    vector_val[3],
    XYZstart[ VIO_MAX_DIMENSIONS ],
    wstart[ VIO_MAX_DIMENSIONS ],
    start[    VIO_MAX_DIMENSIONS ],
    XYZstep[  VIO_MAX_DIMENSIONS ],
    step[     VIO_MAX_DIMENSIONS ],
    step2[    VIO_MAX_DIMENSIONS ],
    s1[       VIO_MAX_DIMENSIONS ],
    voxel[    VIO_MAX_DIMENSIONS ],
    dir[3][3];
  int
    i,
    siz[      VIO_MAX_DIMENSIONS ],
    index[    VIO_MAX_DIMENSIONS ],
    xyzv[     VIO_MAX_DIMENSIONS ],
    XYZcount[ VIO_MAX_DIMENSIONS ],
    count[    VIO_MAX_DIMENSIONS ];
  VIO_General_transform
    *non_lin_part;
  VectorR
    XYZdirections[ VIO_MAX_DIMENSIONS ];
  VIO_Real
    del_x, del_y, del_z, wx, wy,wz;
  VIO_progress_struct
    progress;
  char
    **data_dim_names;


                                /* get the nonlinear part
                                   of the transformation           */

  existing_field = (VIO_Volume)NULL;
  non_lin_part = get_nth_general_transform(globals->trans_info.transformation,
                                           get_n_concated_transforms(
                                               globals->trans_info.transformation)
                                           -1);

  if (get_transform_type( non_lin_part ) == GRID_TRANSFORM){
    existing_field = (VIO_Volume)(non_lin_part->displacement_volume);
  }
  else {
    for(i=0; i<get_n_concated_transforms(globals->trans_info.transformation); i++)
      print ("Transform %d is of type %d\n",i,
             get_transform_type(
                get_nth_general_transform(globals->trans_info.transformation,
                                i) ));

    print_error_and_line_num("Cannot find the deformation field transform to resample",
                             __FILE__, __LINE__);
  }

  /* build a vector volume to store the Grid VIO_Transform */

  new_field = create_volume(4, dim_name_vector_vol, NC_DOUBLE, TRUE, 0.0, 0.0);

  get_volume_XYZV_indices(new_field, xyzv);

  for(i=0; i<VIO_N_DIMENSIONS; i++)
    step2[i] = globals->step[i];
                                /* get new start, count, step and directions,
                                   all returned in X, Y, Z order.          */

  set_up_lattice(existing_field, step2, XYZstart, wstart, XYZcount, XYZstep, XYZdirections);

                                /* reset count and step to be in volume order */
  for(i=0; i<VIO_N_DIMENSIONS; i++) {
    start[      i  ] = wstart[ i ];
    count[ xyzv[i] ] = XYZcount[ i ];
    step[  xyzv[i] ] = XYZstep[  i ];
  }

                                /* add info for the vector dimension */
  count[xyzv[VIO_Z+1]] = 3;
  step[xyzv[VIO_Z+1]] = 0.0;

         /* use the sign of the step returned to set the true step size */
  for(i=0; i<VIO_N_DIMENSIONS; i++) {
    if (step[xyzv[i]]<0)
      step[xyzv[i]] = -1.0 * fabs(globals->step[i]);
    else
      step[xyzv[i]] = fabs(globals->step[i]);
  }

  for(i=0; i<VIO_MAX_DIMENSIONS; i++)  /* set the voxel origin, used in the vol def */
    voxel[i] = 0.0;

  set_volume_sizes(       new_field, count);
  set_volume_separations( new_field, step);

  /*  set_volume_voxel_range( new_field, -MY_MAX_VOX, MY_MAX_VOX);
      set_volume_real_range(  new_field, -1.0*globals->trans_info.max_def_magnitude, globals->trans_info.max_def_magnitude);  - no longer needed, because now using doubles*/

  set_volume_translation( new_field, voxel, start);

  for(i=0; i<VIO_N_DIMENSIONS; i++) {
    dir[VIO_X][i]=XYZdirections[VIO_X].coords[i];
    dir[VIO_Y][i]=XYZdirections[VIO_Y].coords[i];
    dir[VIO_Z][i]=XYZdirections[VIO_Z].coords[i];

  }


  set_volume_direction_cosine(new_field,xyzv[VIO_X],dir[VIO_X]);
  set_volume_direction_cosine(new_field,xyzv[VIO_Y],dir[VIO_Y]);
  set_volume_direction_cosine(new_field,xyzv[VIO_Z],dir[VIO_Z]);


                                /* make sure that the vector dimension
                                   is named! */
  data_dim_names = get_volume_dimension_names(new_field);

  if( strcmp( data_dim_names[ xyzv[VIO_Z+1] ] , MIvector_dimension ) != 0 ) {
    ALLOC((new_field)->dimension_names[xyzv[VIO_Z+1]], \
          strlen(MIvector_dimension  ) + 1 );
    (void) strcpy( (new_field)->dimension_names[xyzv[VIO_Z+1]], MIvector_dimension );
  }

  delete_dimension_names(new_field, data_dim_names);

  if (globals->flags.debug) {
    print ("in resample_deformation_field:\n");
    print ("xyzv[axes] = %d, %d, %d, %d\n",xyzv[VIO_X],xyzv[VIO_Y],xyzv[VIO_Z],xyzv[VIO_Z+1]);

    get_volume_sizes(new_field, siz);
    get_volume_separations(new_field, s1);
    print ("seps: %7.3f %7.3f %7.3f %7.3f %7.3f \n",s1[0],s1[1],s1[2],s1[3],s1[4]);
    print ("size: %7d %7d %7d %7d %7d \n",siz[0],siz[1],siz[2],siz[3],siz[4]);
  }

  alloc_volume_data(new_field);

  if (globals->flags.verbose>0)
    initialize_progress_report( &progress, FALSE, count[xyzv[VIO_X]],
                               "Interpolating new field" );

  /* now resample the values from the input deformation */

  for(i=0; i<VIO_MAX_DIMENSIONS; i++) {
    voxel[i] = 0.0;
    index[i] = 0;
  }

  for(index[xyzv[VIO_X]]=0; index[xyzv[VIO_X]]<count[xyzv[VIO_X]]; index[xyzv[VIO_X]]++) {
    voxel[xyzv[VIO_X]] = (VIO_Real)index[xyzv[VIO_X]];

    for(index[xyzv[VIO_Y]]=0; index[xyzv[VIO_Y]]<count[xyzv[VIO_Y]]; index[xyzv[VIO_Y]]++) {
      voxel[xyzv[VIO_Y]] = (VIO_Real)index[xyzv[VIO_Y]];

      for(index[xyzv[VIO_Z]]=0; index[xyzv[VIO_Z]]<count[xyzv[VIO_Z]]; index[xyzv[VIO_Z]]++) {
        voxel[xyzv[VIO_Z]] = (VIO_Real)index[xyzv[VIO_Z]];

        convert_voxel_to_world(new_field, voxel, &wx,&wy,&wz);


           grid_transform_point(non_lin_part, wx, wy, wz,
                             &del_x, &del_y, &del_z);



                                /* get just the deformation part */
        del_x = del_x - wx;
        del_y = del_y - wy;
        del_z = del_z - wz;


/*        del_x = del_y = del_z = 0.0;
*/
        vector_val[0] = CONVERT_VALUE_TO_VOXEL(new_field, del_x);
        vector_val[1] = CONVERT_VALUE_TO_VOXEL(new_field, del_y);
        vector_val[2] = CONVERT_VALUE_TO_VOXEL(new_field, del_z);

        for(index[xyzv[VIO_Z+1]]=0; index[xyzv[VIO_Z+1]]<3; index[xyzv[VIO_Z+1]]++) {
          SET_VOXEL(new_field, \
                    index[0], index[1], index[2], index[3], index[4], \
                    vector_val[ index[ xyzv[ VIO_Z+1] ] ]);
        }


      }
    }
    if (globals->flags.verbose>0)
      update_progress_report( &progress,index[xyzv[VIO_X]]+1);
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

static void append_new_default_deformation_field(Arg_Data *globals)
{

  VIO_Volume
    new_field;
  VIO_Real
    zero,
    st[VIO_MAX_DIMENSIONS],
    wst[VIO_MAX_DIMENSIONS],
    step[VIO_MAX_DIMENSIONS],
    XYZstart[ VIO_MAX_DIMENSIONS ],
    XYZstep[ VIO_MAX_DIMENSIONS ],
    voxel[VIO_MAX_DIMENSIONS],
    point[VIO_N_DIMENSIONS],
    dir[3][3];

  int
    index[VIO_MAX_DIMENSIONS],
    xyzv[VIO_MAX_DIMENSIONS],
    i,
    count[VIO_MAX_DIMENSIONS],
    XYZcount[VIO_MAX_DIMENSIONS],
    count_extended[VIO_MAX_DIMENSIONS];

  VIO_General_transform
    *grid_trans;

   VectorR
    XYZdirections[ VIO_MAX_DIMENSIONS ];

  /* build a vector volume to store the Grid VIO_Transform */

   /*  ALLOC(new_field,1); not needed since create volume allocs it
       internally and returns a pointer*/

  if (globals->flags.debug) { print ("In append_new_default_deformation_field...\n"); }

  new_field = create_volume(4, dim_name_vector_vol, NC_DOUBLE, TRUE, 0.0, 0.0);

  get_volume_XYZV_indices(new_field, xyzv);

                                /* get the global voxel count and voxel size */
  for(i=0; i<VIO_N_DIMENSIONS; i++) {
    count[xyzv[i]] = globals->count[i];
    count_extended[xyzv[i]] = count[xyzv[i]];
    step[xyzv[i]]  = globals->step[i];
  }
                                /* add info for the vector dimension */
  count[xyzv[VIO_Z+1]] = 3;
  count_extended[xyzv[VIO_Z+1]] = 3;
  step[xyzv[VIO_Z+1]] = 0.0;


  set_volume_sizes(       new_field, count);
  set_volume_separations( new_field, step);
  /*
     set_volume_voxel_range( new_field, -MY_MAX_VOX, MY_MAX_VOX);
     set_volume_real_range(  new_field, -1.0*globals->trans_info.max_def_magnitude, globals->trans_info.max_def_magnitude); no longer needed, now using floats */


  for(i=0; i<VIO_N_DIMENSIONS; i++) {
    dir[VIO_X][i]=globals->directions[VIO_X].coords[i];
    dir[VIO_Y][i]=globals->directions[VIO_Y].coords[i];
    dir[VIO_Z][i]=globals->directions[VIO_Z].coords[i];

  }


  set_volume_direction_cosine(new_field,xyzv[VIO_X],dir[VIO_X]);
  set_volume_direction_cosine(new_field,xyzv[VIO_Y],dir[VIO_Y]);
  set_volume_direction_cosine(new_field,xyzv[VIO_Z],dir[VIO_Z]);


  for(i=0; i<VIO_MAX_DIMENSIONS; i++)  /* set the voxel origin, used in the vol def */
    voxel[i] = 0.0;

  set_volume_translation( new_field, voxel, globals->start);



  if (globals->flags.debug) {
    print("in append new def, the start is: %8.3f %8.3f %8.3f\n", globals->start[VIO_X], globals->start[VIO_Y], globals->start[VIO_Z]);
  }

                   /* now pad the volume along the spatial axis
                      to ensure good coverage of the data space
                      with the deformation field */

  for(i=0; i<VIO_N_DIMENSIONS; i++) {
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
  for(i=0; i<VIO_MAX_DIMENSIONS; i++) count[i] = count_extended[i];

              /* reset the first voxel position with the new origin */
  convert_voxel_to_world(new_field, voxel,
                         &(point[VIO_X]), &(point[VIO_Y]), &(point[VIO_Z]));
  for(i=0; i<VIO_MAX_DIMENSIONS; i++) voxel[i] = 0;

  set_volume_translation(new_field, voxel, point);


  if (globals->flags.debug) {
    print (" point: %8.3f %8.3f %8.3f \n", point[VIO_X], point[VIO_Y], point[VIO_Z]);

    get_volume_starts(new_field, st);
    print (" start: %8.3f %8.3f %8.3f \n", st[xyzv[VIO_X]], st[xyzv[VIO_Y]], st[xyzv[VIO_Z]]);

    voxel[0] = 0;
    voxel[1] = 0;
    voxel[2] = 0;
    get_volume_translation(new_field, voxel, wst);
    print (" wstrt: %8.3f %8.3f %8.3f \n", wst[VIO_X], wst[VIO_Y], wst[VIO_Z]);
    print (" voxel: %8.3f %8.3f %8.3f \n", voxel[xyzv[VIO_X]], voxel[xyzv[VIO_Y]], voxel[xyzv[VIO_Z]]);


    for(i=0; i<3; i++) {
      get_volume_direction_cosine(new_field,xyzv[i], wst);
      print (" dirs: %8.3f %8.3f %8.3f \n", wst[VIO_X], wst[VIO_Y], wst[VIO_Z]);
    }



  }


              /* allocate space for the deformation field data */
  alloc_volume_data(new_field);

              /* Initilize the field to zero deformation */

  /* zero = CONVERT_VALUE_TO_VOXEL(new_field, 0.0); not needed, defs are now doubles */

  for(index[0]=0; index[0]<count[0]; index[0]++)
    for(index[1]=0; index[1]<count[1]; index[1]++)
      for(index[2]=0; index[2]<count[2]; index[2]++)
        for(index[3]=0; index[3]<count[3]; index[3]++)
          {
            SET_VOXEL(new_field, index[0],index[1],index[2],index[3],0, 0.0);  /* was set to 'zero', but now as a double,can be set to 0.0 */
          }

              /* build the new GRID_TRANSFORM */

  ALLOC(grid_trans, 1);

  create_grid_transform(grid_trans, new_field, NULL);

              /* append the deforamation to the current transformation */
  concat_general_transforms(globals->trans_info.transformation, grid_trans,
                            globals->trans_info.transformation);

  delete_volume(new_field);
  delete_general_transform(grid_trans);




}
