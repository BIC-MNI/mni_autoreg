/* ----------------------------- MNI Header -----------------------------------
@NAME       : volume_functions.c
@DESCRIPTION: collection of routines used to manipulate volume data.
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

@CREATED    : Tue Jun 15 08:57:23 EST 1993 LC
@MODIFIED   :  $Log: volume_functions.c,v $
@MODIFIED   :  Revision 96.14  2009-04-03 18:36:59  louis
@MODIFIED   :  made changes to use only DOUBLES for input source and model volumes, and for all estimation of deformation fields
@MODIFIED   :
@MODIFIED   :  Revision 96.13  2006/11/30 09:07:33  rotor
@MODIFIED   :   * many more changes for clean minc 2.0 build
@MODIFIED   :
@MODIFIED   :  Revision 96.11  2005/07/20 20:45:52  rotor
@MODIFIED   :      * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :      * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :      * Removed old VOLUME_IO cruft #defines
@MODIFIED   :      * Fixed up all Makefile.am's in subdirs
@MODIFIED   :      * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :      * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   :  Revision 96.10  2004/02/12 06:04:53  rotor
@MODIFIED   :   * removed /static defs
@MODIFIED   :
@MODIFIED   :  Revision 96.9  2002/03/26 14:15:47  stever
@MODIFIED   :  Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   :  Revision 96.8  2001/01/17 14:26:07  louis
@MODIFIED   :  added a delete_volume() after the copy_volume_definition_no_alloc,
@MODIFIED   :  just to make sure that there is no leak.
@MODIFIED   :
@MODIFIED   :  Revision 96.7  2001/01/17 14:10:38  louis
@MODIFIED   :  Fixed memory leak in the code used to compute the intensity
@MODIFIED   :  normalization function.  I used a copy_volume_definition instead of
@MODIFIED   :  copy_volume_definition_no_alloc to build a volume header structure.
@MODIFIED   :
@MODIFIED   :  Revision 96.6  2000/05/23 16:29:14  louis
@MODIFIED   :  Fixed index ordering in the intensity normalization procedure
@MODIFIED   :
@MODIFIED   :  Revision 96.5  2000/05/16 19:48:09  louis
@MODIFIED   :  adjusting code for optical flow
@MODIFIED   :
@MODIFIED   :  Revision 96.4  2000/05/16 15:50:01  louis
@MODIFIED   :  Now calls set_feature_value_threshold to set default values before calling
@MODIFIED   :  volume intensity normalization.
@MODIFIED   :  Also fixed qsort for vol int normlization
@MODIFIED   :
@MODIFIED   :  Revision 96.3  2000/05/08 17:40:45  louis
@MODIFIED   :  addes a dummy change to VIO_Volume/volume_functions.c to test cvs and its
@MODIFIED   :  suspected removal of all working files.
@MODIFIED   :
@MODIFIED   :  Revision 96.2  2000/05/05 17:57:06  louis
@MODIFIED   :  addes volume intensity normalization code
@MODIFIED   :
@MODIFIED   :  Revision 96.1  1999/10/25 19:59:17  louis
@MODIFIED   :  final checkin before switch to CVS
@MODIFIED   :
 * Revision 96.0  1996/08/21  18:22:15  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.5  1996/08/12  14:16:15  louis
 * Release of MNI_AutoReg version 1.0
 *
 * Revision 1.11  1996/08/12  14:16:13  louis
 * Pre-release
 *
 * Revision 1.10  1995/02/22  08:56:06  collins
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.9  94/04/26  12:54:44  louis
 * updated with new versions of make_rots, extract2_parameters_from_matrix 
 * that include proper interpretation of skew.
 * 
 * Revision 1.8  94/04/06  11:49:00  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.7  94/02/21  16:37:46  louis
 * version before feb 22 changes
 * 
 * Revision 1.6  93/11/15  16:27:13  louis
 * working version, with new library, with RCS revision stuff,
 * before deformations included
 * 
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/minctracc/VIO_Volume/volume_functions.c,v 96.14 2009-04-03 18:36:59 louis Exp $";
#endif

#include <config.h>
#include <float.h>
#include <volume_io.h>
#include "point_vector.h"
#include "constants.h"
#include <Proglib.h>
#include <arg_data.h>                /* definition of the global data struct      */
#include "local_macros.h"

int point_not_masked(VIO_Volume volume, 
                            VIO_Real wx, VIO_Real wy, VIO_Real wz);
VIO_Real get_value_of_point_in_volume(VIO_Real xw, VIO_Real yw, VIO_Real zw, 
          VIO_Volume data);


void  set_feature_value_threshold(VIO_Volume d1, 
                                         VIO_Volume d2,
                                         VIO_Real *global_thres1, 
                                         VIO_Real *global_thres2, 
                                         VIO_Real *threshold1, 
                                         VIO_Real *threshold2);


#define MIN_ZRANGE -5.0
#define MAX_ZRANGE  5.0

void make_zscore_volume(VIO_Volume d1, VIO_Volume m1, 
                               VIO_Real *threshold)
{
  unsigned long
    count;
  int 
    stat_count,
    sizes[VIO_MAX_DIMENSIONS],
    s,r,c;
  VIO_Real
    wx,wy,wz,
    valid_min_dvoxel, valid_max_dvoxel,
    min,max,
    sum, sum2, mean, var, std,
    data_vox,data_val,
    thick[VIO_MAX_DIMENSIONS];

  PointR 
    voxel;

  VIO_Volume 
    vol;

  VIO_progress_struct
    progress;

  /* get default information from data and mask */

  /* build temporary working volume */
 
  vol = copy_volume_definition(d1, NC_UNSPECIFIED, FALSE, 0.0, 0.0);
  set_volume_real_range(vol, MIN_ZRANGE, MAX_ZRANGE);
  get_volume_sizes(d1, sizes);
  get_volume_separations(d1, thick);
  get_volume_voxel_range(d1, &valid_min_dvoxel, &valid_max_dvoxel);

  /* initialize counters and sums */

  count  = 0;
  sum  = 0.0;
  sum2 = 0.0;
  min = 1e38;
  max = -1e38;
  stat_count = 0;

  initialize_progress_report(&progress, FALSE, sizes[0]*sizes[1]*sizes[2] + 1,
                             "Tally stats" );

                                /* do first pass, to get mean and std */
  for(s=0; s<sizes[0]; s++) {
    for(r=0; r<sizes[1]; r++) {
      for(c=0; c<sizes[2]; c++) {

        stat_count++;
        update_progress_report( &progress, stat_count);
        convert_3D_voxel_to_world(d1, (VIO_Real)s, (VIO_Real)r, (VIO_Real)c, &wx, &wy, &wz);

        if (m1 != NULL) {
          convert_3D_world_to_voxel(m1, wx, wy, wz, &Point_x(voxel), &Point_y(voxel), &Point_z(voxel));
        }
        else {
          wx = 0.0; wy = 0.0; wz = 0.0;
        }

        if (point_not_masked(m1, wx,wy,wz)) {
          
          GET_VOXEL_3D( data_vox,  d1 , s, r, c );

          if (data_vox >= valid_min_dvoxel && data_vox <= valid_max_dvoxel) { 

            data_val = CONVERT_VOXEL_TO_VALUE(d1, data_vox);
            
            if (data_val > *threshold) {
              sum  += data_val;
              sum2 += data_val*data_val;
              
              count++;
              
              if (data_val < min)
                min = data_val;
              else
                if (data_val > max)
                  max = data_val;
            }
          }
        }
      }
    }
  }
  terminate_progress_report( &progress );

  stat_count = 0;
  initialize_progress_report(&progress, FALSE, sizes[0]*sizes[1]*sizes[2] + 1,
                             "Zscore convert" );

                                /* calc mean and std */
  mean = sum / (float)count;
  var  = ((float)count*sum2 - sum*sum) / ((float)count*((float)count-1));
  std  = sqrt(var);

  min = 1e38;
  max = -1e38;

                                /* replace the voxel values */
  for(s=0; s<sizes[0]; s++) {
    for(r=0; r<sizes[1]; r++) {
      for(c=0; c<sizes[2]; c++) {
        
        stat_count++;
        update_progress_report( &progress, stat_count);

        GET_VOXEL_3D( data_vox,  d1, s, r, c );
        
        if (data_vox >= valid_min_dvoxel && data_vox <= valid_max_dvoxel) { 
          
          data_val = CONVERT_VOXEL_TO_VALUE(d1, data_vox);
          
          if (data_val > *threshold) {

                                /* instead of   
                                   data_val = CONVERT_VALUE_TO_VOXEL(d1, data_vox);
                                   i will use
                                   data_val = CONVERT_VALUE_TO_VOXEL(d1, vol);

                                   since the values in vol are changed with respect to the
                                   new z-score volume */

            data_val = (data_val - mean) / std;
            if (data_val< MIN_ZRANGE) data_val = MIN_ZRANGE;
            if (data_val> MAX_ZRANGE) data_val = MAX_ZRANGE;

            data_vox = CONVERT_VALUE_TO_VOXEL( vol, data_val);
            

            if (data_val < min) {
              min = data_val;
            }
            else {
              if (data_val > max)
                max = data_val;
            }
          }
          else
            data_vox = -DBL_MAX;   /* should be fill_value! */
          
          SET_VOXEL_3D( d1 , s, r, c, data_vox );
        }
        
      }
    }
  }

  terminate_progress_report( &progress );

  set_volume_real_range(d1, MIN_ZRANGE, MAX_ZRANGE);        /* reset the data volume's range */

  *threshold = (*threshold - mean) / std;

  delete_volume(vol);
  
}

void add_speckle_to_volume(VIO_Volume d1, 
                                  float speckle,
                                  double  *start, int *count, VectorR directions[])
{
  VectorR
    vector_step;

  PointR 
    starting_position,
    slice,
    row,
    col;
  VIO_Real valid_min_voxel, valid_max_voxel;

  double
    tx,ty,tz,
    voxel_value;
  int
    xi,yi,zi,
    flip_flag,r,c,s;


  flip_flag = FALSE;

  get_volume_voxel_range(d1, &valid_min_voxel, &valid_max_voxel);
  fill_Point( starting_position, start[0], start[1], start[2]);
  
  for(s=0; s<count[SLICE_IND]; s++) {

    SCALE_VECTOR( vector_step, directions[SLICE_IND], s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    for(r=0; r<count[ROW_IND]; r++) {

      SCALE_VECTOR( vector_step, directions[ROW_IND], r);
      ADD_POINT_VECTOR( row, slice, vector_step );

      SCALE_POINT( col, row, 1.0); /* init first col position */
      for(c=0; c<count[COL_IND]; c++) {

        convert_3D_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);

        xi = (int)floor(tx + 0.5);
        yi = (int)floor( ty + 0.5);
        zi = (int)floor( tz + 0.5 );

        GET_VOXEL_3D( voxel_value, d1 , xi, yi, zi ); 


        if (voxel_value >= valid_min_voxel && voxel_value <= valid_max_voxel) {
          if (flip_flag)
            voxel_value = voxel_value * (1 + 0.01*speckle);
          else
            voxel_value = voxel_value * (1 - 0.01*speckle);

          SET_VOXEL_3D( d1 , xi, yi, zi, voxel_value );
        }

        flip_flag = !flip_flag;


        ADD_POINT_VECTOR( col, col, directions[COL_IND] );

      }
        


    }
  }


}


void save_volume(VIO_Volume d, char *filename)
{
  VIO_Status status;

  status = output_volume(filename,NC_UNSPECIFIED, FALSE, 0.0, 0.0, d, (char *)NULL,
                         (minc_output_options *)NULL);

  if (status != VIO_OK)
    print_error_and_line_num("problems writing  volume `%s'.",__FILE__, __LINE__, filename);

}


static int compare(VIO_Real  *a, VIO_Real *b)
{
  if ( *a > *b)  return(1);
  else if ( *a < *b) return (-1);
  else  return(0);
}

void qs_list(float *item2, int left, int right)
{
    register int i,j;
    float x,y;
    
    i=left;
    j=right;
    x=item2[(left+right)/2];
    
    do
    {
        while (item2[i]<x && i<right) i++;
        while (x<item2[j] && j>left) j--;
        
        if (i<=j)
        {
            y=item2[i];
            item2[i]=item2[j];
            item2[j]=y;
            i++;
            j--;
        }
    } while (i<=j);
    
    if (left<j) qs_list(item2,left,j);
    if (i<right) qs_list(item2,i,right);
}


void normalize_data_to_match_target(VIO_Volume d1, VIO_Volume m1, VIO_Real thresh1,
                                           VIO_Volume d2, VIO_Volume m2, VIO_Real thresh2,
                                           Arg_Data *globals)
{

  VectorR
    vector_step;

  PointR
    starting_position,
    slice,
    row,
    col,
    pos2,
    voxel;

  double
    tx,ty,tz;
  int
    i,j,k,
    r,c,s;

  VIO_Real
    min_range, max_range,
    data_vox, data_val,
    value1, value2;
  
  VIO_Real
    t1,t2,                        /* temporary threshold values     */
    s1,s2,s3;                   /* to store the sums for f1,f2,f3 */
  float 
    *ratios,
    result;                                /* the result */
  int 
    sizes[VIO_MAX_DIMENSIONS],ratios_size,count1,count2;

  VIO_Volume 
    vol;

  VIO_progress_struct
    progress;

  VIO_Data_types 
    data_type;


  set_feature_value_threshold(d1,d2, 
                              &thresh1, &thresh2,
                              &t1,      &t2);                              

  if (globals->flags.debug) {
    print ("In normalize_data_to_match_target, thresh = %10.3f %10.3f\n",t1,t2) ;
  }

  ratios_size = globals->count[ROW_IND] * globals->count[COL_IND] * globals->count[SLICE_IND];

  ALLOC(ratios, ratios_size);  

  fill_Point( starting_position, globals->start[VIO_X], globals->start[VIO_Y], globals->start[VIO_Z]);

  s1 = s2 = s3 = 0.0;
  count1 = count2 = 0;

  for(s=0; s<=globals->count[SLICE_IND]; s++) {

    SCALE_VECTOR( vector_step, globals->directions[SLICE_IND], s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    for(r=0; r<=globals->count[ROW_IND]; r++) {
      
      SCALE_VECTOR( vector_step, globals->directions[ROW_IND], r);
      ADD_POINT_VECTOR( row, slice, vector_step );
      
      SCALE_POINT( col, row, 1.0); /* init first col position */
      for(c=0; c<=globals->count[COL_IND]; c++) {
        
        convert_3D_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);
        
        fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
        
        if (point_not_masked(m1, Point_x(col), Point_y(col), Point_z(col))) {

          value1 = get_value_of_point_in_volume( Point_x(col), Point_y(col), Point_z(col), d1);

          if ( value1 > t1 ) {

            count1++;

            DO_TRANSFORM(pos2, globals->trans_info.transformation, col);
            
            convert_3D_world_to_voxel(d2, Point_x(pos2), Point_y(pos2), Point_z(pos2), &tx, &ty, &tz);
            
            fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
        
            if (point_not_masked(m2, Point_x(pos2), Point_y(pos2), Point_z(pos2))) {

              value2 = get_value_of_point_in_volume( Point_x(pos2), Point_y(pos2), Point_z(pos2), d2);

              if ( (value2 > t2)  && 
                   ((value2 < -1e-15) || (value2 > 1e-15)) ) {
                  
                ratios[count2++] = value1 / value2 ;

                s1 += value1*value2;
                s2 += value1*value1;
                s3 += value2*value2;
                  
                
              } /* if voxel in d2 */
            } /* if point in mask volume two */
          } /* if voxel in d1 */
        } /* if point in mask volume one */
        
        ADD_POINT_VECTOR( col, col, globals->directions[COL_IND] );
        
      } /* for c */
    } /* for r */
  } /* for s */
  

  if (count2 > 0) {

    if (globals->flags.debug) (void)print ("Starting qsort of ratios...");


    qs_list (ratios,0,count2);

    if (globals->flags.debug) (void)print ("Done.\n");

    result = ratios[ (int)(count2/2) ];        /* the median value */

    if (globals->flags.debug) (void)print ("Normalization: %7d %7d -> %10.8f\n",count1,count2,result);

    if ( fabs(result) < 1e-15) {
      print_error_and_line_num("Error computing normalization ratio `%f'.",__FILE__, __LINE__, result);
    }
    else {

      data_type = get_volume_data_type(d1);



      
      switch( data_type ) {
      case VIO_SIGNED_BYTE: 
      case VIO_UNSIGNED_BYTE: 
      case VIO_SIGNED_SHORT: 
      case VIO_UNSIGNED_SHORT: 
	
	/* build temporary working volume */
	
	vol = copy_volume_definition_no_alloc(d1, NC_UNSPECIFIED, FALSE, 0.0, 0.0);
	get_volume_minimum_maximum_real_value(d1, &min_range, &max_range);
	min_range /= result;
	max_range /= result;
	set_volume_real_range(vol, min_range, max_range);
	get_volume_sizes(d1, sizes);
	
	initialize_progress_report(&progress, FALSE, sizes[0]*sizes[1]*sizes[2] + 1,
				   "Normalizing source data" );
	count1 = 0;
	
	/* reset values in the data volume */
	
	for(i=0; i<sizes[0]; i++)
	  for(j=0; j<sizes[1]; j++) {
	    count1++;
	    update_progress_report( &progress, count1);
	    for(k=0; k<sizes[2]; k++) {
	      GET_VOXEL_3D( data_vox,  d1, i, j, k );
	      data_val = CONVERT_VOXEL_TO_VALUE(d1, data_vox);
	      data_val /= result;
	      data_vox = CONVERT_VALUE_TO_VOXEL( vol, data_val);
	      SET_VOXEL_3D( d1 , i, j, k, data_vox );
	    }
	  }
	
	terminate_progress_report( &progress );
	
	set_volume_real_range(d1, min_range, max_range);
	
	if (globals->flags.debug) (void)print ("After normalization min,max, thresh = %f %f %f\n",
					       min_range, max_range, t1/result);
	
	delete_volume(vol);
	break;
	
      default:			/* then volume should be either float or double */
	
	get_volume_sizes(d1, sizes);
	initialize_progress_report(&progress, FALSE, sizes[0]*sizes[1]*sizes[2] + 1,
				   "Normalizing source data" );
	count1 = 0;
	
	/* nomalize the values in the data volume */
	
	for(i=0; i<sizes[0]; i++)
	  for(j=0; j<sizes[1]; j++) {
	    count1++;
	    update_progress_report( &progress, count1);
	    
	    for(k=0; k<sizes[2]; k++) {	/* it should be possible to directly stream through the voxels, without indexing... */
	      GET_VOXEL_3D( data_vox,  d1, i, j, k );
	      data_val = CONVERT_VOXEL_TO_VALUE(d1, data_vox);
	      data_val /= result;
	      data_vox = CONVERT_VALUE_TO_VOXEL( d1, data_val);
	      SET_VOXEL_3D( d1 , i, j, k, data_vox );
	    }
	  }
	terminate_progress_report( &progress );
      }
    }
    
  }
  FREE(ratios);

  
}

