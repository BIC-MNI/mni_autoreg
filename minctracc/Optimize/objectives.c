/* ----------------------------- MNI Header -----------------------------------
@NAME       : objectives.c
@DESCRIPTION: File containing objective function routines.
@METHOD     : 
@GLOBALS    : 
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

@CREATED    : Wed Jun  9 12:56:08 EST 1993 LC
@MODIFIED   :  $Log: objectives.c,v $
@MODIFIED   :  Revision 1.8  1995-02-22 08:56:06  louis
@MODIFIED   :  Montreal Neurological Institute version.
@MODIFIED   :  compiled and working on SGI.  this is before any changes for SPARC/
@MODIFIED   :  Solaris.
@MODIFIED   :
 * Revision 1.7  94/04/06  11:48:44  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.6  94/02/21  16:35:56  louis
 * version before feb 22 changes
 * 
 * Revision 1.5  93/11/15  16:27:07  louis
 * working version, with new library, with RCS revision stuff,
 * before deformations included
 * 
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Optimize/objectives.c,v 1.8 1995-02-22 08:56:06 louis Exp $";
#endif


#include <volume_io.h>
#include <limits.h>
#include <recipes.h>
#include "constants.h"
#include "arg_data.h"

#include "segment_table.h"
#include "local_macros.h"
#include <print_error.h>

extern Arg_Data main_args;

extern Segment_Table *segment_table;

/* ----------------------------- MNI Header -----------------------------------
@NAME       : xcorr_objective
@INPUT      : volumetric data, for use in correlation.  

	      Correlation is accomplished under the currently defined
	      affine transformation in `globals->trans_info.transformation.trans_data'

	      The cross correlation is calculated using sub-sampling
	      in each of the volumes.

@OUTPUT     : none

@RETURNS    : the normalized cross correlation value (min=0.0 max = 1.0)
              0.0 means perfect fit!

@DESCRIPTION: this routine does volumetric subsampling in world space.
              These world coordinates are mapped back into each 
	      data volume to get the actual value at the given location.
	      The globals lattice info is used to calculate
              positions of the sub-samples in each vol.

	      The correlation is equal to:

	      r =    f1 / ( sqrt(f2) * sqrt(f3) )
              
	      where:          
	             f1 = sum( d1 .* d2)    (point to point multiply)
		     f2 = sum( d1 .^ 2 )    (sum of all points squared)
		     f3 = sum( d2 .^ 2 )    (sum of all points squared)

		    note: -point refers to the value at a sub-sample
		          -the xmat is applied to d2, before calculating r
		              
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb 4, 1992 lc original in fit_vol
@MODIFIED   : Wed Jun 16 10:08:46 EST 1993 LC
        new code for minc tracc, copied from routines in fit_vol.
---------------------------------------------------------------------------- */
public float xcorr_objective(Volume d1,
			     Volume d2,
			     Volume m1,
			     Volume m2, 
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
    r,c,s;

  Real
    value1, value2;
  
  Real
    s1,s2,s3;                   /* to store the sums for f1,f2,f3 */
  float 
    result;				/* the result */
  int 
    count1,count2;


  fill_Point( starting_position, globals->start[X], globals->start[Y], globals->start[Z]);

  s1 = s2 = s3 = 0.0;
  count1 = count2 = 0;

  for_inclusive(s,0,globals->count[SLICE_IND]) {

    SCALE_VECTOR( vector_step, globals->directions[SLICE_IND], s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    for_inclusive(r,0,globals->count[ROW_IND]) {
      
      SCALE_VECTOR( vector_step, globals->directions[ROW_IND], r);
      ADD_POINT_VECTOR( row, slice, vector_step );
      
      SCALE_POINT( col, row, 1.0); /* init first col position */
      for_inclusive(c,0,globals->count[COL_IND]) {
	
	convert_3D_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);
	
	fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
	
	if (point_not_masked(m1, Point_x(col), Point_y(col), Point_z(col))) {
	  
	  if (INTERPOLATE_TRUE_VALUE( d1, &voxel, &value1 )) {

	    count1++;

	    DO_TRANSFORM(pos2, globals->trans_info.transformation, col);
	    
	    convert_3D_world_to_voxel(d2, Point_x(pos2), Point_y(pos2), Point_z(pos2), &tx, &ty, &tz);
	    
	    fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
	
	    if (point_not_masked(m2, Point_x(pos2), Point_y(pos2), Point_z(pos2))) {
	      
	      if (INTERPOLATE_TRUE_VALUE( d2, &voxel, &value2 )) {

		count2++;

		if (value1 > globals->threshold[0] || value2 > globals->threshold[1] ) {
		  
		  s1 += value1*value2;
		  s2 += value1*value1;
		  s3 += value2*value2;
		  
		} 
		
	      } /* if voxel in d2 */
	    } /* if point in mask volume two */
	  } /* if voxel in d1 */
	} /* if point in mask volume one */
	
	ADD_POINT_VECTOR( col, col, globals->directions[COL_IND] );
	
      } /* for c */
    } /* for r */
  } /* for s */
  
  result = 1.0 - s1 / (sqrt((double)s2)*sqrt((double)s3));
  
  if (globals->flags.debug) (void)print ("%7d %7d -> %10.8f\n",count1,count2,result);
  
  return (result);
  
}


public float ssc_objective(Volume d1,
			   Volume d2,
			   Volume m1,
			   Volume m2, 
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
    r,c,s;

  Real
    value1, value2;
  
  float 
    result;				/* the result */
  int 
    count1, count2,
    greater;
  unsigned  long
    zero_crossings;



  fill_Point( starting_position, globals->start[X], globals->start[Y], globals->start[Z]);

  greater = TRUE;
  zero_crossings = count1 = count2 = 0;


  /* ------------------------  count along rows (fastest=col) first ------------------- */


  for_inclusive(s,0,globals->count[SLICE_IND]) {

    SCALE_VECTOR( vector_step, globals->directions[SLICE_IND], s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    for_inclusive(r,0,globals->count[ROW_IND]) {
      
      SCALE_VECTOR( vector_step, globals->directions[ROW_IND], r);
      ADD_POINT_VECTOR( row, slice, vector_step );
      
      SCALE_POINT( col, row, 1.0); /* init first col position */
      for_inclusive(c,0,globals->count[COL_IND]) {
	
	convert_3D_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);
	
	fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
	
	if (point_not_masked(m1, Point_x(col), Point_y(col), Point_z(col))) {
	  
	  if (INTERPOLATE_TRUE_VALUE( d1, &voxel, &value1 )) {

	    count1++;

	    DO_TRANSFORM(pos2, globals->trans_info.transformation, col);
	    
	    convert_3D_world_to_voxel(d2, Point_x(pos2), Point_y(pos2), Point_z(pos2), &tx, &ty, &tz);
	    
	    fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
	
	    if (point_not_masked(m2, Point_x(pos2), Point_y(pos2), Point_z(pos2))) {
	      
	      if (INTERPOLATE_TRUE_VALUE( d2, &voxel, &value2 )) {

		count2++;

		if (!((greater && value1>value2) || (!greater && value1<value2))) {
		  greater = !greater;
		  zero_crossings++;
		} 
		
	      } /* if voxel in d2 */
	    } /* if point in mask volume two */
	  } /* if voxel in d1 */
	} /* if point in mask volume one */
	
	ADD_POINT_VECTOR( col, col, globals->directions[COL_IND] );
	
      } /* for c */
    } /* for r */
  } /* for s */

  /* ------------------------  count along cols second  --(fastest=row)--------------- */

  for_inclusive(s,0,globals->count[SLICE_IND]) {

    SCALE_VECTOR( vector_step, globals->directions[SLICE_IND], s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    for_inclusive(c,0,globals->count[COL_IND]) {
      
      SCALE_VECTOR( vector_step, globals->directions[COL_IND], c);
      ADD_POINT_VECTOR( col, slice, vector_step );
      
      SCALE_POINT( row, col, 1.0); /* init first row position */

      for_inclusive(r,0,globals->count[ROW_IND]) {

	convert_3D_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);
	
	fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
	
	if (point_not_masked(m1, Point_x(col), Point_y(col), Point_z(col))) {
	  
	  if (INTERPOLATE_TRUE_VALUE( d1, &voxel, &value1 )) {

	    count1++;

	    DO_TRANSFORM(pos2, globals->trans_info.transformation, col);
	    
	    convert_3D_world_to_voxel(d2, Point_x(pos2), Point_y(pos2), Point_z(pos2), &tx, &ty, &tz);
	    
	    fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
	
	    if (point_not_masked(m2,Point_x(pos2), Point_y(pos2), Point_z(pos2))) {
	      
	      if (INTERPOLATE_TRUE_VALUE( d2, &voxel, &value2 )) {

		count2++;

		if (!((greater && value1>value2) || (!greater && value1<value2))) {
		  greater = !greater;
		  zero_crossings++;
		} 
		
	      } /* if voxel in d2 */
	    } /* if point in mask volume two */
	  } /* if voxel in d1 */
	} /* if point in mask volume one */
	
	ADD_POINT_VECTOR( row, row, globals->directions[ROW_IND] );
	
      } /* for c */
    } /* for r */
  } /* for s */




  /* ------------------------  count along slices last ------------------------ */


  for_inclusive(c,0,globals->count[COL_IND]) {
    
    SCALE_VECTOR( vector_step, globals->directions[COL_IND], c);
    ADD_POINT_VECTOR( col, starting_position, vector_step );

    for_inclusive(r,0,globals->count[ROW_IND]) {
      
      SCALE_VECTOR( vector_step, globals->directions[ROW_IND], r);
      ADD_POINT_VECTOR( row, col, vector_step );
      
      SCALE_POINT( slice, row, 1.0); /* init first col position */

      for_inclusive(s,0,globals->count[SLICE_IND]) {
	
	convert_3D_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);
	
	fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
	
	if (point_not_masked(m1, Point_x(col), Point_y(col), Point_z(col))) {
	  
	  if (INTERPOLATE_TRUE_VALUE( d1, &voxel, &value1 )) {

	    count1++;

	    DO_TRANSFORM(pos2, globals->trans_info.transformation, col);
	    
	    convert_3D_world_to_voxel(d2, Point_x(pos2), Point_y(pos2), Point_z(pos2), &tx, &ty, &tz);
	    
	    fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
	
	    if (point_not_masked(m2, Point_x(pos2), Point_y(pos2), Point_z(pos2))) {
	      
	      if (INTERPOLATE_TRUE_VALUE( d2, &voxel, &value2 )) {

		count2++;

		if (!((greater && value1>value2) || (!greater && value1<value2))) {
		  greater = !greater;
		  zero_crossings++;
		} 
		
	      } /* if voxel in d2 */
	    } /* if point in mask volume two */
	  } /* if voxel in d1 */
	} /* if point in mask volume one */
	
	ADD_POINT_VECTOR( slice, slice, globals->directions[SLICE_IND] );
	
      } /* for c */
    } /* for r */
  } /* for s */


  result = -1.0 * (float)zero_crossings;

  if (globals->flags.debug) (void)print ("%7d %7d -> %10.8f\n",count1,count2,result);

  return (result);
  
}

public float zscore_objective(Volume d1,
			   Volume d2,
			   Volume m1,
			   Volume m2, 
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
    r,c,s;

  Real
    value1, value2;
  
  Real
    z2_sum;
  float 
    result;				/* the result */
  int 
    count1,count2,count3;



  fill_Point( starting_position, globals->start[X], globals->start[Y], globals->start[Z]);

  z2_sum = 0.0;
  count1 = count2 = count3 = 0;

  for_inclusive(s,0,globals->count[SLICE_IND]) {

    SCALE_VECTOR( vector_step, globals->directions[SLICE_IND], s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    for_inclusive(r,0,globals->count[ROW_IND]) {
      
      SCALE_VECTOR( vector_step, globals->directions[ROW_IND], r);
      ADD_POINT_VECTOR( row, slice, vector_step );
      
      SCALE_POINT( col, row, 1.0); /* init first col position */
      for_inclusive(c,0,globals->count[COL_IND]) {
	
	convert_3D_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);
	
	fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
	
	if (point_not_masked(m1, Point_x(col), Point_y(col), Point_z(col))) {
	  
	  if (INTERPOLATE_TRUE_VALUE( d1, &voxel, &value1 )) {

	    count1++;

	    DO_TRANSFORM(pos2, globals->trans_info.transformation, col);
	    
	    convert_3D_world_to_voxel(d2, Point_x(pos2), Point_y(pos2), Point_z(pos2), &tx, &ty, &tz);
	    
	    fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
	
	    if (point_not_masked(m2, Point_x(pos2), Point_y(pos2), Point_z(pos2))) {
	      
	      if (INTERPOLATE_TRUE_VALUE( d2, &voxel, &value2 )) {

		count2++;

		if (ABS(value1) > globals->threshold[0] || ABS(value2) > globals->threshold[1] ) {
		  count3++;
		  z2_sum +=  (value1-value2)*(value1-value2);
		} 
	
		
	      } /* if voxel in d2 */
	    } /* if point in mask volume two */
	  } /* if voxel in d1 */
	} /* if point in mask volume one */
	
	ADD_POINT_VECTOR( col, col, globals->directions[COL_IND] );
	
      } /* for c */
    } /* for r */
  } /* for s */

  if (count3 > 0)
    result = sqrt((double)z2_sum) / count3;
  else
    result = sqrt((double)z2_sum);

  if (globals->flags.debug) (void)print ("%7d %7d %7d -> %10.8f\n",count1,count2,count3,result);
  
  return (result);
  
}



public float vr_objective(Volume d1,
			  Volume d2,
			  Volume m1,
			  Volume m2, 
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
    r,c,s;

  Real
    value1, value2,
    voxel_value1;
  
  float
    rat,
    total_variance,
    *limits,
    *rat_sum,
    *rat2_sum,
    *var;
  unsigned long
    total_count,
    *count3;

  float 
    result;				/* the result */
  int 
    index,i,count1,count2;


				/* build segmentation info!  */

  rat_sum  = vector(1,segment_table->groups);
  rat2_sum = vector(1,segment_table->groups);
  var      = vector(1,segment_table->groups);
  limits   = vector(1,segment_table->groups);

  ALLOC(count3, (segment_table->groups+1));

				/* build world lattice info */


  fill_Point( starting_position, globals->start[X], globals->start[Y], globals->start[Z]);

				/* init running sums and counters. */
  for_inclusive( i, 1, segment_table->groups) {
    rat_sum[i] = 0.0;
    rat2_sum[i] = 0.0;
    count3[i]   = 0;
    var[i] = 0.0;
  }
  count1 = count2 = 0;

				/* loop through each node of lattice */
  for_inclusive(s,0,globals->count[SLICE_IND]) {

    SCALE_VECTOR( vector_step, globals->directions[SLICE_IND], s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    for_inclusive(r,0,globals->count[ROW_IND]) {
      
      SCALE_VECTOR( vector_step, globals->directions[ROW_IND], r);
      ADD_POINT_VECTOR( row, slice, vector_step );
      
      SCALE_POINT( col, row, 1.0); /* init first col position */
      for_inclusive(c,0,globals->count[COL_IND]) {
	
	convert_3D_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);
	
	fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
	
	if (point_not_masked(m1, Point_x(col), Point_y(col), Point_z(col))) {
	  
	  if (INTERPOLATE_TRUE_VALUE( d1, &voxel, &value1 )) {

	    count1++;
	    voxel_value1 = CONVERT_VALUE_TO_VOXEL(d1,value1 );


	    DO_TRANSFORM(pos2, globals->trans_info.transformation, col);
	    
	    convert_3D_world_to_voxel(d2, Point_x(pos2), Point_y(pos2), Point_z(pos2), &tx, &ty, &tz);
	    
	    fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
	
	    if (point_not_masked(m2,Point_x(pos2), Point_y(pos2), Point_z(pos2) )) {
	      
	      if (INTERPOLATE_TRUE_VALUE( d2, &voxel, &value2 )) {

		count2++;
		/* voxel_value2 = CONVERT_VALUE_TO_VOXEL(d1,value2 ); */

	        if (ABS(value1) > globals->threshold[0] && ABS(value2) > globals->threshold[1] ) {

		  index = (*segment_table->segment)( voxel_value1, segment_table);

		  if (index>0) {
		    count3[index]++;
		    rat = value1 / value2;
		    rat_sum[index] += rat;
		    rat2_sum[index] +=  rat*rat;
		  }
		  else {
		    print_error("Cannot segment voxel value %d into one of %d groups.", 
				__FILE__, __LINE__, voxel_value1,segment_table->groups );
		    exit(EXIT_FAILURE);

		  }
		} 
		
		
	      } /* if voxel in d2 */
	    } /* if point in mask volume two */
	  } /* if voxel in d1 */
	} /* if point in mask volume one */
	
	ADD_POINT_VECTOR( col, col, globals->directions[COL_IND] );
	
      } /* for c */
    } /* for r */
  } /* for s */


  total_variance = 0.0;
  total_count = 0;

  for_inclusive( index, 1, segment_table->groups) {
    if (count3[index] > 1) 
      total_count += count3[index];
  }

  if (total_count > 1) {
    for_inclusive( index, 1, segment_table->groups) {
      if (count3[index] > 1) {
	var[index]  = ((double)count3[index]*rat2_sum[index] - rat_sum[index]*rat_sum[index]) / 
	  ((double)count3[index]*((double)count3[index]-1.0));
	
	total_variance += ((double)count3[index]/(double)total_count) * var[index];
      }
      else
	var[index] = 0.0;
    }
  }
  else
    total_variance = FLT_MAX;

  result = total_variance;

  if (globals->flags.debug) print ("%7d %7d %7d -> %10.8f\n",count1,count2,count3[1],result);

  free_vector(rat_sum,  1, segment_table->groups);
  free_vector(rat2_sum, 1, segment_table->groups);
  free_vector(var,      1, segment_table->groups);
  free_vector(limits,   1, segment_table->groups);

  FREE(count3);



  return (result);
  
}

