/* ----------------------------- MNI Header -----------------------------------
@NAME       : obj_fn_mutual_info.c

@DESCRIPTION: collection of routines necessary to calculate the
              mutual information similarity criteria as described by 
	      Collignon et al (IPMI 95) and by Maes (submitted 96)
@METHOD     : the `trick' in the method is the use of `partial volume
              interpolation' (implemented for trilinear interpolation
	      only here).  See descriptions below.
              
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

@CREATED    : Tue Mar 12 09:37:44 MET 1996
@MODIFIED   : $Log: obj_fn_mutual_info.c,v $
@MODIFIED   : Revision 96.2  1997-11-12 21:07:43  louis
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
 * Revision 9.6  1996/08/21  18:22:04  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:16:03  louis
 * Never released version 0.95
 *
 * Revision 1.2  1996/08/12  14:15:56  louis
 * Pre-release
 *
 * Revision 1.1  1996/03/25  10:33:15  collins
 * Initial revision
 *
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Optimize/obj_fn_mutual_info.c,v 96.2 1997-11-12 21:07:43 louis Exp $";
#endif

#include <internal_volume_io.h>
#include "constants.h"
#include "local_macros.h"
#include "arg_data.h"
#include "vox_space.h"

extern Arg_Data main_args;

			/* these are defined/alloc'd in optimize.c  */

extern Real            **prob_hash_table; 
extern Real            *prob_fn1;         
extern Real            *prob_fn2;         

public int point_not_masked(Volume volume, Real wx, Real wy, Real wz);
public int voxel_point_not_masked(Volume volume, 
                                  Real vx, Real vy, Real vz);

				

/* ----------------------------- MNI Header -----------------------------------
@NAME       : partial_volume_interpolation
@INPUT      : volume           - pointer to volume data
              coord[]          - voxel-coordinates of point to interpolate
@OUTPUT     : intensity_vals[] - array of 8 voxel values from corners of
                                 interpolation cube
              intensity_vals[] - array of 8 fractional used used to interpolate
	                         the desired value.
              result           - the interpolated TRUE value.
@RETURNS    : TRUE if wx, wy,wz is within the volume and can be interpolated, 
              FALSE otherwise.
@DESCRIPTION: procedure to compute the partial volume interpolation required
              to evaluate the mutual information objective function.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Tue Mar 12 09:37:44 MET 1996
@MODIFIED   : 
---------------------------------------------------------------------------- */
public BOOLEAN partial_volume_interpolation(Volume data,
					    Real coord[],
					    Real intensity_vals[],
					    Real fractional_vals[],
					    Real *result)
{
  long ind0, ind1, ind2, max[3];
  int sizes[3];
  static double f0, f1, f2, r0, r1, r2, r1r2, r1f2, f1r2, f1f2;
  
  /* Check that the coordinate is inside the volume */
  
  get_volume_sizes(data, sizes);
  max[0]=sizes[0];
  max[1]=sizes[1];
  max[2]=sizes[2];
  
  if (( coord[X]  < 0) || ( coord[X]  >= max[0]-1) ||
      ( coord[Y]  < 0) || ( coord[Y]  >= max[1]-1) ||
      ( coord[Z]  < 0) || ( coord[Z]  >= max[2]-1)) {
    
    return(FALSE);
  }
    
  /* Get the whole part of the coordinate */ 
  ind0 = (long)  coord[X] ;
  ind1 = (long)  coord[Y] ;
  ind2 = (long)  coord[Z] ;
  if (ind0 >= max[0]-1) ind0 = max[0]-1;
  if (ind1 >= max[1]-1) ind1 = max[1]-1;
  if (ind2 >= max[2]-1) ind2 = max[2]-1;
  
  /* Get the relevant voxels */
  GET_VOXEL_3D( intensity_vals[0] ,  data, ind0  , ind1  , ind2   ); 
  GET_VOXEL_3D( intensity_vals[1] ,  data, ind0  , ind1  , ind2+1 ); 
  GET_VOXEL_3D( intensity_vals[2] ,  data, ind0  , ind1+1, ind2   ); 
  GET_VOXEL_3D( intensity_vals[3] ,  data, ind0  , ind1+1, ind2+1 ); 
  GET_VOXEL_3D( intensity_vals[4] ,  data, ind0+1, ind1  , ind2   ); 
  GET_VOXEL_3D( intensity_vals[5] ,  data, ind0+1, ind1  , ind2+1 ); 
  GET_VOXEL_3D( intensity_vals[6] ,  data, ind0+1, ind1+1, ind2   ); 
  GET_VOXEL_3D( intensity_vals[7] ,  data, ind0+1, ind1+1, ind2+1 ); 

  /* Get the fraction parts */
  f0 =  coord[X]  - ind0;
  f1 =  coord[Y]  - ind1;
  f2 =  coord[Z]  - ind2;
  r0 = 1.0 - f0;
  r1 = 1.0 - f1;
  r2 = 1.0 - f2;
  
  r1r2 = r1 * r2;
  r1f2 = r1 * f2;
  f1r2 = f1 * r2;
  f1f2 = f1 * f2;

  fractional_vals[0] = r0 * r1r2;
  fractional_vals[1] = r0 * r1f2;
  fractional_vals[2] = r0 * f1r2;
  fractional_vals[3] = r0 * f1f2;
  fractional_vals[4] = f0 * r1r2;
  fractional_vals[5] = f0 * r1f2;
  fractional_vals[6] = f0 * f1r2;
  fractional_vals[7] = f0 * f1f2;

  /* Do the interpolation */

  *result =
    r0 *  (r1r2 * intensity_vals[0] +
	   r1f2 * intensity_vals[1] +
	   f1r2 * intensity_vals[2] +
	   f1f2 * intensity_vals[3]);
  *result +=
    f0 *  (r1r2 * intensity_vals[4] +
	   r1f2 * intensity_vals[5] +
	   f1r2 * intensity_vals[6] +
	   f1f2 * intensity_vals[7]);

  *result = CONVERT_VOXEL_TO_VALUE(data, *result);

  return TRUE;
  
}
	

private blur_pdf( Real *pdf, int blur_size, int pdf_length) {

    Real *temp_pdf;
    int  i,j,blur_by2;

    if (blur_size > 1) 
    {
	
        ALLOC(temp_pdf, pdf_length);

				/* copy the pdf into the temporary pdf */
	for_less(i,0,pdf_length)	   
	    temp_pdf[i] = pdf[i];

	if (blur_size==3) 
	{
	   for_less(i,1,pdf_length-1)  
	       pdf[i] = (temp_pdf[i-1] + temp_pdf[i] + temp_pdf[i+1])/3.0;	   
	}
	else 
	{
	    
				/* now blur it, storing result in pdf */

	    blur_by2 = (int)(blur_size / 2);
	    
	    blur_size = blur_by2 * 2 + 1;
	    
				/* blur the starting end */
	    if (blur_by2 > 1 && pdf_length > blur_size)	
	    {
		for_less(i,0,blur_by2) 
		{
		    pdf[i] = 0.0;
		    for_less(j,0,2*i)
			pdf[i] += temp_pdf[j];
		    pdf[i] /= (2.0*i);		
		}	    
	    }
	    
	    
	    for_less(i,blur_by2,pdf_length-blur_by2) 
	    {
		pdf[i] = 0.0;
		
		for_inclusive(j, -blur_by2, blur_by2)
		    pdf[i] += temp_pdf[i+j];
		
		pdf[i] /= (Real)blur_size;	    
	    }	
			   
				/* blur the ending end */
	    if (blur_by2 > 1 && pdf_length > blur_size)	
	    {
		for_less(i,0,blur_by2) 
		{
		    pdf[pdf_length-i-1] = 0.0;
		    for_less(j,0,2*i)
			pdf[pdf_length-i-1] += temp_pdf[pdf_length-1-j];
		    pdf[pdf_length-i-1] /= (2.0*i);		
		}	    
	    }
	    
	}	
	FREE(temp_pdf);
    }
    
}    



private blur_jpdf (Real **hist, int blur_size, int pdf_length) {

    Real **temp_hist, *temp_col;
    int i,j;
    
    if (blur_size > 1) 
    {
		    
	ALLOC2D(temp_hist, pdf_length, pdf_length);
	ALLOC  (temp_col, pdf_length);

	for_less(i,0,pdf_length)   /* copy the histogram */
	    for_less(j,0,pdf_length)	   
		temp_hist[i][j] = hist[i][j];

	for_less(i,1,pdf_length-1)   /* blur the rows */
	    blur_pdf(temp_hist[i], blur_size, pdf_length);
	
	for_less(i,1,pdf_length-1)  { /* blur the cols */

	    for_less(j,0,pdf_length)
		temp_col[j] = temp_hist[j][i];	    

	    blur_pdf(temp_col, blur_size, pdf_length);

	    for_less(j,0,pdf_length)
		hist[j][i] = temp_col[j];	    
	}
	

	FREE  (temp_col);
	FREE2D(temp_hist);
    }

}


/* this function will calculate the mutual information similarity
   value based on the paper by Collignon, IPMI95, p 266 

   limits/constraints/caveats:
   - this version is limited to manipulation of byte volumes only
   - globals->groups must = 256 (since only byte data supported)
     (this is forced in optimize.c)
   - all calculations are computed on voxel (byte) intensity values
     (and _not_ the REAL value, as is done with all other obj functions)
   - because of this, the input data vols should be appropriately set up
     by the user
   - ONLY partial volume interpolation is used: there is no support for
     other interpolation methods.

*/

public float mutual_information_objective(Volume d1,
					  Volume d2,
					  Volume m1,
					  Volume m2, 
					  Arg_Data *globals)
{

  VectorR			/* these variables are used to step through */
    vector_step;		/* the 3D lattice                           */
  PointR
    starting_position,
    slice,
    row,
    col,
    pos2;
  Real
    voxel_coord[3];
  int
    i,j,
    count1,count2,		/* number of nodes in first vol, second vol */
    index1[8],
    index2[8],
    r,c,s;
  
  Real
    *temp_pdf,
    **temp_hist,
    min_range1, max_range1, range1,
    min_range2, max_range2, range2,
    intensity_vals1[8],		/* voxel values to index into histogram */
    intensity_vals2[8],
    fractional_vals1[8],	/* fractional values to add to histo */
    fractional_vals2[8],
    value1, value2;
  float 
    mutual_info_result;			

  Voxel_space_struct *vox_space;
  Transform          *trans;


				/* init any objective function specific
				   stuff here                           */
  count1 = count2 = 0;
  mutual_info_result = 0.0;

  for_less(i,0,globals->groups) {
    prob_fn1[i] = 0.0;
    prob_fn2[i] = 0.0;
  }

  for_less(i,0,globals->groups) 
    for_less(j,0,globals->groups) 
      prob_hash_table[i][j] = 0.0;

  get_volume_real_range(d1, &min_range1, &max_range1);
  get_volume_real_range(d2, &min_range2, &max_range2);

  range1 = max_range1 - min_range1;
  range2 = max_range2 - min_range2;

                                /* prepare data for the voxel-to-voxel
                                   space transformation (instead of the
                                   general but inefficient world-world
                                   computations. */

  vox_space = new_voxel_space_struct();
  get_into_voxel_space(globals, vox_space, d1, d2);
  trans = get_linear_transform_ptr(vox_space->voxel_to_voxel_space);

				/* get ready to step though the 3D lattice
                                   */

  fill_Point( starting_position, vox_space->start[X], vox_space->start[Y], vox_space->start[Z]);

  /* ---------- step through all slices of lattice ------------- */
  for_less(s,0,globals->count[SLICE_IND]) {

    SCALE_VECTOR( vector_step, vox_space->directions[SLICE_IND], s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    /* ---------- step through all rows of lattice ------------- */
    for_less(r,0,globals->count[ROW_IND]) {
      
      SCALE_VECTOR( vector_step, vox_space->directions[ROW_IND], r);
      ADD_POINT_VECTOR( row, slice, vector_step );
      
      SCALE_POINT( col, row, 1.0); /* init first col position */

      /* ---------- step through all cols of lattice ------------- */
      for_less(c,0,globals->count[COL_IND]) {
	
				   /* get the node value in volume 1,
				      if it falls within the volume    */

	if (voxel_point_not_masked(m1, Point_x(col), Point_y(col), Point_z(col))) {
	  
           voxel_coord[X] = Point_x(col);
           voxel_coord[Y] = Point_y(col);
           voxel_coord[Z] = Point_z(col);
           

	  if (partial_volume_interpolation(d1, 
					   voxel_coord, 
					   intensity_vals1,
					   fractional_vals1,
					   &value1 )) {

	    if (value1 > globals->threshold[0]) { /* is the voxel in the thresholded region? */

	      count1++;
				/* transform the node coordinate into
				   volume 2                             */

              my_homogenous_transform_point(trans,
                                            Point_x(col), Point_y(col), Point_z(col), 1.0,
                                            &Point_x(pos2), &Point_y(pos2), &Point_z(pos2));
	      
	      /* get the node value in volume 2,
		 if it falls within the volume    */
	      
	      if (voxel_point_not_masked(m2,Point_x(pos2), Point_y(pos2), Point_z(pos2) )) {
                 
                 voxel_coord[X] = Point_x(pos2);
                 voxel_coord[Y] = Point_y(pos2);
                 voxel_coord[Z] = Point_z(pos2);
                 
                 if (partial_volume_interpolation(d2, 
                                                  voxel_coord, 
                                                  intensity_vals2,
                                                  fractional_vals2,
                                                  &value2 )) {
		  
                    if (value2 > globals->threshold[1]) { /* is the voxel in the thresholded region? */

                       count2++;
                       
                       for_less(i,0,8) {
                          index1[i] = ROUND( intensity_vals1[i] );
                          index2[i] = ROUND( intensity_vals2[i] );
                          prob_fn1[ index1[i] ] += fractional_vals1[i];
                          prob_fn2[ index2[i] ] += fractional_vals2[i];
                       }
                       for_less(i,0,8) 
                          for_less(j,0,8) {
                             prob_hash_table[ index1[i] ][ index2[j] ] += 
                                fractional_vals1[i]*fractional_vals2[j];
                          }
                       
                       
                    } /* if value2>thres */
                 } /* if voxel in d2 */
	      } /* if point in mask volume two */
           } /* if value1>thres */
         } /* if voxel in d1 */
	} /* if point in mask volume one */
	
	ADD_POINT_VECTOR( col, col, vox_space->directions[COL_IND] );
	
      } /* for c */
    } /* for r */
  } /* for s */



  /* now that the data for the objective function has been accumulated
     over the lattice nodes, blur the probability distribution functions
  */

  blur_pdf (prob_fn1,        globals->blur_pdf, globals->groups);
  blur_pdf (prob_fn2,        globals->blur_pdf, globals->groups);
  blur_jpdf(prob_hash_table, globals->blur_pdf, globals->groups);  

  /* now finish the objective function calculation, 
     placing the final objective function value in  'mutual_info_result' */

      
  mutual_info_result = 0.0;
  if (count2>0) {

				/* normalize to count2  */
    for_less(i,0,globals->groups) {
      prob_fn1[i] /= count2;
      prob_fn2[i] /= count2;
    }
    
    for_less(i,0,globals->groups) 
      for_less(j,0,globals->groups) 
	prob_hash_table[i][j] /= count2;
    

    for_less(i,0,globals->groups) 
      for_less(j,0,globals->groups) {
	
	if ( prob_fn1[i] > 0.0 &&  prob_fn2[j] > 0.0 && prob_hash_table[i][j]>0.0)

	  mutual_info_result += prob_hash_table[i][j] * 
	         log( prob_hash_table[i][j] / (prob_fn1[i] * prob_fn2[j])
		    );
      }
    
    mutual_info_result *= -1.0;
  }

  if (globals->flags.debug) (void)print ("%7d %7d -> %f\n",count1,count2,mutual_info_result);


  /* don't forget to free up any variables you declared above */
    
  return (mutual_info_result);
  
}

