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
@MODIFIED   : Revision 96.11  2009-04-03 18:36:59  louis
@MODIFIED   : made changes to use only DOUBLES for input source and model volumes, and for all estimation of deformation fields
@MODIFIED   :
@MODIFIED   : Revision 96.10  2008/10/08 15:17:49  louis
@MODIFIED   : added -nmi option for linear normalized mutual information
@MODIFIED   :
@MODIFIED   : Revision 96.9  2006/11/30 09:07:32  rotor
@MODIFIED   :  * many more changes for clean minc 2.0 build
@MODIFIED   :
@MODIFIED   : Revision 96.8  2006/11/29 09:09:34  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 96.7  2005/07/20 20:45:50  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 96.6  2004/02/12 06:08:21  rotor
@MODIFIED   :  * removed /static defs
@MODIFIED   :
@MODIFIED   : Revision 96.5  2002/03/26 14:15:44  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 96.4  2000/03/15 08:42:47  stever
@MODIFIED   : Code cleanup: all functions prototyped (except ParseArgs.c), no useless declarations, etc
@MODIFIED   :
@MODIFIED   : Revision 96.3  1999/10/25 19:59:08  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :
 * Revision 96.2  1997/11/12  21:07:43  louis
 * no changes, other than rcsid...
 *
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
static char rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/minctracc/Optimize/obj_fn_mutual_info.c,v 96.11 2009-04-03 18:36:59 louis Exp $";
#endif

#include <volume_io.h>
#include "constants.h"
#include "local_macros.h"
#include "minctracc_arg_data.h"
#include "vox_space.h"
#include "objectives.h"
#include <math.h>

extern Arg_Data *main_args;

                        /* these are defined/alloc'd in optimize.c  */

extern VIO_Real            **prob_hash_table; 
extern VIO_Real            *prob_fn1;         
extern VIO_Real            *prob_fn2;         

int point_not_masked(VIO_Volume volume, VIO_Real wx, VIO_Real wy, VIO_Real wz);
int voxel_point_not_masked(VIO_Volume volume, 
                                  VIO_Real vx, VIO_Real vy, VIO_Real vz);

                                

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
VIO_BOOL partial_volume_interpolation(VIO_Volume data,
                                            VIO_Real coord[],
                                            VIO_Real intensity_vals[],
                                            VIO_Real fractional_vals[],
                                            VIO_Real *result)
{
  long ind0, ind1, ind2, max[3];
  int sizes[3];
  static double f0, f1, f2, r0, r1, r2, r1r2, r1f2, f1r2, f1f2;
  
  /* Check that the coordinate is inside the volume */
  
  get_volume_sizes(data, sizes);
  max[0]=sizes[0];
  max[1]=sizes[1];
  max[2]=sizes[2];
  
  if (( coord[VIO_X]  < 0) || ( coord[VIO_X]  >= max[0]-1) ||
      ( coord[VIO_Y]  < 0) || ( coord[VIO_Y]  >= max[1]-1) ||
      ( coord[VIO_Z]  < 0) || ( coord[VIO_Z]  >= max[2]-1)) {
    
    return(FALSE);
  }
    
  /* Get the whole part of the coordinate */ 
  ind0 = (long)  coord[VIO_X] ;
  ind1 = (long)  coord[VIO_Y] ;
  ind2 = (long)  coord[VIO_Z] ;
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
  f0 =  coord[VIO_X]  - ind0;
  f1 =  coord[VIO_Y]  - ind1;
  f2 =  coord[VIO_Z]  - ind2;
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
        

static void blur_pdf( VIO_Real *pdf, int blur_size, int pdf_length) {

    VIO_Real *temp_pdf;
    int  i,j,blur_by2;

    if (blur_size > 1) 
    {
        
        ALLOC(temp_pdf, pdf_length);

                                /* copy the pdf into the temporary pdf */
        for(i=0; i<pdf_length; i++)           
            temp_pdf[i] = pdf[i];

        if (blur_size==3) 
        {
           for(i=1; i<pdf_length-1; i++)  
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
                for(i=0; i<blur_by2; i++) 
                {
                    pdf[i] = 0.0;
                    for(j=0; j<2*i; j++)
                        pdf[i] += temp_pdf[j];
                    pdf[i] /= (2.0*i);                
                }            
            }
            
            
            for(i=blur_by2; i<pdf_length-blur_by2; i++) 
            {
                pdf[i] = 0.0;
                
                for(j=-blur_by2; j<=blur_by2; j++)
                    pdf[i] += temp_pdf[i+j];
                
                pdf[i] /= (VIO_Real)blur_size;            
            }        
                           
                                /* blur the ending end */
            if (blur_by2 > 1 && pdf_length > blur_size)        
            {
                for(i=0; i<blur_by2; i++) 
                {
                    pdf[pdf_length-i-1] = 0.0;
                    for(j=0; j<2*i; j++)
                        pdf[pdf_length-i-1] += temp_pdf[pdf_length-1-j];
                    pdf[pdf_length-i-1] /= (2.0*i);                
                }            
            }
            
        }        
        FREE(temp_pdf);
    }
    
}    



static void blur_jpdf (VIO_Real **hist, int blur_size, int pdf_length) {

    VIO_Real **temp_hist, *temp_col;
    int i,j;
    
    if (blur_size > 1) 
    {
                    
        VIO_ALLOC2D(temp_hist, pdf_length, pdf_length);
        ALLOC  (temp_col, pdf_length);

        for(i=0; i<pdf_length; i++)   /* copy the histogram */
            for(j=0; j<pdf_length; j++)           
                temp_hist[i][j] = hist[i][j];

        for(i=1; i<pdf_length-1; i++)   /* blur the rows */
            blur_pdf(temp_hist[i], blur_size, pdf_length);
        
        for(i=1; i<pdf_length-1; i++)  { /* blur the cols */

            for(j=0; j<pdf_length; j++)
                temp_col[j] = temp_hist[j][i];            

            blur_pdf(temp_col, blur_size, pdf_length);

            for(j=0; j<pdf_length; j++)
                hist[j][i] = temp_col[j];            
        }
        

        FREE  (temp_col);
        VIO_FREE2D(temp_hist);
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

float mutual_information_objective(VIO_Volume d1,
                                          VIO_Volume d2,
                                          VIO_Volume m1,
                                          VIO_Volume m2, 
                                          Arg_Data *globals)
{

  VectorR                        /* these variables are used to step through */
    vector_step;                /* the 3D lattice                           */
  PointR
    starting_position,
    slice,
    row,
    col,
    pos2;
  VIO_Real
    voxel_coord[3];
  int
    i,j,
    count1,count2,                /* number of nodes in first vol, second vol */
    index1[8],
    index2[8],
    r,c,s;
  
  VIO_Real
    min_range1, max_range1, range1,
    min_range2, max_range2, range2,
    intensity_vals1[8],                /* voxel values to index into histogram */
    intensity_vals2[8],
    fractional_vals1[8],        /* fractional values to add to histo */
    fractional_vals2[8],
    value1, value2;
  double
    Hy, Hx, Ixy;		/* entropies */
  double
    product, Redundancy;
  float 
    mutual_info_result;                        

  Voxel_space_struct *vox_space;
  VIO_Transform          *trans;


                                /* init any objective function specific
                                   stuff here                           */
  count1 = count2 = 0;
  mutual_info_result = 0.0;

  for(i=0; i<globals->groups; i++) {
    prob_fn1[i] = 0.0;
    prob_fn2[i] = 0.0;
  }

  for(i=0; i<globals->groups; i++) 
    for(j=0; j<globals->groups; j++) 
      prob_hash_table[i][j] = 0.0;

  /*
    this was here, but appears to be useless!  dlc 04/2009

    get_volume_minimum_maximum_real_value(d1, &min_range1, &max_range1);
    get_volume_minimum_maximum_real_value(d2, &min_range2, &max_range2);
    
    range1 = max_range1 - min_range1;
    range2 = max_range2 - min_range2;

  */


                                /* prepare data for the voxel-to-voxel
                                   space transformation (instead of the
                                   general but inefficient world-world
                                   computations. */

  vox_space = new_voxel_space_struct();
  get_into_voxel_space(globals, vox_space, d1, d2);
  trans = get_linear_transform_ptr(vox_space->voxel_to_voxel_space);

                                /* get ready to step though the 3D lattice
                                   */

  fill_Point( starting_position, vox_space->start[VIO_X], vox_space->start[VIO_Y], vox_space->start[VIO_Z]);

  /* ---------- step through all slices of lattice ------------- */
  for(s=0; s<globals->count[SLICE_IND]; s++) {

    SCALE_VECTOR( vector_step, vox_space->directions[SLICE_IND], s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    /* ---------- step through all rows of lattice ------------- */
    for(r=0; r<globals->count[ROW_IND]; r++) {
      
      SCALE_VECTOR( vector_step, vox_space->directions[ROW_IND], r);
      ADD_POINT_VECTOR( row, slice, vector_step );
      
      SCALE_POINT( col, row, 1.0); /* init first col position */

      /* ---------- step through all cols of lattice ------------- */
      for(c=0; c<globals->count[COL_IND]; c++) {
        
                                   /* get the node value in volume 1,
                                      if it falls within the volume    */

        if (voxel_point_not_masked(m1, Point_x(col), Point_y(col), Point_z(col))) {
          
           voxel_coord[VIO_X] = Point_x(col);
           voxel_coord[VIO_Y] = Point_y(col);
           voxel_coord[VIO_Z] = Point_z(col);
           

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
                 
                 voxel_coord[VIO_X] = Point_x(pos2);
                 voxel_coord[VIO_Y] = Point_y(pos2);
                 voxel_coord[VIO_Z] = Point_z(pos2);
                 
                 if (partial_volume_interpolation(d2, 
                                                  voxel_coord, 
                                                  intensity_vals2,
                                                  fractional_vals2,
                                                  &value2 )) {
                  
                    if (value2 > globals->threshold[1]) { /* is the voxel in the thresholded region? */

                       count2++;
                       
                       for(i=0; i<8; i++) {
                          index1[i] = VIO_ROUND( intensity_vals1[i] );
                          index2[i] = VIO_ROUND( intensity_vals2[i] );
                          prob_fn1[ index1[i] ] += fractional_vals1[i];
                          prob_fn2[ index2[i] ] += fractional_vals2[i];
                       }
                       for(i=0; i<8; i++) 
                          for(j=0; j<8; j++) {
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
  Hx = 0.0;
  Hy = 0.0;
  Ixy = 0.0;
  


  if (count2>0) {

                                /* normalize to count2  */
    for(i=0; i<globals->groups; i++) {
      prob_fn1[i] /= count2;
      prob_fn2[i] /= count2;
    }
    
    for(i=0; i<globals->groups; i++) 
      for(j=0; j<globals->groups; j++) 
        prob_hash_table[i][j] /= count2;
    


    if ( globals->obj_function == normalized_mutual_information_objective ) { /* ie, -nmi option */

      /* mutual information of X and Y is defined as
	 I(X;Y) = sum_x ( sum_y ( p(x,y) * log[ p(x,y) / ( p1(x)*p2(y) )  ]  )
	 where p(x,y) is joint pobability distribution of X and Y
	 and p1(x) and p2(y) are the marginal probability distibutions
	     I(X;Y) is computed for -mi option (below in the else)

	 note that 
	 I(X;Y) = I(Y,X) -> symmetic
	 I(X,Y) >=0      -> non-negative
     
	 if X and Y are inpendent, then 
	    ->  p(x,y) = p(x)*p(z) and 
	    -> log [ p(x,y) / (p(x)*p(y) ] = log (1) = 0
	 
	 normalized MI = nMI = redundancy:

	 R = I(X;Y) / (H(X) + H(Y))
	 attains a minimum of 0; 
	 and a max of min( H(X),H(Y) ) /  (H(X) + H(Y))

	 where H(x) = - sum_i p(x_i)*log p(x_i)

	 so here, we compute the normalized symmetric redudancy based on MI

	 so we need 
	    H(X) and H(Y) (these are variables Hx and Hy below
	    I(X;Y) (stored in Ixy below)
      */



      for(i=0; i<globals->groups; i++) {	/* compute marginal entropies */
	if (prob_fn1[i]>0.0) Hx += -1.0 * (double)prob_fn1[i] * log((double)prob_fn1[i]);
	if (prob_fn2[i]>0.0) Hy += -1.0 * (double)prob_fn2[i] * log((double)prob_fn2[i]);	
	
      }
      
      for(i=0; i<globals->groups; i++) {        /* compute mutual information */
	for(j=0; j<globals->groups; j++) {
	  product = prob_fn1[i]*prob_fn2[j] ;
	  if (prob_hash_table[i][j]>0.0 && product>0.0) 
	    Ixy += (double)prob_hash_table[i][j] *  log( (double)( prob_hash_table[i][j]/product));
	}
      }
	     
      if (( Hx + Hy) > 0.0) {	/* compute redundancy, and return as the normalized mutual info */
	Redundancy = Ixy / (Hx + Hy);
	mutual_info_result = Redundancy;
      }
      else {
	mutual_info_result = 0.0;
      }

    } else {			/* this is the standard MI computation pre Oct 2008 (the -mi option for linear reg) */
      for(i=0; i<globals->groups; i++) 
	for(j=0; j<globals->groups; j++) {
	  
	  if ( prob_fn1[i] > 0.0 &&  prob_fn2[j] > 0.0 && prob_hash_table[i][j]>0.0)
	    /* this is the same as Ixy, just above */
	    mutual_info_result += prob_hash_table[i][j] * 
	      log( prob_hash_table[i][j] / (prob_fn1[i] * prob_fn2[j]) );
	}
    }

    mutual_info_result *= -1.0;
  }

  if (globals->flags.debug) {
    (void)print ("%7d %7d -> %f ( %f %f %f )\n",count1,count2,mutual_info_result, Hx, Hy, Ixy);
  }


  /* don't forget to free up any variables you declared above */
    
  return (mutual_info_result);
  
}


float normalized_mutual_information_objective(VIO_Volume d1,
                                          VIO_Volume d2,
                                          VIO_Volume m1,
                                          VIO_Volume m2, 
                                          Arg_Data *globals)
{
  return (mutual_information_objective(d1,d2,m1,m2,globals));
}
