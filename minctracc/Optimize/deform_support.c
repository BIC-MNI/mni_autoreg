/* ----------------------------- MNI Header -----------------------------------
@NAME       : deform_support.c
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
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

@CREATED    : Tue Feb 22 08:37:49 EST 1994
@MODIFIED   : $Log: deform_support.c,v $
@MODIFIED   : Revision 1.6  1995-05-02 11:30:27  louis
@MODIFIED   : started clean up of code, separation of non used procedures into
@MODIFIED   : old_methods.c.  This version was working, but I am now going to
@MODIFIED   : rewrite everything to use GRID_TRANSFORM.
@MODIFIED   :
 * Revision 1.5  1995/02/22  08:56:06  louis
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.4  94/06/21  10:59:28  louis
 * working optimized version of program.  when compiled with -O, this
 * code is approximately 60% faster than the previous version.
 * 
 * 
 * Revision 1.3  94/06/06  18:46:53  louis
 * working version: clamp and blur of deformation lattice now ensures
 * a smooth recovered deformation.  Unfortunately, the r = cost-similarity
 * function used in the optimization is too heavy on the cost_fn.  This has
 * to get fixed...
 * 
 * 

 * Revision 1.2  94/06/02  20:15:56  louis
 * made modifications to allow deformations to be calulated in 2D on slices. 
 * changes had to be made in set_up_lattice, init_lattice when defining
 * the special case of a single slice....
 * Build_default_deformation_field also had to reflect these changes.
 * do_non-linear-optimization also had to check if one of dimensions had
 * a single element.
 * All these changes were made, and slightly tested.  Another type of
 * deformation strategy will be necessary (to replace the deformation 
 * perpendicular to the surface, since it does not work well).
 * 
 * Revision 1.1  94/04/06  11:47:27  louis
 * Initial revision
 * 

---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Optimize/deform_support.c,v 1.6 1995-05-02 11:30:27 louis Exp $";
#endif

#include <volume_io.h>
#include <louis_splines.h>

#include <limits.h>
#include "line_data.h"
#include "point_vector.h"



#define DERIV_FRAC      0.6
#define FRAC1           0.5
#define FRAC2           0.0833333

extern double smoothing_weight;


public void set_up_lattice(Volume data, 
			   double *user_step, /* user requested spacing for lattice */
			   double *start,     /* world starting position of lattice */
			   int    *count,     /* number of steps in each direction */
			   double *step,      /* step size in each direction */
			   VectorR directions[]);/* array of vector directions for each index*/


public void get_average_warp_of_neighbours(Volume dx, Volume dy, Volume dz, 
					   int i,int j,int k,
					   Real *mx, Real *my, Real *mz)
{
  int 
    count, sizes[3];
  Real 
    ri,rj,rk,
    voxel,
    val_x, val_y, val_z,
    tx, ty, tz;
  
  get_volume_sizes(dx, sizes);

  *mx = 0.0; *my = 0.0; *mz = 0.0;
  ri = (Real)i; rj = (Real)j; rk = (Real)k;
  count = 0;

  if (sizes[0]>2) {
    convert_3D_voxel_to_world(dx, ri+1, rj, rk, &tx, &ty, &tz);
    GET_VOXEL_3D(voxel, dx, i+1,j,k); val_x = CONVERT_VOXEL_TO_VALUE(dx , voxel ); 
    GET_VOXEL_3D(voxel, dy, i+1,j,k); val_y = CONVERT_VOXEL_TO_VALUE(dy , voxel ); 
    GET_VOXEL_3D(voxel, dz, i+1,j,k); val_z = CONVERT_VOXEL_TO_VALUE(dz , voxel ); 
    *mx += tx+val_x; *my += ty+val_y; *mz+= tz+val_z;
    
    convert_3D_voxel_to_world(dx, ri-1, rj, rk, &tx, &ty, &tz);
    GET_VOXEL_3D(voxel, dx, i-1,j,k); val_x = CONVERT_VOXEL_TO_VALUE(dx , voxel ); 
    GET_VOXEL_3D(voxel, dy, i-1,j,k); val_y = CONVERT_VOXEL_TO_VALUE(dy , voxel ); 
    GET_VOXEL_3D(voxel, dz, i-1,j,k); val_z = CONVERT_VOXEL_TO_VALUE(dz , voxel ); 
    *mx += tx+val_x; *my += ty+val_y; *mz+= tz+val_z;
    count += 2;
  }

  if (sizes[1]>2) {
    convert_3D_voxel_to_world(dy, ri, rj+1, rk, &tx, &ty, &tz);
    GET_VOXEL_3D(voxel, dx, i,j+1,k); val_x = CONVERT_VOXEL_TO_VALUE(dx , voxel ); 
    GET_VOXEL_3D(voxel, dy, i,j+1,k); val_y = CONVERT_VOXEL_TO_VALUE(dy , voxel ); 
    GET_VOXEL_3D(voxel, dz, i,j+1,k); val_z = CONVERT_VOXEL_TO_VALUE(dz , voxel ); 
    *mx += tx+val_x; *my += ty+val_y; *mz+= tz+val_z;
    
    convert_3D_voxel_to_world(dy, ri, rj-1, rk, &tx, &ty, &tz);
    GET_VOXEL_3D(voxel, dx, i,j-1,k); val_x = CONVERT_VOXEL_TO_VALUE(dx , voxel ); 
    GET_VOXEL_3D(voxel, dy, i,j-1,k); val_y = CONVERT_VOXEL_TO_VALUE(dy , voxel ); 
    GET_VOXEL_3D(voxel, dz, i,j-1,k); val_z = CONVERT_VOXEL_TO_VALUE(dz , voxel ); 
    *mx += tx+val_x; *my += ty+val_y; *mz+= tz+val_z;
    count += 2;
  }

  if (sizes[2]>2) {
    convert_3D_voxel_to_world(dz, ri, rj, rk+1, &tx, &ty, &tz);
    GET_VOXEL_3D(voxel, dx, i,j,k+1); val_x = CONVERT_VOXEL_TO_VALUE(dx , voxel ); 
    GET_VOXEL_3D(voxel, dy, i,j,k+1); val_y = CONVERT_VOXEL_TO_VALUE(dy , voxel ); 
    GET_VOXEL_3D(voxel, dz, i,j,k+1); val_z = CONVERT_VOXEL_TO_VALUE(dz , voxel ); 
    *mx += tx+val_x; *my += ty+val_y; *mz+= tz+val_z;
    
    convert_3D_voxel_to_world(dz, ri, rj, rk-1, &tx, &ty, &tz);
    GET_VOXEL_3D(voxel, dx, i,j,k-1); val_x = CONVERT_VOXEL_TO_VALUE(dx , voxel ); 
    GET_VOXEL_3D(voxel, dy, i,j,k-1); val_y = CONVERT_VOXEL_TO_VALUE(dy , voxel ); 
    GET_VOXEL_3D(voxel, dz, i,j,k-1); val_z = CONVERT_VOXEL_TO_VALUE(dz , voxel ); 
    *mx += tx+val_x; *my += ty+val_y; *mz+= tz+val_z;
    count += 2;
  }
  
  *mx /= count; *my /= count; *mz /= count; /* average warp of neighbours, in target space  */
}

				/* add additional to current, return answer in additional */

public void add_additional_warp_to_current(Volume adx, Volume ady, Volume adz,
					   Volume dx, Volume dy, Volume dz,
					   Real weight)
{
  int
    loop_start[3], loop_end[3],
    i,j,k,sizes[3];
  Real 
    value, value2, voxel;

  get_volume_sizes(adx, sizes);

  for_less(i,0,3) {
    if (sizes[i]>3) {
      loop_start[i] = 1;
      loop_end[i] = sizes[i]-1;
    }
    else {
      loop_start[i]=0;
      loop_end[i] = sizes[i];
    }
  }


  for_less(i,loop_start[0],loop_end[0]) {
    for_less(j,loop_start[1],loop_end[1]) {
      for_less(k,loop_start[2],loop_end[2]){
	
	
	GET_VOXEL_3D(voxel, adx, i,j,k); value = CONVERT_VOXEL_TO_VALUE( adx, voxel ); 
	GET_VOXEL_3D(voxel, dx,  i,j,k); value2 = CONVERT_VOXEL_TO_VALUE( dx, voxel ); value2 += value*weight;
	voxel = CONVERT_VALUE_TO_VOXEL(adx, value2);
	SET_VOXEL_3D(adx, i,j,k, voxel); 

	GET_VOXEL_3D(voxel, ady, i,j,k); value = CONVERT_VOXEL_TO_VALUE( ady, voxel ); 
	GET_VOXEL_3D(voxel, dy,  i,j,k); value2 = CONVERT_VOXEL_TO_VALUE( dy, voxel ); value2 += value*weight;
	voxel = CONVERT_VALUE_TO_VOXEL(ady, value2);
	SET_VOXEL_3D(ady, i,j,k, voxel); 
		
	GET_VOXEL_3D(voxel, adz, i,j,k); value = CONVERT_VOXEL_TO_VALUE( adz, voxel ); 
	GET_VOXEL_3D(voxel, dz,  i,j,k); value2 = CONVERT_VOXEL_TO_VALUE( dz, voxel ); value2 += value*weight;
	voxel = CONVERT_VALUE_TO_VOXEL(adz, value2);
	SET_VOXEL_3D(adz, i,j,k, voxel); 
		
      }
    }
  }
  
}
					    


/*******************************************************************************/
/*  procedure: smooth_the_warp

    desc: this procedure will smooth the current warp stored in warp_d? and 
          return the smoothed warp in smooth_d?

    meth: smoothing is accomplished by averaging the 6 neighbour of each node
          with the value at that node.

	  new_val = FRAC1*old_val + (1-FRAC1)*sum_6_neighbours(val);
    
*/
public void smooth_the_warp(Volume smooth_dx, Volume smooth_dy, Volume smooth_dz, 
			    Volume warp_dx, Volume warp_dy, Volume warp_dz,
			    Volume warp_mag, Real thres) 
{
  
  int 
    i,j,k,loop_start[3], loop_end[3],sizes[3];
  Real 
    mag,
    wx, wy, wz, mx, my, mz, tx,ty,tz,
    voxel;
  progress_struct
    progress;

  get_volume_sizes(warp_dx, sizes);
  
  for_less(i,0,3) {
    if (sizes[i]>3) {
      loop_start[i] = 1;
      loop_end[i] = sizes[i]-1;
    }
    else {
      loop_start[i]=0;
      loop_end[i] = sizes[i];
    }
  }


  initialize_progress_report( &progress, FALSE, 
			     (loop_end[0]-loop_start[0])*(loop_end[1]-loop_start[1]) + 1,
			     "Smoothing deformations" );

  for_less(i,loop_start[0],loop_end[0]) {
    for_less(j,loop_start[1],loop_end[1]) {
      for_less(k,loop_start[2],loop_end[2]){

	GET_VOXEL_3D(voxel, warp_mag, i,j,k); 
	mag = CONVERT_VOXEL_TO_VALUE(warp_mag , voxel ); 
	if (mag >= thres) {
	  
	  convert_3D_voxel_to_world(warp_dx, i, j, k, &tx, &ty, &tz);
	  
	  GET_VOXEL_3D(voxel, warp_dx, i,j,k); 
	  wx = CONVERT_VOXEL_TO_VALUE( warp_dx, voxel ); 
	  GET_VOXEL_3D(voxel, warp_dy, i,j,k); 
	  wy = CONVERT_VOXEL_TO_VALUE( warp_dy, voxel ); 
	  GET_VOXEL_3D(voxel, warp_dz, i,j,k); 
	  wz = CONVERT_VOXEL_TO_VALUE( warp_dz, voxel ); 
	  
	  tx += wx; ty += wy; tz += wz;
	  
	  get_average_warp_of_neighbours(warp_dx, warp_dy, warp_dz, 
					 i,j,k,
					 &mx, &my, &mz);
	  
	  wx += smoothing_weight*(mx - tx);
	  wy += smoothing_weight*(my - ty);
	  wz += smoothing_weight*(mz - tz);
	  
	  voxel = CONVERT_VALUE_TO_VOXEL(smooth_dx, wx);
	  SET_VOXEL_3D(smooth_dx, i,j,k, voxel); 
	  voxel = CONVERT_VALUE_TO_VOXEL(smooth_dy, wy);
	  SET_VOXEL_3D(smooth_dy, i,j,k, voxel); 
	  voxel = CONVERT_VALUE_TO_VOXEL(smooth_dz, wz);
	  SET_VOXEL_3D(smooth_dz, i,j,k, voxel); 
	}
	else {
	  GET_VOXEL_3D(voxel, warp_dx, i,j,k); SET_VOXEL_3D(smooth_dx, i,j,k, voxel); 
	  GET_VOXEL_3D(voxel, warp_dy, i,j,k); SET_VOXEL_3D(smooth_dy, i,j,k, voxel); 
	  GET_VOXEL_3D(voxel, warp_dz, i,j,k); SET_VOXEL_3D(smooth_dz, i,j,k, voxel); 
	}
	  
      }
      update_progress_report( &progress,
			     (loop_end[1]-loop_start[1])*(i-loop_start[0])+(j-loop_start[1])+1  );

    }
    terminate_progress_report( &progress );
  }
}

public void clamp_warp_deriv(Volume dx, Volume dy, Volume dz)
{
  int clamp_needed, clamped_once, i,j,k,sizes[3];
  Real steps[3], diff, voxel, value1, value2, old2 ;
  progress_struct
    progress;

  get_volume_sizes(dx, sizes);
  get_volume_separations(dx, steps);
  
  initialize_progress_report( &progress, FALSE, sizes[0]+sizes[1]+sizes[2] + 1,
			     "Deriv check" );

  clamped_once = FALSE;

  if (sizes[0]>2)
  for_less(j,0,sizes[1]) {
    for_less(k,0,sizes[2]){      

      
      do {
	
	clamp_needed = FALSE;
	
	GET_VOXEL_3D(voxel, dz, 0 ,j,k); value1 = CONVERT_VOXEL_TO_VALUE( dz, voxel ); 
	
	for_less(i,1,sizes[0]) {
	  GET_VOXEL_3D(voxel, dz, i,j,k);   value2 = CONVERT_VOXEL_TO_VALUE( dz, voxel ); 
	  old2 = value2;
	  
	  diff = value1 - (value2+steps[0]);
	  
	  if (diff>0.0) {
	    clamp_needed = TRUE; clamped_once = TRUE;
	    
	    value1 -= diff * DERIV_FRAC;
	    value2 += diff * DERIV_FRAC;
	    
	    voxel = CONVERT_VALUE_TO_VOXEL( dz, value1); SET_VOXEL_3D(dz, i-1,j,k, voxel);
	    voxel = CONVERT_VALUE_TO_VOXEL( dz, value2); SET_VOXEL_3D(dz, i,j,k, voxel);
	  }
	  
	  value1 = old2;
	}

      } while (clamp_needed); 
    }

    update_progress_report( &progress, j+1 );
  }


  if (sizes[1]>2)
  for_less(k,0,sizes[2]){      
    for_less(i,0,sizes[0]) {

      
      do { 
	
      clamp_needed = FALSE;

	GET_VOXEL_3D(voxel, dy, i ,0,k); value1 = CONVERT_VOXEL_TO_VALUE( dy, voxel ); 
	
	for_less(j,1,sizes[1]) {
	  GET_VOXEL_3D(voxel, dy, i,j,k);   value2 = CONVERT_VOXEL_TO_VALUE( dy, voxel ); 
	  old2 = value2;
	  
	  diff = value1 - (value2+steps[1]);
	  
	  if (diff>0.0) {
	    clamp_needed = TRUE; clamped_once = TRUE;

	    value1 -= diff * DERIV_FRAC;
	    value2 += diff * DERIV_FRAC;
	    
	    voxel = CONVERT_VALUE_TO_VOXEL( dy, value1); SET_VOXEL_3D(dy, i,j-1,k, voxel);
	    voxel = CONVERT_VALUE_TO_VOXEL( dy, value2); SET_VOXEL_3D(dy, i,j,k, voxel);
	  }
	  
	  value1 = old2;
	}

      } while (clamp_needed); 
    }
    update_progress_report( &progress, sizes[1]+k+1 );
  }


  if (sizes[2]>2)
  for_less(i,0,sizes[0]) {
    for_less(j,0,sizes[1]) {

      
      do {
	
      clamp_needed = FALSE;

	GET_VOXEL_3D(voxel, dx, i ,j,0); value1 = CONVERT_VOXEL_TO_VALUE( dx, voxel ); 
	
	for_less(k,0,sizes[2]){      
	  GET_VOXEL_3D(voxel, dx, i,j,k);   value2 = CONVERT_VOXEL_TO_VALUE( dx, voxel ); 
	  old2 = value2;
	  
	  diff = value1 - (value2+steps[2]);
	  
	  if (diff>0.0) {
	    clamp_needed = TRUE; clamped_once = TRUE;

	    value1 -= diff * DERIV_FRAC;
	    value2 += diff * DERIV_FRAC;
	    
	    voxel = CONVERT_VALUE_TO_VOXEL( dx, value1); SET_VOXEL_3D(dx, i,j,k-1, voxel);
	    voxel = CONVERT_VALUE_TO_VOXEL( dx, value2); SET_VOXEL_3D(dx, i,j,k, voxel);
	  }
	  
	  value1 = old2;
	}

      } while (clamp_needed);
    }
    update_progress_report( &progress, sizes[2]+sizes[1]+i+1 );

  }


  if (clamped_once)
    print ("Clamped once!\n");

}


/*******************************************************************************/
/*  procedure: smooth_the_warp

    desc: this procedure will smooth the current warp stored in warp_d? and 
          return the smoothed warp in smooth_d?

    meth: smoothing is accomplished by averaging the 6 neighbour of each node
          with the value at that node.

	  new_val = FRAC1*old_val + (1-FRAC1)*sum_6_neighbours(val);
    
*/
public void smooth2_the_warp(Volume smooth_dx, Volume smooth_dy, Volume smooth_dz, 
			     Volume warp_dx, Volume warp_dy, Volume warp_dz) {

  int i,j,k,sizes[3];
  Real voxel, value, value2;
  progress_struct
    progress;

  get_volume_sizes(warp_dx, sizes);
  
  initialize_progress_report( &progress, FALSE, sizes[0]*sizes[1] + 1,
			     "Smoothing deformations" );

  for_less(i,1,sizes[0]-1)
    for_less(j,1,sizes[1]-1) {
      for_less(k,1,sizes[2]-1){
	
	GET_VOXEL_3D(voxel, warp_dx, i,j,k); 
	value = CONVERT_VOXEL_TO_VALUE( warp_dx, voxel ); value *= FRAC1;
	GET_VOXEL_3D(voxel, warp_dx, i+1,j,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dx, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dx, i-1,j,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dx, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dx, i,j+1,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dx, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dx, i,j-1,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dx, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dx, i,j,k+1); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dx, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dx, i,j,k-1); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dx, voxel ); value += value2*FRAC2; 
	voxel = CONVERT_VALUE_TO_VOXEL(smooth_dx, value);
	SET_VOXEL_3D(smooth_dx, i,j,k, voxel); 
	
	GET_VOXEL_3D(voxel, warp_dy, i,j,k); 
	value = CONVERT_VOXEL_TO_VALUE( warp_dy, voxel ); value *= FRAC1;
	GET_VOXEL_3D(voxel, warp_dy, i+1,j,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dy, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dy, i-1,j,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dy, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dy, i,j+1,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dy, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dy, i,j-1,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dy, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dy, i,j,k+1); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dy, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dy, i,j,k-1); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dy, voxel ); value += value2*FRAC2; 
	voxel = CONVERT_VALUE_TO_VOXEL(smooth_dy, value);
	SET_VOXEL_3D(smooth_dy, i,j,k, voxel); 
	
	GET_VOXEL_3D(voxel, warp_dz, i,j,k); 
	value = CONVERT_VOXEL_TO_VALUE( warp_dz, voxel ); value *= FRAC1;
	GET_VOXEL_3D(voxel, warp_dz, i+1,j,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dz, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dz, i-1,j,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dz, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dz, i,j+1,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dz, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dz, i,j-1,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dz, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dz, i,j,k+1); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dz, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dz, i,j,k-1); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dz, voxel ); value += value2*FRAC2; 
	voxel = CONVERT_VALUE_TO_VOXEL(smooth_dz, value);
	SET_VOXEL_3D(smooth_dz, i,j,k, voxel); 
	
      }
      update_progress_report( &progress, sizes[1]*i+j+1 );

    }
    terminate_progress_report( &progress );


  
}

/*   get the value of the maximum gradient magnitude */

public Real get_maximum_magnitude(Volume dxyz)
{
  int i,j,k,sizes[3];
  Real voxel, val, max;

  max = -DBL_MAX;
  get_volume_sizes(dxyz, sizes);
  
  for_less(i,0,sizes[0])
    for_less(j,0,sizes[1])
      for_less(k,0,sizes[2]){
	
	GET_VOXEL_3D(voxel, dxyz, i,j,k); val = CONVERT_VOXEL_TO_VALUE(dxyz, voxel);

	if (val > max) max = val;

      }

  return(max);

}

/* copied from tricubic_interpolant in interpolation.c */

public int tricubic_slice_interpolant(Volume volume, 
				      PointR *coord, double *result)
{
  Real
    v00,v01,v02,v03, 
    v10,v11,v12,v13, 
    v20,v21,v22,v23, 
    v30,v31,v32,v33;
   long 
     ind0, ind1, ind2, max[3];
   double 
     frac[3];
   int 
     sizes[3];
   int 
     flag;
   double temp_result;

   /* Check that the coordinate is inside the volume */

   get_volume_sizes(volume, sizes);
   max[0] = sizes[0];
   max[1] = sizes[1];
   max[2] = sizes[2];
   
   if ((Point_y( *coord ) < 0) || (Point_y( *coord ) >= max[1]-1) ||
       (Point_z( *coord ) < 0) || (Point_z( *coord ) >= max[2]-1)) {

     flag = nearest_neighbour_interpolant(volume, coord, &temp_result) ;
     *result = temp_result;
     return(flag);
   }


   /* Get the whole and fractional part of the coordinate */
   ind0 = (long) Point_x( *coord );
   ind1 = (long) Point_y( *coord );
   ind2 = (long) Point_z( *coord );
   frac[1] = Point_y( *coord ) - ind1;
   frac[2] = Point_z( *coord ) - ind2;
   ind1--;
   ind2--;

   /* Check for edges - do linear interpolation at edges */
   if ((ind1 >= max[1]-3) || (ind1 < 0) ||
       (ind2 >= max[2]-3) || (ind2 < 0)) {
      return trilinear_interpolant(volume, coord, result);
   }
   
   GET_VOXEL_3D(v00, volume, ind0, ind1, ind2);
   GET_VOXEL_3D(v01, volume, ind0, ind1, ind2+1);
   GET_VOXEL_3D(v02, volume, ind0, ind1, ind2+2);
   GET_VOXEL_3D(v03, volume, ind0, ind1, ind2+3);

   GET_VOXEL_3D(v10, volume, ind0, ind1+1, ind2);
   GET_VOXEL_3D(v11, volume, ind0, ind1+1, ind2+1);
   GET_VOXEL_3D(v12, volume, ind0, ind1+1, ind2+2);
   GET_VOXEL_3D(v13, volume, ind0, ind1+1, ind2+3);

   GET_VOXEL_3D(v20, volume, ind0, ind1+2, ind2);
   GET_VOXEL_3D(v21, volume, ind0, ind1+2, ind2+1);
   GET_VOXEL_3D(v22, volume, ind0, ind1+2, ind2+2);
   GET_VOXEL_3D(v23, volume, ind0, ind1+2, ind2+3);

   GET_VOXEL_3D(v30, volume, ind0, ind1+3, ind2);
   GET_VOXEL_3D(v31, volume, ind0, ind1+3, ind2+1);
   GET_VOXEL_3D(v32, volume, ind0, ind1+3, ind2+2);
   GET_VOXEL_3D(v33, volume, ind0, ind1+3, ind2+3);

   CUBIC_BIVAR( v00,v01,v02,v03, v10,v11,v12,v13, v20,v21,v22,v23, v30,v31,v32,v33, frac[1],frac[2], temp_result);


  *result = CONVERT_VOXEL_TO_VALUE(volume,temp_result);

   return(TRUE);

}



public void make_super_sampled_data(Volume orig_dx, Volume *dx, Volume *dy, Volume *dz)

{

  int 
    i,
    sizes[3], sizes2[3];
  Real
    voxel[3],
    start[3],
    steps[3],
    steps1[3],
    steps2[3];
  VectorR 
    directions[3];

    
  get_volume_sizes(orig_dx, sizes);
  get_volume_separations(orig_dx, steps);
  
  for_less(i,0,3) {
    steps1[i] = steps[i]/2;
    if(sizes[i]==1) {
      steps2[i] = 1000;
    }
    else {
      steps2[i] = steps[i]/2;
    }
  }
  
  *dx = copy_volume_definition(orig_dx, NC_UNSPECIFIED, FALSE, 0.0,0.0);

  set_up_lattice(*dx, steps2, start, sizes2, steps, directions);

  for_less(i,0,3)		/* use the sign of the step returned to set the true step size */
    if (steps[i]<0) steps[i] = -1.0 * ABS(steps1[i]); else steps[i] = ABS(steps1[i]);
  
  set_volume_sizes(*dx, sizes2);
  set_volume_separations(*dx,steps);
  voxel[0] = 0.0;  voxel[1] = 0.0;  voxel[2] = 0.0;
  set_volume_translation(*dx, voxel, start);
  
  *dy = copy_volume_definition(*dx, NC_UNSPECIFIED, FALSE, 0.0,0.0);
  *dz = copy_volume_definition(*dx, NC_UNSPECIFIED, FALSE, 0.0,0.0);
  
  alloc_volume_data(*dx);
  alloc_volume_data(*dy);
  alloc_volume_data(*dz);
  

}

public void interpolate_super_sampled_data(Volume orig_dx, Volume orig_dy, Volume orig_dz,
					   Volume dx, Volume dy, Volume dz,
					   int dim)

{

  int 
    i,j,k,
    sizes2[3];
  Real
    wx,wy,wz,
    nx,ny,nz,
    voxel_val,
    val_x, val_y, val_z;
  PointR voxel_point; 
  progress_struct
    progress;

    
  get_volume_sizes(dx, sizes2);
  
  if (dim==3) {
    initialize_progress_report( &progress, FALSE, sizes2[0],
			       "Super-sampling defs:" );
    for_less(i,0,sizes2[0]) {
      for_less(j,0,sizes2[1]) {
	for_less(k,0,sizes2[2]){
	  convert_3D_voxel_to_world(dx, (Real)i, (Real)j, (Real)k,
				    &wx,&wy,&wz);
	  convert_3D_world_to_voxel(orig_dx,wx,wy,wz,&nx,&ny,&nz);
	  
	  fill_Point( voxel_point, nx, ny, nz  ); /* build the voxel POINT */
	  
	  if ( !tricubic_interpolant( orig_dx, &voxel_point, &val_x ))
	    val_x = 0.0;
	  if ( !tricubic_interpolant( orig_dy, &voxel_point, &val_y ))
	    val_y = 0.0;
	  if ( !tricubic_interpolant( orig_dz, &voxel_point, &val_z ))
	    val_z = 0.0;
	  voxel_val = CONVERT_VALUE_TO_VOXEL(dx, val_x); SET_VOXEL_3D(dx, i,j,k, voxel_val);
	  voxel_val = CONVERT_VALUE_TO_VOXEL(dy, val_y); SET_VOXEL_3D(dy, i,j,k, voxel_val);
	  voxel_val = CONVERT_VALUE_TO_VOXEL(dz, val_z); SET_VOXEL_3D(dz, i,j,k, voxel_val);
	
	  
	}
      }
      update_progress_report( &progress, i+1 );
    }
    terminate_progress_report( &progress );
  }
  else {
    initialize_progress_report( &progress, FALSE, sizes2[1],
			       "Super-sampling defs:" );

    for_less(j,0,sizes2[1]) {
      for_less(k,0,sizes2[2]){
	convert_3D_voxel_to_world(dx, 0.0, (Real)j, (Real)k,
				  &wx,&wy,&wz);
	convert_3D_world_to_voxel(orig_dx,wx,wy,wz,&nx,&ny,&nz);
	
	fill_Point( voxel_point, nx, ny, nz  ); /* build the voxel POINT */
	
	if ( !tricubic_slice_interpolant( orig_dx, &voxel_point, &val_x ))
	  val_x = 0.0;
	if ( !tricubic_slice_interpolant( orig_dy, &voxel_point, &val_y ))
	  val_y = 0.0;
	if ( !tricubic_slice_interpolant( orig_dz, &voxel_point, &val_z ))
	  val_z = 0.0;
	voxel_val = CONVERT_VALUE_TO_VOXEL(dx, val_x); SET_VOXEL_3D(dx, 0,j,k, voxel_val);
	voxel_val = CONVERT_VALUE_TO_VOXEL(dy, val_y); SET_VOXEL_3D(dy, 0,j,k, voxel_val);
	voxel_val = CONVERT_VALUE_TO_VOXEL(dz, val_z); SET_VOXEL_3D(dz, 0,j,k, voxel_val);
	
	
      }
      update_progress_report( &progress, j+1 );
    }
    
    terminate_progress_report( &progress );
  }


}


public Real get_value_of_point_in_volume(Real xw, Real yw, Real zw, 
					  Volume data)
     
{

  Real
    mag,
    xvox, yvox, zvox;
  
  PointR 
    voxel;

  convert_3D_world_to_voxel(data, xw, yw, zw,
			    &zvox, &yvox, &xvox);
  fill_Point( voxel, zvox, yvox, xvox );

  if (!trilinear_interpolant(data, &voxel, &mag)) 
    return(-DBL_MAX); 
  else 
    return(mag);
}

/*******************************************************************************/
/* debugging procedure called to save the deformation at each iterative step  */

public void save_data(char *basename, int i, int j,
		      General_transform *transform)
{

  Status 
    status;
  STRING 
    comments,name;
  FILE *file;

  (void)sprintf(comments,"step %d of %d of the non-linear estimation",i,j);

  status = open_file_with_default_suffix(basename,
					 get_default_transform_file_suffix(),
					 WRITE_FILE, ASCII_FORMAT, &file );
  
  if( status == OK )
    status = output_transform(file,
			      basename,
			      &i,
			      comments,
			      transform);
  
  if( status == OK )
    status = close_file( file );
  
  if (status!=OK)
    print ("Error saving %s\n",basename);
}

public void    build_source_lattice(Real x, Real y, Real z,
				    float PX[], float PY[], float PZ[],
				    Real width_x, Real width_y, Real width_z, 
				    int nx, int ny, int nz,
				    int ndim, int *length)
{
  int 
    c, 
    i,j,k;
  float 
    radius_squared,
    tx,ty,tz;

  *length = 0; 
  c = 1;
  radius_squared = 0.55 * 0.55;	/* a bit bigger than .5^2 */
  
  if (ndim==2) {
    for_less(i,0,nx)
      for_less(j,0,ny) {
	tx = -0.5 + (float)(i)/(float)(nx-1);
	ty = -0.5 + (float)(j)/(float)(ny-1);
	if ((tx*tx + ty*ty) <= radius_squared) {
	  PX[c] = (float)(x + width_x * tx);
	  PY[c] = (float)(y + width_y * ty);
	  PZ[c] = (float)z;
	  c++;
	  (*length)++;
	}
      }
  }
  else {
    for_less(i,0,nx)
      for_less(j,0,ny)
	for_less(k,0,nz) {
	  tx = -0.5 + (float)(i)/(float)(nx-1);
	  ty = -0.5 + (float)(j)/(float)(ny-1);
	  tz = -0.5 + (float)(k)/(float)(nz-1);
	  if ((tx*tx + ty*ty + tz*tz) <= radius_squared) {
	    PX[c] = (float)(x + width_x * tx);
	    PY[c] = (float)(y + width_y * ty);
	    PZ[c] = (float)(z + width_z * tz);
	    c++;
	    (*length)++;
	  }
	}
  }
}

public void go_get_samples_in_source(Volume data,
				     float x[], float y[], float z[],
				     float samples[],
				     int len,
				     int inter_type) 
{
  int 
    flag, c;
  Real 
    val;
  PointR
    point;


  for_inclusive(c,1,len) {
    convert_3D_world_to_voxel(data, (Real)x[c], (Real)y[c], (Real)z[c], 
			      &Point_x(point), &Point_y(point), &Point_z(point));

    if (inter_type==0)
	flag = nearest_neighbour_interpolant(data,  &point , &val);
    else {
      if (inter_type==1)
	flag = trilinear_interpolant(data,  &point , &val);
      else
	flag = tricubic_interpolant(data,  &point , &val);
    }

    if (flag)
      samples[c] = (float)val;
    else
      samples[c] = 0.0;
  }

}


public float go_get_samples_with_offset(Volume data,
				       float *x, float *y, float *z,
				       Real dx, Real dy, Real dz,
				       int len, float sqrt_s1, float *a1)
{
  float
    sample, r,
    s1,s3;			/* to store the sums for f1,f2,f3 */
  int 
    sizes[3],
    ind0, ind1, ind2, 
    c;  
  int xs,ys,zs;
  float
    f_trans, f_scale;

  unsigned char  ***byte_ptr;
  unsigned short ***ushort_ptr;
           short ***sshort_ptr;
  
  get_volume_sizes(data, sizes);  
  xs = sizes[0];  
  ys = sizes[1];  
  zs = sizes[2];

  f_trans = data->real_value_translation;
  f_scale = data->real_value_scale;

  dx += 0.5;
  dy += 0.5;
  dz += 0.5;

  s1 = 0.0;
  s3 = 0.0;

  ++a1; 

  switch( data->data_type ) {
  case UNSIGNED_BYTE:  

    byte_ptr = data->data;

    ++x; ++y; ++z; 

    for_inclusive(c,1,len) {
      
      ind0 = (int) ( *x++ + dz );
      ind1 = (int) ( *y++ + dy );
      ind2 = (int) ( *z++ + dx );
      
      if (ind0>=0 && ind0<xs &&
	  ind1>=0 && ind1<ys &&
	  ind2>=0 && ind2<zs) {
	sample = byte_ptr[ind0][ind1][ind2] * f_scale + f_trans;
      }
      else
	sample = 0.0;
      
      s1 += *a1++ * sample;
      s3 += sample * sample;
      

    }
    break;
  case SIGNED_SHORT:  

    sshort_ptr = data->data;

    ++x; ++y; ++z; 

    for_inclusive(c,1,len) {
      
      ind0 = (int) ( *x++ + dz );
      ind1 = (int) ( *y++ + dy );
      ind2 = (int) ( *z++ + dx );
      
      if (ind0>=0 && ind0<xs &&
	  ind1>=0 && ind1<ys &&
	  ind2>=0 && ind2<zs) {
	sample = sshort_ptr[ind0][ind1][ind2] * f_scale + f_trans;
      }
      else
	sample = 0.0;
      
      s1 += *a1++ * sample;
      s3 += sample * sample;
      

    }
    break;
  case UNSIGNED_SHORT:  

    ushort_ptr = data->data;

    ++x; ++y; ++z; 

    for_inclusive(c,1,len) {
      
      ind0 = (int) ( *x++ + dz );
      ind1 = (int) ( *y++ + dy );
      ind2 = (int) ( *z++ + dx );
      
      if (ind0>=0 && ind0<xs &&
	  ind1>=0 && ind1<ys &&
	  ind2>=0 && ind2<zs) {
	sample = ushort_ptr[ind0][ind1][ind2] * f_scale + f_trans;
      }
      else
	sample = 0.0;
      
      s1 += *a1++ * sample;
      s3 += sample * sample;
      

    }
    break;
  default:
    print_error_and_line_num("Data type not supported in go_get_samples_with_offset",__FILE__, __LINE__);
  }



  if ( sqrt_s1 < 0.01 && s3 < 0.0001) {
    r = 1.0;
  }
  else {
    if ( sqrt_s1 < 0.01 || s3 < 0.0001) {
      r = 0.0;
    }
    else {
      r = s1 / (sqrt_s1*sqrt((double)s3));
    }
  }

  return(r);
}


