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
@MODIFIED   : Revision 1.3  1994-06-06 18:46:53  louis
@MODIFIED   : working version: clamp and blur of deformation lattice now ensures
@MODIFIED   : a smooth recovered deformation.  Unfortunately, the r = cost-similarity
@MODIFIED   : function used in the optimization is too heavy on the cost_fn.  This has
@MODIFIED   : to get fixed...
@MODIFIED   :
@MODIFIED   :
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
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Optimize/deform_support.c,v 1.3 1994-06-06 18:46:53 louis Exp $";
#endif

#include "volume_io.h"
#include "limits.h"
#include "line_data.h"

#define DERIV_FRAC      0.6
#define FRAC1           0.5
#define FRAC2           0.0833333

typedef Real (*Difference_Function) 
     (Line_data *m, Line_data *d, Real offset);



/*
% return the value of the 1st or second derivative of the data, based on cubic
% interpolation through the four equally spaced data function value
% points.
% xpos is the fractional distance from the x position of f2 for which the
% value is to be returned.
function v = cubic_p(f1,f2,f3,f4,xpos)

%  for f1 = f(-1),
%      f2 = f(0),
%      f3 = f(1),
%      f4 = f(2), 
%
%      a1 =                f2;
%      a2 =  -1.0*f1/3.0 - f2/2.0 + f3     - f4/6.0;
%      a3 =       f1/2.0 - f2     + f3/2.0 ;
%      a4 =  -1.0*f1/6.0 + f2/2.0 - f3/2.0 + f4/6.0;
%
%  to estimate function f() at xpos, using cubic interpolation:
%
%      f   = a1 + a2*xpos + a3*xpos^2 + a4*xpos^3
%      fp  = a2 + 2*a3*xpos + 3*a4*xpos*xpos;      first derivitive
%      fpp = 2*a3 + 6*a4*xpos;	                   second derivitive 
*/

public Real CUBIC_PP(Real f1, Real f2, Real f3, Real f4, Real xpos) {
  Real v,a3,a4;
  
  a3 =      f1/2.0 - f2     + f3/2.0 ;
  a4 = -1.0*f1/6.0 + f2/2.0 - f3/2.0 + f4/6.0;
  
  v =  2*a3 + 6*a4*xpos;
  
  return(v);
}

public Real CUBIC_P(Real f1, Real f2, Real f3, Real f4, Real xpos) {

  Real v,a2,a3,a4;
  
  a2 = -1.0*f1/3.0 - f2/2.0 + f3     - f4/6.0;
  a3 =      f1/2.0 - f2     + f3/2.0 ;
  a4 = -1.0*f1/6.0 + f2/2.0 - f3/2.0 + f4/6.0;
  
  v =  a2 + 2*a3*xpos + 3*a4*xpos*xpos;
  
  return(v);
}



/*
  real space 2nd derivative value interpolation
*/
public Real inter_2p(Line_data *line, Real real_p) {

  Real pos, u, v, v0,v1,v2,v3;
  int p,m;
				/* establish voxel position */
  
  pos = (real_p - line->start) / line->step;
  m = line->count;
  
  p = FLOOR(pos);
  u = pos-p;
  
  v=0;
  
  if (pos>=1 &&  pos<(m-2)) {
    v0 = line->data[p-1];
    v1 = line->data[p];
    v2 = line->data[p+1];
    v3 = line->data[p+2];
    
    v = CUBIC_PP(v0,v1,v2,v3,u);
    
  } else if (pos>=0 &&  pos<1) {
    v0 = line->data[p];
    v1 = line->data[p+1];
    v2 = line->data[p+2];
    v3 = line->data[p+3];
    u=u+1;

    v = CUBIC_PP(v0,v1,v2,v3,u);

  } else if (pos>=(m-2) &&  pos<(m-1)) {
    v0 = line->data[p-2];
    v1 = line->data[p-1];
    v2 = line->data[p];
    v3 = line->data[p+1];
    u=u-1;
    
    v = CUBIC_PP(v0,v1,v2,v3,u);
  }

  return(v);
}




/* 
  real space 1st derivative value interpolation
*/ 
public Real inter_1p(Line_data *line, Real real_p) {

  Real pos, u, v, v0,v1,v2,v3;
  int m,p;
				/* establish voxel position */
  
  pos = (real_p - line->start) / line->step;
  m = line->count;
  
  p = FLOOR(pos);
  u = pos-p;
  
  v=0;
  
  if (pos>=1 &&  pos<(m-2)) {
    v0 = line->data[p-1];
    v1 = line->data[p];
    v2 = line->data[p+1];
    v3 = line->data[p+2];
    
    v = CUBIC_P(v0,v1,v2,v3,u);
    
  } else if (pos>=0 &&  pos<1) {
    v0 = line->data[p];
    v1 = line->data[p+1];
    v2 = line->data[p+2];
    v3 = line->data[p+3];
    u=u+1;

    v = CUBIC_P(v0,v1,v2,v3,u);

  } else if (pos>=(m-2) &&  pos<(m-1)) {
    v0 = line->data[p-2];
    v1 = line->data[p-1];
    v2 = line->data[p];
    v3 = line->data[p+1];
    u=u-1;
    
    v = CUBIC_P(v0,v1,v2,v3,u);
  }

  return(v);
}







/*
  voxel space value interpolation
*/
public Real inter(Line_data *line, Real pos) {

  Real v,u,v0,v1,v2,v3,f;
  int m,p;

  m = line->count;
  v = 0;

  if (pos>=1 &&  pos<(m-2)) {
    p = FLOOR(pos);
    u = pos-p;
    
    v0 = line->data[p-1];
    v1 = line->data[p];
    v2 = line->data[p+1];
    v3 = line->data[p+2];
    
    v = ( (v1) + (u) * ( 0.5 * ((v2)-(v0)) + (u) * ( (v0) - 2.5 * (v1) + 2.0 * (v2) - 0.5 * (v3) + (u) * ( -0.5 * (v0) + 1.5 * (v1) - 1.5 * (v2) + 0.5 * (v3)  ) ) ) );
  }
  else {
    if (pos<-0.5 || pos>=(m-0.5))
      v = 0;
    else {
      if (pos<0 || pos>(m-1)) {
	if (pos<0)
	  v = 0.5 * line->data[0];
	else
	  v = 0.5 * line->data[m-1];
      }
      else {
	if (pos<1) {
	  f = pos;
	  v = (1.0-f)*line->data[0] + f*line->data[1];
	} else
	  if (pos>=m-2) {
	    f = pos-(m-2);
	    v = (1.0-f)*line->data[m-2] + f*line->data[m-1];
	  }
      }
    }
  }
  
  return(v);
}

/*
  real space value interpolation
*/
public Real inter_p(Line_data *data, Real real_p) {

  Real pos,v;

  pos = (real_p - data->start) / data->step;

  v = inter(data, pos);

  return(v);
}



/*  return the difference between the overlapping area
    of the two functions of model+offset wrt data.
*/
public Real difference(Line_data *model,
		      Line_data *data,
		      Real offset) {

  Real t, d, val, minx, maxx;
  int i;
  
  val = 0;
  minx = data->start;
  maxx = data->start + data->count*data->step;

  for_less(i,0,model->count) {

    t = model->start + i*model->step + offset; /* new position */

    if (t>=minx && t <=maxx) {
      d = model->data[i] - inter_p(data, t);
      val = val + d*d;
    }
  }
  return(val);
}

/* procedure: brute
     use brute force to search for the offset that minimizes the difference
     between the data in m and d, where m is a spatial subset of d.
   method:
     search with offsets from -limit to limit, stepping by the value in step
     and return the best offset position.
*/
public Real brute(Difference_Function differ,
		  Line_data *m,
		  Line_data *d,
		  Real step,
		  Real limit)
{

  Real 
    v_min, v, pos,
    dist;

  v_min = FLT_MAX;

  dist = 0;

  pos = -limit;

  while (pos<=limit) {


    v = (*(differ))(m,d,pos);

/*
print ("%6.2f %10.8f\n",pos,v*1000);
*/
    if (v<v_min) {
      v_min = v;
      dist = pos;
    }
    pos += step;
  }

  pos = dist;

  if (dist<limit) {
    v = (*(differ))(m,d,pos+step/2);

    if (v<v_min) {
      v_min = v;
      dist = pos+step/2;
    }
  }

  if (dist> -limit) {
    v = (*(differ))(m,d,pos-step/2);
    if (v<v_min) {
      v_min = v;
      dist = pos-step/2;
    }
  }

  return(dist);
}

#define R  0.61803399
#define C  1.0-R

public Real gold(Difference_Function differ,
		 Line_data *m,
		 Line_data *d,
		 Real limit)
{

  Real /*f0,*/f1,f2,/*f3,*/d1,d2,x0,x3,x1,x2,ax,bx,cx,dist;
  
  int i;

  i = 0;
  dist = 0;

  ax = -0.6*limit;
  cx = 0.6*limit;


  
  x0 = ax;
  bx = 0;
  x3 = cx;

  d1=cx-bx;
  d2=bx=ax;

  if ( ABS(d1) > ABS(d2)) {
    x1=bx;
    x2=bx+C*(cx-bx);
  }
  else {
    x2=bx;
    x1=bx-C*(bx-ax);
  }

  f1= (*(differ))(m,d,x1);
  f2= (*(differ))(m,d,x2);

  while (abs(x3-x0) > limit/16 && i<20) {
    i++;
    if (f2 < f1) {
      x0=x1; x1=x2; x2=R*x1+C*x3;
      /* f0=f1;  */
      f1=f2; f2= (*(differ))(m,d,x2);
    }
    else {
      x3=x2; x2=x1; x1=R*x2+C*x0;
      /* f3=f2; */
      f2=f1; f1= (*(differ))(m,d,x1);
    }
  }
  
  if (f1 < f2)
    dist=x1;
  else
    dist=x2;
  
  return(dist);
  
  
}


/*
  find the  offset to best match model (m) with data (d)

  offset returned is constrained to: 
       | offset | <= limit

  if op_type = 1, do brute force, otherwise use secant
*/
public Real find_offset_to_match(Line_data *m,
				 Line_data *d,
				 Real      limit,
				 int       op_type) 
{

  Real dist,step;

  if (op_type==1) {
    step = limit/10;
    dist = brute(difference,m,d,step,limit);
  }
  else
    dist = gold(difference,m,d,limit);

  return(dist);
}



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

  if (sizes[0]>1) {
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

  if (sizes[1]>1) {
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

  if (sizes[2]>1) {
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
public void smooth_the_warp(Volume smooth_dx, Volume smooth_dy, Volume smooth_dz, 
			    Volume warp_dx, Volume warp_dy, Volume warp_dz) 
{
  
  int 
    i,j,k,loop_start[3], loop_end[3],sizes[3];
  Real 
    wx, wy, wz, mx, my, mz, tx,ty,tz,
    voxel, value, value2;
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

	wx += (1.0 - FRAC1)*(mx - tx);
	wy += (1.0 - FRAC1)*(my - ty);
	wz += (1.0 - FRAC1)*(mz - tz);

	voxel = CONVERT_VALUE_TO_VOXEL(smooth_dx, wx);
	SET_VOXEL_3D(smooth_dx, i,j,k, voxel); 
	voxel = CONVERT_VALUE_TO_VOXEL(smooth_dy, wy);
	SET_VOXEL_3D(smooth_dy, i,j,k, voxel); 
	voxel = CONVERT_VALUE_TO_VOXEL(smooth_dz, wz);
	SET_VOXEL_3D(smooth_dz, i,j,k, voxel); 

      }
      update_progress_report( &progress, sizes[1]*i+j+1 );

    }
    terminate_progress_report( &progress );
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

public Real get_maximum_magnitude(Volume dx, Volume dy, Volume dz)
{
  int i,j,k,sizes[3];
  Real vx, vy, vz, x, y, z, val, max;

  max = -DBL_MAX;
  get_volume_sizes(dx, sizes);

  for_less(i,0,sizes[0])
    for_less(j,0,sizes[1])
      for_less(k,0,sizes[2]){

	GET_VOXEL_3D(vx, dx, i,j,k); x = CONVERT_VOXEL_TO_VALUE(dx, vx);
	GET_VOXEL_3D(vy, dy, i,j,k); y = CONVERT_VOXEL_TO_VALUE(dy, vy);
	GET_VOXEL_3D(vz, dz, i,j,k); z = CONVERT_VOXEL_TO_VALUE(dz, vz);
	val = x*x + y*y + z*z;

	if (val > max) max = val;

      }

  if (max>0.0)
    max = sqrt(max);

  return(max);

}
