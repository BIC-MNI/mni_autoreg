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
@MODIFIED   : Revision 1.1  1994-04-06 11:47:27  louis
@MODIFIED   : Initial revision
@MODIFIED   :

---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Optimize/deform_support.c,v 1.1 1994-04-06 11:47:27 louis Exp $";
#endif

#include "volume_io.h"
#include "limits.h"
#include "line_data.h"

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

private Real CUBIC_PP(Real f1, Real f2, Real f3, Real f4, Real xpos) {
  Real v,a3,a4;
  
  a3 =      f1/2.0 - f2     + f3/2.0 ;
  a4 = -1.0*f1/6.0 + f2/2.0 - f3/2.0 + f4/6.0;
  
  v =  2*a3 + 6*a4*xpos;
  
  return(v);
}

private Real CUBIC_P(Real f1, Real f2, Real f3, Real f4, Real xpos) {

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
				 int       op_type) {


  Real dist,step;

  if (op_type==1) {
    step = limit/10;
    dist = brute(difference,m,d,step,limit);
  }
  else
    dist = gold(difference,m,d,limit);

  return(dist);
}



