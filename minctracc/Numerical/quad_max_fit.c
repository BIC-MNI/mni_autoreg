/* ----------------------------- MNI Header -----------------------------------
@NAME       : quad_max_fit.c
@DESCRIPTION: routines to find maximum of a 3x3 or 3x3x3 region using
              2D or 3D quadratic estimation
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

@CREATED    : Mon May 22 14:14:50 MET DST 1995
@MODIFIED   : $Log: quad_max_fit.c,v $
@MODIFIED   : Revision 1.1  1995-06-12 14:30:23  louis
@MODIFIED   : Initial revision
@MODIFIED   :

---------------------------------------------------------------------------- */

#include <volume_io.h>		

    /* local structures */

typedef struct {
  Real 
    u,v,w,
    uu,vv,ww,
    uv,uw,vw;
} deriv_3D_struct;


typedef struct {
  Real 
    u,v,
    uu,vv,
    uv;
} deriv_2D_struct;

#define MINIMUM_DET_ALLOWED 0.000001

    /* local prototypes */

private BOOLEAN negative_2D_definite(deriv_2D_struct *c);
private BOOLEAN negative_3D_definite(deriv_3D_struct *c);
private BOOLEAN positive_3D_definite(deriv_3D_struct *c);
public   void    estimate_3D_derivatives(Real r[3][3][3], 
					deriv_3D_struct *c);	     
private void    estimate_2D_derivatives(Real r[3][3], 
					deriv_2D_struct *c);	     

/* this procedure will return TRUE with the dx,dy,dz (offsets) that 
   correspond to the maximum value of the quadratic function fit through 
   data points represented in r[u][v][w], unless the maxtrix represented in 
   r[u][v][w] is negative definite
   in which case, the function will return FALSE.                         */

public BOOLEAN return_3D_disp_from_quad_fit(Real r[3][3][3], 
					    Real *dispu, 
					    Real *dispv, 
					    Real *dispw)	
{
  deriv_3D_struct 
    d;			/* the 1st and second order derivatives */
  Real
    inv[3][3],A[3][3],
    a[3][3],
    detA;
  int 
    i,j,k,return_flag;

  *dispu = *dispv = *dispw = 0.0;

  estimate_3D_derivatives(r,&d);

  /*    /          \    the values that form the matrix A come from the  */
  /*    | uu uv uw |    second order derivatives in 'd'.  If this matrix */
  /* A= | uv vv vw |    is not negative definite, then we can not find   */  
  /*    | uw vw ww |    a maximum in the region defined by r[][][].      */
  /*    \          /                                                     */


  if ( !negative_3D_definite(&d) )  {
print ("+");
    return(FALSE);		/* no disp can be calculated */
  }
  else {			/* otherwise, I have the derivatives, now get the 
				   displacements that correspond to a maximum
				   in r[][][]                                     */

    detA = d.uu * (d.vv*d.ww - d.vw*d.vw) -            
           d.uv * (d.uv*d.ww - d.vw*d.uw) + 	       
	   d.uw * (d.uv*d.vw - d.vv*d.uw) ;	       
						       
    if ( ABS( detA ) <= MINIMUM_DET_ALLOWED ) {
print ("0");
      return(FALSE);
    }
    else {
				/* a = inv(A) */

      a[0][0] = (d.vv*d.ww - d.vw*d.vw) / detA;
      a[1][0] = (d.uw*d.vw - d.uv*d.ww) / detA;
      a[2][0] = (d.uv*d.vw - d.uw*d.vv) / detA;
      a[0][1] = (d.vw*d.uw - d.uv*d.ww) / detA;
      a[1][1] = (d.uu*d.ww - d.uw*d.uw) / detA;
      a[2][1] = (d.uw*d.uv - d.uu*d.vw) / detA;
      a[0][2] = (d.uv*d.vw - d.vv*d.uw) / detA;
      a[1][2] = (d.uw*d.uv - d.uu*d.vw) / detA;
      a[2][2] = (d.uu*d.vv - d.uv*d.uv) / detA;

				/* Ax = -b  -->  x = inv(A) (-b) ,
				   where b is the vector of first derivatives*/

      *dispu = -a[0][0]*d.u - a[0][1]*d.v - a[0][2]*d.w;
      *dispv = -a[1][0]*d.u - a[1][1]*d.v - a[1][2]*d.w;
      *dispw = -a[2][0]*d.u - a[2][1]*d.v - a[2][2]*d.w;

/*
   debug test to see if a * A = identity
*/
      A[0][0] = d.uu ; A[0][1] = d.uv ; A[0][2] = d.uw ;
      A[1][0] = d.uv ; A[1][1] = d.vv ; A[1][2] = d.vw ;
      A[2][0] = d.uw ; A[2][1] = d.vw ; A[2][2] = d.ww ;

      for_less(i,0,3)
	for_less(j,0,3) {
	  inv[i][j] = 0.0;
	  for_less(k,0,3)
	    inv[i][j] += a[i][k] * A[k][j];
	}

      
      for_less(i,0,3)
	for_less(j,0,3) {
	  if (((i==j) && ABS( 1.0 - inv[i][j] ) > 0.000001) ||
	      ((i!=j) && ABS( inv[i][j] ) > 0.000001)) {
	    print ("error in inverse for return_3D_disp_from_quad_fit\n");
	    print ("%8.5f %8.5f %8.5f\n",  inv[0][0],inv[0][1],inv[0][2]);
	    print ("%8.5f %8.5f %8.5f\n",  inv[1][0],inv[1][1],inv[1][2]);
	    print ("%8.5f %8.5f %8.5f\n\n",inv[2][0],inv[2][1],inv[2][2]);
	    i=3; j=3;
	  }
	}



      

      if ( ABS( *dispu ) > 2.0 || ABS( *dispv ) > 2.0 || ABS( *dispw ) > 2.0 ) {	

      print ("%8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f\n",
	     r[0][0][0],r[1][0][0],r[2][0][0], 
	     r[0][0][1],r[1][0][1],r[2][0][1], r[0][0][2],r[1][0][2],r[2][0][2]);
      print ("%8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f\n",
	     r[0][1][0],r[1][1][0],r[2][1][0], 
	     r[0][1][1],r[1][1][1],r[2][1][1], r[0][1][2],r[1][1][2],r[2][1][2]);
      print ("%8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f\n",
	     r[0][2][0],r[1][2][0],r[2][2][0], 
	     r[0][2][1],r[1][2][1],r[2][2][1], r[0][2][2],r[1][2][2],r[2][2][2]);


      print ("detA = %f\n", detA);
      print ("/%8.5f %8.5f %8.5f\\  /%8.5f %8.5f %8.5f\\  /%8.5f\\  ->  du = %8.5f\n",  
	     d.uu, d.uv, d.uw, a[0][0],a[0][1],a[0][2], d.u, *dispu);
      print ("|%8.5f %8.5f %8.5f|  |%8.5f %8.5f %8.5f|  |%8.5f|  ->  dv = %8.5f\n",  
	     d.uv, d.vv, d.vw, a[1][0],a[1][1],a[1][2], d.v, *dispv);
      print ("\\%8.5f %8.5f %8.5f/  \\%8.5f %8.5f %8.5f/  \\%8.5f/  ->  dw = %8.5f\n\n",
	     d.uw, d.vw, d.ww, a[2][0],a[2][1],a[2][2], d.w, *dispw);
      
      
	*dispu = *dispv = *dispw = 0.0;
	return(FALSE);
      }
      else {
print ("!");
	return(TRUE);
      }
    }
  }

}

public BOOLEAN return_2D_disp_from_quad_fit(Real r[3][3], /* the values used in the quad fit */
					    Real *dispu, /* the displacements returned */
					    Real *dispv)	
{
  deriv_2D_struct 
    d;			/* the 1st and second order derivatives */
  Real
    a[2][2],
    detA;

  *dispu = *dispv = 0.0;
  estimate_2D_derivatives(r,&d);

  /*    /       \    the values forming the matrix A come from the    */
  /*    | uu uv |    second order derivatives in 'd'.  If this matrix */
  /* A= | uv vv |    is not negative definite, then we can not find   */  
  /*    \       /    a maximum in the region defined by r[][].        */

  if ( !negative_2D_definite(&d) )  {

/*    print ("  not neg\n\n"); */
    return(FALSE);		/* no disp can be calculated */

  }
  else {			/* otherwise, I have the derivatives, now get 
				   the displacements that correspond to a 
				   maximum in r[][]                        */

    detA = d.uu * d.vv - d.uv*d.uv;
						       
    if (ABS( detA) < MINIMUM_DET_ALLOWED) {

/*      print (" (det=%f)\n\n", detA); */
      return(FALSE);

    }
    else {
				/* a = inv(A) */
/*      print (" (det=%f)\n", detA); */

      a[0][0] =  d.vv / detA;
      a[0][1] = -d.uv / detA;
      a[1][0] = -d.uv / detA;
      a[1][1] =  d.uu / detA;

				/* Ax = -b  -->  x = inv(A) (-b),
				   where b is the vector of first derivatives*/

      *dispu = -a[0][0]*d.u - a[0][1]*d.v;
      *dispv = -a[1][0]*d.u - a[1][1]*d.v;

/*
 print ("       %8.5f %8.5f %8.5f \n",	r[0][0], r[1][0], r[2][0]);
 print ("r =    %8.5f %8.5f %8.5f  du = %12.7f\n",r[0][1], r[1][1], r[2][1], *dispu);
 print ("       %8.5f %8.5f %8.5f  dv = %12.7f\n\n",r[0][2], r[1][2], r[2][2], *dispv);
*/

/*
 print ("  / %8.5f %8.5f \\ -1 / %8.5f %8.5f  \\   du = %12.7f\n",   
	d.uu, d.uv, a[0][0], a[0][1], *dispu);
 print ("A=\\ %8.5f %8.5f / A =\\ %8.5f %8.5f /   dv = %12.7f\n\n", 
	d.uv, d.vv, a[0][1], a[1][1], *dispv);
*/

      if ( ABS( *dispu ) > 3.0 || ABS( *dispv ) > 3.0 ) {
	*dispu = *dispv = 0.0;
	return(FALSE);
      }
      else
	return(TRUE);


    }
  }

}



/********************************************************************
   estimate first and second order derivatives by fitting a quadratic
   to given small neighborhood of values stored in r[u][v][w]. 
 ********************************************************************/

public void estimate_3D_derivatives(Real r[3][3][3], 
				   deriv_3D_struct *d) 

{
  int i;
  Real	*p11, *p12, *p13;
  Real	*p21, *p22, *p23;
  Real	*p31, *p32, *p33;

  Real	slice_u1, slice_u2, slice_u3;
  Real	slice_v1, slice_v2, slice_v3;
  Real	slice_w1, slice_w2, slice_w3;
  
  Real	edge_u1_v1, edge_u1_v2, edge_u1_v3;
  Real	edge_u2_v1, edge_u2_v2, edge_u2_v3;
  Real	edge_u3_v1, edge_u3_v2, edge_u3_v3;
  Real	edge_u1_w1, edge_u1_w2, edge_u1_w3;
  Real	edge_u2_w1, edge_u2_w2, edge_u2_w3;
  Real	edge_u3_w1, edge_u3_w2, edge_u3_w3;
  Real	edge_v1_w1, edge_v1_w2, edge_v1_w3;
  Real	edge_v2_w1, edge_v2_w2, edge_v2_w3;
  Real	edge_v3_w1, edge_v3_w2, edge_v3_w3;
  
  /* --- 3x3x3 [u][v][w] --- */
  
  p11 = r[0][0]; 
  p12 = r[0][1]; 
  p13 = r[0][2]; 
  p21 = r[1][0]; 
  p22 = r[1][1]; 
  p23 = r[1][2]; 
  p31 = r[2][0]; 
  p32 = r[2][1]; 
  p33 = r[2][2]; 
  
				/* lines varying along w */
  edge_u1_v1 = ( *p11     + *(p11+1) + *(p11+2));
  edge_u1_v2 = ( *p12     + *(p12+1) + *(p12+2));
  edge_u1_v3 = ( *p13     + *(p13+1) + *(p13+2));
  edge_u2_v1 = ( *p21     + *(p21+1) + *(p21+2));
  edge_u2_v2 = ( *p22     + *(p22+1) + *(p22+2));
  edge_u2_v3 = ( *p23     + *(p23+1) + *(p23+2));
  edge_u3_v1 = ( *p31     + *(p31+1) + *(p31+2));
  edge_u3_v2 = ( *p32     + *(p32+1) + *(p32+2));
  edge_u3_v3 = ( *p33     + *(p33+1) + *(p33+2));
  
				/* lines varying along v */
  edge_u1_w1 = (  *p11    +  *p12    +  *p13   );
  edge_u1_w2 = ( *(p11+1) + *(p12+1) + *(p13+1));
  edge_u1_w3 = ( *(p11+2) + *(p12+2) + *(p13+2));
  edge_u2_w1 = (  *p21    +  *p22    +  *p23   );
  edge_u2_w2 = ( *(p21+1) + *(p22+1) + *(p23+1));
  edge_u2_w3 = ( *(p21+2) + *(p22+2) + *(p23+2));
  edge_u3_w1 = (  *p31    +  *p32    +  *p33   );
  edge_u3_w2 = ( *(p31+1) + *(p32+1) + *(p33+1));
  edge_u3_w3 = ( *(p31+2) + *(p32+2) + *(p33+2));
  
				/* lines varying along u */
  edge_v1_w1 = (  *p11    +  *p21    +  *p31   );
  edge_v1_w2 = ( *(p11+1) + *(p21+1) + *(p31+1));
  edge_v1_w3 = ( *(p11+2) + *(p21+2) + *(p31+2));
  edge_v2_w1 = (  *p12    +  *p22    +  *p32   );
  edge_v2_w2 = ( *(p12+1) + *(p22+1) + *(p32+1));
  edge_v2_w3 = ( *(p12+2) + *(p22+2) + *(p32+2));
  edge_v3_w1 = (  *p13    +  *p23    +  *p33   );
  edge_v3_w2 = ( *(p13+1) + *(p23+1) + *(p33+1));
  edge_v3_w3 = ( *(p13+2) + *(p23+2) + *(p33+2));
  
  slice_u1 =  (edge_u1_w1 + edge_u1_w2 + edge_u1_w3);
  slice_u2 =  (edge_u2_w1 + edge_u2_w2 + edge_u2_w3);
  slice_u3 =  (edge_u3_w1 + edge_u3_w2 + edge_u3_w3);
  slice_v1 =  (edge_v1_w1 + edge_v1_w2 + edge_v1_w3);
  slice_v2 =  (edge_v2_w1 + edge_v2_w2 + edge_v2_w3);
  slice_v3 =  (edge_v3_w1 + edge_v3_w2 + edge_v3_w3);
  slice_w1 =  (edge_u1_w1 + edge_u2_w1 + edge_u3_w1);
  slice_w2 =  (edge_u1_w2 + edge_u2_w2 + edge_u3_w2);
  slice_w3 =  (edge_u1_w3 + edge_u2_w3 + edge_u3_w3);
  
  d->u  = (slice_u3 - slice_u1) / 18.0;                          
  d->v  = (slice_v3 - slice_v1) / 18.0;                          
  d->w  = (slice_w3 - slice_w1) / 18.0;                          
  d->uu = (slice_u3 + slice_u1 - 2*slice_u2) / 9.0;                   
  d->vv = (slice_v3 + slice_v1 - 2*slice_v2) / 9.0;                   
  d->ww = (slice_w3 + slice_w1 - 2*slice_w2) / 9.0;                   
  d->uv = (edge_u3_v3 + edge_u1_v1 - edge_u3_v1 - edge_u1_v3) / 12.0;
  d->uw = (edge_u3_w3 + edge_u1_w1 - edge_u3_w1 - edge_u1_w3) / 12.0;  
  d->vw = (edge_v3_w3 + edge_v1_w1 - edge_v3_w1 - edge_v1_w3) / 12.0;  
  
} /* estimate_3D_derivatives */


private void estimate_2D_derivatives(Real r[3][3], 
				     deriv_2D_struct *d)

{
  Real	p11, p12, p13;
  Real	p21, p22, p23;
  Real	p31, p32, p33;

  /* --- 3x3 ---   [u][v]   */
  
  p11 = r[0][0];   p21 = r[1][0];   p31 = r[2][0]; 
  p12 = r[0][1];   p22 = r[1][1];   p32 = r[2][1]; 
  p13 = r[0][2];   p23 = r[1][2];   p33 = r[2][2]; 


/*
   old (pre 6.8.95) working weights:
  d->u  = (p31 + p32 + p33 - p11 - p12 - p13) / 2.0;                          
  d->v  = (p13 + p23 + p33 - p11 - p21 - p31) / 2.0;                          
  d->uu = (p31 + p32 + p33 + p11 + p12 + p13 - 2*(p21 + p22 + p23))/2.0;
  d->vv = (p13 + p23 + p33 + p11 + p21 + p31 - 2*(p12 + p22 + p32))/2.0;
  d->uv = (p11 + p33 - p13 - p31 ) /  2.0;  
*/

  d->u  = (p31 + p32 + p33 - p11 - p12 - p13) / 6.0;                          
  d->v  = (p13 + p23 + p33 - p11 - p21 - p31) / 6.0;                          
  d->uu = (p31 + p32 + p33 + p11 + p12 + p13 - 2*(p21 + p22 + p23))/3.0;
  d->vv = (p13 + p23 + p33 + p11 + p21 + p31 - 2*(p12 + p22 + p32))/3.0;
  d->uv = (p11 + p33 - p13 - p31 ) /  4.0;  

} /* estimate_2D_derivatives */



/* 
   negative_3D_definite()

   if the 3D quadratic represented by deriv_struc is negative definite
   then      return TRUE
   otherwise return FALSE.

   strang p 338, if A is neg def, then (-A) if pos def.
*/
private BOOLEAN negative_3D_definite(deriv_3D_struct *c)

{

  deriv_3D_struct d;

  d.uu = -1.0 * c->uu;
  d.vv = -1.0 * c->vv;
  d.ww = -1.0 * c->ww;
  d.uv = -1.0 * c->uv;
  d.uw = -1.0 * c->uw;
  d.vw = -1.0 * c->vw;
  

  return ( positive_3D_definite( &d ) );

}

/*    /          \                                                     
      | uu uv uw |   the matrix A is positive definite is all upper-   
   A= | uv vv vw |   left submatrices have positive determinants        
      | uw vw ww |                                                     
      \          /            strang, p 331                            */

private BOOLEAN positive_3D_definite(deriv_3D_struct *c)
{
  return ((c->uu > 0.0) &&

	  ((c->uu*c->vv - c->uv*c->uv) > 0.0 ) &&

	  ((c->uu * (c->vv*c->ww - c->vw*c->vw) -            
	    c->uv * (c->uv*c->ww - c->vw*c->uw) + 	       
	    c->uw * (c->uv*c->vw - c->vv*c->uw)) >0.0)
	  );
}

private BOOLEAN negative_2D_definite(deriv_2D_struct *c)

{

  return ((c->uu < 0) &&
	  ((c->uu*c->vv - c->uv*c->uv) > 0)
	  );

}

    


