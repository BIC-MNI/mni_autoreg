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
@MODIFIED   : Revision 96.7  2006-11-30 09:07:32  rotor
@MODIFIED   :  * many more changes for clean minc 2.0 build
@MODIFIED   :
@MODIFIED   : Revision 96.6  2006/11/29 09:09:33  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 96.5  2005/07/20 20:45:49  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 96.4  2004/02/13 00:16:48  rotor
@MODIFIED   :  * removed /static defs
@MODIFIED   :
@MODIFIED   : Revision 96.3  2002/03/26 14:15:41  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 96.2  2000/03/15 08:42:43  stever
@MODIFIED   : Code cleanup: all functions prototyped (except ParseArgs.c), no useless declarations, etc
@MODIFIED   :
@MODIFIED   : Revision 96.1  1997/11/03 19:59:49  louis
@MODIFIED   : - now include volume_io/internal_volume_io.h instead of volume_io.h
@MODIFIED   : - now include quad_max_fit.h for prototypes and struc defs
@MODIFIED   :
 * Revision 96.0  1996/08/21  18:21:58  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:21:54  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:15:50  louis
 * Never released version 0.95
 *
 * Revision 1.4  1996/08/12  14:15:47  louis
 * Pre-release
 *
 * Revision 1.3  1995/09/11  12:37:16  collins
 * this is a working version. While this file is not used in mni_reg-0.1g,
 * it will be part of the non-linear release.
 *
 * Revision 1.2  1995/08/21  11:31:14  collins
 * This version of the quadratic fit seems to work reasonably well in 3D
 * and compares favorably to the simplex optimization for hoge/jacob fitting
 * at 16mm.
 *
 * The use of quadratic fitting is twice as fast as simplex (6.4 vs 15.2 min)
 * for the 16mm fit (again with {hoge,jacob}_16_dxyz.mnc and -step 8 8 8.
 *
 * Revision 1.1  1995/06/12  14:30:23  collins
 * Initial revision
 *

---------------------------------------------------------------------------- */

#include <volume_io.h>                
#include <math.h>
#include <quad_max_fit.h>

#define SMALL_EPS 0.000000001

#define MINIMUM_DET_ALLOWED 0.00000001

extern int stat_quad_total;
extern int stat_quad_zero;
extern int stat_quad_two;
extern int stat_quad_plus;
extern int stat_quad_minus;
extern int stat_quad_semi;

    /* local prototypes */

static VIO_BOOL negative_2D_definite(deriv_2D_struct *c);
static VIO_BOOL negative_3D_definite(deriv_3D_struct *c);
static VIO_BOOL positive_3D_semidefinite(deriv_3D_struct *c);
static VIO_BOOL positive_3D_definite(deriv_3D_struct *c);
static void    estimate_2D_derivatives(VIO_Real r[3][3], 
                                        deriv_2D_struct *c);             

/* this procedure will return TRUE with the dx,dy,dz (offsets) that 
   correspond to the MAXIMUM value of the quadratic function fit through 
   data points represented in r[u][v][w], unless the maxtrix represented in 
   r[u][v][w] is negative definite
   in which case, the function will return FALSE.                         */

VIO_BOOL return_3D_disp_from_quad_fit(VIO_Real r[3][3][3], 
                                            VIO_Real *dispu, 
                                            VIO_Real *dispv, 
                                            VIO_Real *dispw)        
{
  deriv_3D_struct 
    d;                        /* the 1st and second order derivatives */
  VIO_Real
    du,dv,dw,
    a[3][3],
    detA;
  int 
    try_2d_def,
    count_u, count_v, count_w,
    i,j,k;
  char
    res;

  try_2d_def = FALSE;
  *dispu = *dispv = *dispw = 0.0;
  estimate_3D_derivatives(r,&d);

  /*    /          \    the values that form the matrix A come from the  */
  /*    | uu uv uw |    second order derivatives in 'd'.  If this matrix */
  /* A= | uv vv vw |    is not negative definite, then we can not find   */  
  /*    | uw vw ww |    a maximum in the region defined by r[][][].      */
  /*    \          /                                                     */



  if ( !negative_3D_definite(&d) )  {
    try_2d_def = TRUE; res = '+'; /* no 3D disp can be calculated */
  }
  else {                        /* otherwise, I have the derivatives, now get the 
                                   displacements that correspond to a maximum
                                   in r[][][]                                     */

    detA = d.uu * (d.vv*d.ww - d.vw*d.vw) -            
           d.uv * (d.uv*d.ww - d.vw*d.uw) +                
           d.uw * (d.uv*d.vw - d.vv*d.uw) ;               
                                                       
    if ( fabs( detA ) <= MINIMUM_DET_ALLOWED ) {
      try_2d_def = TRUE; res = '0';
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

      if ( fabs( *dispu ) < 2.0 && fabs( *dispv ) < 2.0 && fabs( *dispw ) < 2.0 ) {        
        print ("!"); return (TRUE);
      }
      else {
        try_2d_def = TRUE; res = '2';
      }
    }
  }

        
  if (try_2d_def) {
    
    *dispu = *dispv = *dispw = 0.0;
    count_u = count_v = count_w = 0;
    
    for(i=0; i<3; i++)
      for(j=0; j<3; j++)
        a[i][j] = r[1][i][j];
    
    if (return_2D_disp_from_quad_fit(a, &dv, &dw) ) {
      *dispv += dv;
      *dispw += dw;
      count_v++;
      count_w++;
    }
    
    for(i=0; i<3; i++)
      for(j=0; j<3; j++)
        a[i][j] = r[i][1][j];
    
    if (return_2D_disp_from_quad_fit(a, &du, &dw) ) {
      *dispu += du;
      *dispw += dw;
      count_u++;
      count_w++;
    }
    
    for(i=0; i<3; i++)
      for(j=0; j<3; j++)
        a[i][j] = r[i][j][1];
    
    if (return_2D_disp_from_quad_fit(a, &du, &dv) ) {
      *dispu += du;
      *dispv += dv;
      count_u++;
      count_v++;
    }
    
    if (count_u == 0 && count_v == 0 && count_w ==0) {
/*      print ("\n%8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f\n",
             r[0][0][0],r[1][0][0],r[2][0][0], 
             r[0][0][1],r[1][0][1],r[2][0][1], r[0][0][2],r[1][0][2],r[2][0][2]);
      print ("%8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f\n",
             r[0][1][0],r[1][1][0],r[2][1][0], 
             r[0][1][1],r[1][1][1],r[2][1][1], r[0][1][2],r[1][1][2],r[2][1][2]);
      print ("%8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f\n",
             r[0][2][0],r[1][2][0],r[2][2][0], 
             r[0][2][1],r[1][2][1],r[2][2][1], r[0][2][2],r[1][2][2],r[2][2][2]);
      
      print ("/%8.5f %8.5f %8.5f\\  /%8.5f %8.5f %8.5f\\  /%8.5f\\  ->  du = %8.5f\n",  
             d.uu, d.uv, d.uw, a[0][0],a[0][1],a[0][2], d.u, *dispu);
      print ("|%8.5f %8.5f %8.5f|  |%8.5f %8.5f %8.5f|  |%8.5f|  ->  dv = %8.5f\n",  
             d.uv, d.vv, d.vw, a[1][0],a[1][1],a[1][2], d.v, *dispv);
      print ("\\%8.5f %8.5f %8.5f/  \\%8.5f %8.5f %8.5f/  \\%8.5f/  ->  dw = %8.5f\n\n",
             d.uw, d.vw, d.ww, a[2][0],a[2][1],a[2][2], d.w, *dispw);
  */

      print ("%c",res);        
                                /* In this case, no deformation can be
                                   estimated when using the 2D quad
                                   fit.  So, just search out the
                                   maximum value within the 3x3 array,
                                   and return half the displacement
                                   needed to get there.*/

      count_u=1;count_v=1;count_w=1;
      for(i=0; i<3; i++)
        for(j=0; j<3; j++)
          for(k=0; k<3; k++)
            if (r[i][j][k]>r[count_u][count_v][count_w]) {
              count_u=i; count_v=j; count_w=k;
            }

      *dispu = (count_u - 1.0) /2.0;
      *dispv = (count_v - 1.0) /2.0;
      *dispw = (count_w - 1.0) /2.0;

      return(TRUE);
    }
    else {
      if ( !(count_u == 0)) *dispu /= count_u ;
      if ( !(count_v == 0)) *dispv /= count_v ;
      if ( !(count_w == 0)) *dispw /= count_w ;
      
      print ("*");
      
      return(TRUE);
    }
  }
  else
    return(FALSE);
}


/* this procedure will return TRUE with the dx,dy,dz (offsets) that
   correspond to the MINIMUM of the quadratic function fit through
   data points represented in r[u][v][w], which must be at least
   positve semi-definite.  
   if positive-definite, 
       then the inverse of the covariance matrix is used to
       directly calculate the MINIMUM of the quadratic
       and return TRUE
   if positive semi-definite,
       then find the MINIMUM only along the directions that 
       have NON-NULL eigenvectors, and return TRUE
   otherwise,
       return half the distance required to get the minimum value
       represented in r[][][], and return FALSE */


VIO_BOOL return_3D_disp_from_min_quad_fit(VIO_Real r[3][3][3], 
                                                VIO_Real *dispu, 
                                                VIO_Real *dispv, 
                                                VIO_Real *dispw)        
{
  deriv_3D_struct 
    d;                        /* the 1st and second order derivatives */
  VIO_Real
    a[3][3],
    detA;
  int 
    hack_to_min_val,
    count_u, count_v, count_w,
    i,j,k;

  *dispu = *dispv = *dispw = 0.0;
  estimate_3D_derivatives_new(r,&d);
  hack_to_min_val = FALSE;
  stat_quad_total++;

  /*    /          \    the values that form the matrix A come from the  */
  /*    | uu uv uw |    second order derivatives in 'd'.  If this matrix */
  /* A= | uv vv vw |    is not positive definite, then we can not find   */  
  /*    | uw vw ww |    a minimum in the region defined by r[][][].      */
  /*    \          /                                                     */

  if (positive_3D_definite(&d)) {
                        /* I have the derivatives, now get the 
                           displacements that correspond to a minimum
                           in r[][][]                                     */

    detA = d.uu * (d.vv*d.ww - d.vw*d.vw) -            
           d.uv * (d.uv*d.ww - d.vw*d.uw) +                
           d.uw * (d.uv*d.vw - d.vv*d.uw) ;               
                                                       
    if ( fabs( detA ) <= MINIMUM_DET_ALLOWED ) {
      hack_to_min_val = TRUE; stat_quad_zero++;
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

      if ( fabs( *dispu ) < 2.0 && fabs( *dispv ) < 2.0 && fabs( *dispw ) < 2.0 ) {        
        stat_quad_plus++; return (TRUE);
      }
      else {
        hack_to_min_val = TRUE; stat_quad_two++;    
      }
    }
  }
  else {
    if (positive_3D_semidefinite(&d)) {
      stat_quad_semi++;
      /* get_disp_from_positive_semidefinite(&d, dispu, dispv, dispw); */
    }
    else{
      stat_quad_minus++;;
      hack_to_min_val = TRUE;    
    }
  }

  hack_to_min_val = FALSE;

  if (hack_to_min_val) {
    
    *dispu = *dispv = *dispw = 0.0;
    count_u = count_v = count_w = 0;
    
    count_u=1;count_v=1;count_w=1;
    for(i=0; i<3; i++)
      for(j=0; j<3; j++)
        for(k=0; k<3; k++)
          if (r[i][j][k]<r[count_u][count_v][count_w]) {
            count_u=i; count_v=j; count_w=k;
          }
    
    *dispu = (count_u - 1.0) /2.0;
    *dispv = (count_v - 1.0) /2.0;
    *dispw = (count_w - 1.0) /2.0;
    
    return(TRUE);
  }

  return(FALSE);
}

VIO_BOOL return_2D_disp_from_quad_fit(VIO_Real r[3][3], /* the values used in the quad fit */
                                            VIO_Real *dispu, /* the displacements returned */
                                            VIO_Real *dispv)        
{
  deriv_2D_struct 
    d;                        /* the 1st and second order derivatives */
  VIO_Real
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
    return(FALSE);                /* no disp can be calculated */

  }
  else {                        /* otherwise, I have the derivatives, now get 
                                   the displacements that correspond to a 
                                   maximum in r[][]                        */

    detA = d.uu * d.vv - d.uv*d.uv;
                                                       
    if (fabs( detA) < MINIMUM_DET_ALLOWED) {

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

      if ( fabs( *dispu ) > 1.5 || fabs( *dispv ) > 1.5 ) {
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

void estimate_3D_derivatives(VIO_Real r[3][3][3], 
                                    deriv_3D_struct *d) 

{

  VIO_Real        *p11, *p12, *p13;
  VIO_Real        *p21, *p22, *p23;
  VIO_Real        *p31, *p32, *p33;

  VIO_Real        slice_u1, slice_u2, slice_u3;
  VIO_Real        slice_v1, slice_v2, slice_v3;
  VIO_Real        slice_w1, slice_w2, slice_w3;
  
  VIO_Real        edge_u1_v1, /* edge_u1_v2,*/ edge_u1_v3;
/*  VIO_Real        edge_u2_v1, edge_u2_v2, edge_u2_v3; */
  VIO_Real        edge_u3_v1,  /* edge_u3_v2,*/ edge_u3_v3;
  VIO_Real        edge_u1_w1, edge_u1_w2, edge_u1_w3;
  VIO_Real        edge_u2_w1, edge_u2_w2, edge_u2_w3;
  VIO_Real        edge_u3_w1, edge_u3_w2, edge_u3_w3;
  VIO_Real        edge_v1_w1, edge_v1_w2, edge_v1_w3;
  VIO_Real        edge_v2_w1, edge_v2_w2, edge_v2_w3;
  VIO_Real        edge_v3_w1, edge_v3_w2, edge_v3_w3;
  
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
/*  edge_u1_v2 = ( *p12     + *(p12+1) + *(p12+2)); */
  edge_u1_v3 = ( *p13     + *(p13+1) + *(p13+2));
/*  edge_u2_v1 = ( *p21     + *(p21+1) + *(p21+2)); */
/*  edge_u2_v2 = ( *p22     + *(p22+1) + *(p22+2)); */
/*  edge_u2_v3 = ( *p23     + *(p23+1) + *(p23+2)); */
  edge_u3_v1 = ( *p31     + *(p31+1) + *(p31+2));
/*  edge_u3_v2 = ( *p32     + *(p32+1) + *(p32+2)); */
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

#define INV_SQRT2 1.0
/*0.70710678*/
#define INV_SQRT3 1.0
/*0.57735027*/
void estimate_3D_derivatives_weighted(VIO_Real r[3][3][3], 
                                             deriv_3D_struct *d) 

{

  VIO_Real        *p11, *p12, *p13;
  VIO_Real        *p21, *p22, *p23;
  VIO_Real        *p31, *p32, *p33;

  VIO_Real        slice_u1, slice_u2, slice_u3;
  VIO_Real        slice_v1, slice_v2, slice_v3;
  VIO_Real        slice_w1, slice_w2, slice_w3;
  
  VIO_Real        edge_u1_v1, /* edge_u1_v2,*/ edge_u1_v3;
/*  VIO_Real        edge_u2_v1, edge_u2_v2, edge_u2_v3; */
  VIO_Real        edge_u3_v1,  /* edge_u3_v2,*/ edge_u3_v3;
  VIO_Real        edge_u1_w1, edge_u1_w2, edge_u1_w3;
  VIO_Real        edge_u2_w1, edge_u2_w2, edge_u2_w3;
  VIO_Real        edge_u3_w1, edge_u3_w2, edge_u3_w3;
  VIO_Real        edge_v1_w1, edge_v1_w2, edge_v1_w3;
  VIO_Real        edge_v2_w1, edge_v2_w2, edge_v2_w3;
  VIO_Real        edge_v3_w1, edge_v3_w2, edge_v3_w3;
  
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
  edge_u1_v1 = ( INV_SQRT3 * *p11     + INV_SQRT2 * *(p11+1) + INV_SQRT3 * *(p11+2));
  edge_u1_v3 = ( INV_SQRT3 * *p13     + INV_SQRT2 * *(p13+1) + INV_SQRT3 * *(p13+2));
  edge_u3_v1 = ( INV_SQRT3 * *p31     + INV_SQRT2 * *(p31+1) + INV_SQRT3 * *(p31+2));
  edge_u3_v3 = ( INV_SQRT3 * *p33     + INV_SQRT2 * *(p33+1) + INV_SQRT3 * *(p33+2));
  
                                /* lines varying along v */
  edge_u1_w1 = ( INV_SQRT3 * *p11    +  INV_SQRT2 * *p12    +  INV_SQRT3 * *p13   );
  edge_u1_w2 = ( INV_SQRT3 * *(p11+1) + INV_SQRT2 * *(p12+1) + INV_SQRT3 * *(p13+1));
  edge_u1_w3 = ( INV_SQRT3 * *(p11+2) + INV_SQRT2 * *(p12+2) + INV_SQRT3 * *(p13+2));
  edge_u2_w1 = ( INV_SQRT2 * *p21    +             *p22    +  INV_SQRT2 * *p23   );
  edge_u2_w2 = ( INV_SQRT2 * *(p21+1) +            *(p22+1) + INV_SQRT2 * *(p23+1));
  edge_u2_w3 = ( INV_SQRT2 * *(p21+2) +            *(p22+2) + INV_SQRT2 * *(p23+2));
  edge_u3_w1 = ( INV_SQRT3 * *p31    +  INV_SQRT2 * *p32    +  INV_SQRT3 * *p33   );
  edge_u3_w2 = ( INV_SQRT3 * *(p31+1) + INV_SQRT2 * *(p32+1) + INV_SQRT3 * *(p33+1)); 
  edge_u3_w3 = ( INV_SQRT3 * *(p31+2) + INV_SQRT2 * *(p32+2) + INV_SQRT3 * *(p33+2));
  
                                /* lines varying along u */
  edge_v1_w1 = ( INV_SQRT3 * *p11    +  INV_SQRT2 * *p21    +  INV_SQRT3 * *p31   );
  edge_v1_w2 = ( INV_SQRT3 * *(p11+1) + INV_SQRT2 * *(p21+1) + INV_SQRT3 * *(p31+1));
  edge_v1_w3 = ( INV_SQRT3 * *(p11+2) + INV_SQRT2 * *(p21+2) + INV_SQRT3 * *(p31+2));
  edge_v2_w1 = ( INV_SQRT2 * *p12    +             *p22    +  INV_SQRT2 * *p32   );
  edge_v2_w2 = ( INV_SQRT2 * *(p12+1) +            *(p22+1) + INV_SQRT2 * *(p32+1));
  edge_v2_w3 = ( INV_SQRT2 * *(p12+2) +            *(p22+2) + INV_SQRT2 * *(p32+2));
  edge_v3_w1 = ( INV_SQRT3 * *p13    +  INV_SQRT2 * *p23    +  INV_SQRT3 * *p33   );
  edge_v3_w2 = ( INV_SQRT3 * *(p13+1) + INV_SQRT2 * *(p23+1) + INV_SQRT3 * *(p33+1));
  edge_v3_w3 = ( INV_SQRT3 * *(p13+2) + INV_SQRT2 * *(p23+2) + INV_SQRT3 * *(p33+2));
  
  slice_u1 =  (edge_u1_w1 + edge_u1_w2 + edge_u1_w3);
  slice_u2 =  (edge_u2_w1 + edge_u2_w2 + edge_u2_w3);
  slice_u3 =  (edge_u3_w1 + edge_u3_w2 + edge_u3_w3);
  slice_v1 =  (edge_v1_w1 + edge_v1_w2 + edge_v1_w3);
  slice_v2 =  (edge_v2_w1 + edge_v2_w2 + edge_v2_w3);
  slice_v3 =  (edge_v3_w1 + edge_v3_w2 + edge_v3_w3);
  slice_w1 =  (edge_u1_w1 + edge_u2_w1 + edge_u3_w1);
  slice_w2 =  (edge_u1_w2 + edge_u2_w2 + edge_u3_w2);
  slice_w3 =  (edge_u1_w3 + edge_u2_w3 + edge_u3_w3);
  
  d->u  = (slice_u3 - slice_u1) / ((4*(INV_SQRT3 + INV_SQRT2) +1)*2.0);
  d->v  = (slice_v3 - slice_v1) / ((4*(INV_SQRT3 + INV_SQRT2) +1)*2.0);
  d->w  = (slice_w3 - slice_w1) / ((4*(INV_SQRT3 + INV_SQRT2) +1)*2.0);
  d->uu = (slice_u3 + slice_u1 - 2*slice_u2) / (4*(INV_SQRT3 + INV_SQRT2) +1);
  d->vv = (slice_v3 + slice_v1 - 2*slice_v2) / (4*(INV_SQRT3 + INV_SQRT2) +1);
  d->ww = (slice_w3 + slice_w1 - 2*slice_w2) / (4*(INV_SQRT3 + INV_SQRT2) +1);
  d->uv = (edge_u3_v3 + edge_u1_v1 - edge_u3_v1 - edge_u1_v3) / (4.0 * (2.0*INV_SQRT3 + INV_SQRT2));
  d->uw = (edge_u3_w3 + edge_u1_w1 - edge_u3_w1 - edge_u1_w3) / (4.0 * (2.0*INV_SQRT3 + INV_SQRT2));
  d->vw = (edge_v3_w3 + edge_v1_w1 - edge_v3_w1 - edge_v1_w3) / (4.0 * (2.0*INV_SQRT3 + INV_SQRT2));


/*

  d->u  = (slice_u3 - slice_u1) /  18.0;
  d->v  = (slice_v3 - slice_v1) /  18.0;
  d->w  = (slice_w3 - slice_w1) /  18.0;
  d->uu = (slice_u3 + slice_u1 - 2*slice_u2) / 9.0; 
  d->vv = (slice_v3 + slice_v1 - 2*slice_v2) / 9.0; 
  d->ww = (slice_w3 + slice_w1 - 2*slice_w2) / 9.0; 
  d->uv = (edge_u3_v3 + edge_u1_v1 - edge_u3_v1 - edge_u1_v3) / 12.0;
  d->uw = (edge_u3_w3 + edge_u1_w1 - edge_u3_w1 - edge_u1_w3) / 12.0;
  d->vw = (edge_v3_w3 + edge_v1_w1 - edge_v3_w1 - edge_v1_w3) / 12.0;
*/
  
} /* estimate_3D_derivatives */
#undef INV_SQRT2 
#undef INV_SQRT3 


/* the procedure above 'smooths' the derivative estimates,
   the procedure here will use the minimum amount of info
   to estimate the derivatives */
void estimate_3D_derivatives_new(VIO_Real r[3][3][3], 
                                        deriv_3D_struct *d) 

{

  d->u  = (r[2][1][1] - r[0][1][1] ) / 2.0;
  d->v  = (r[1][2][1] - r[1][0][1] ) / 2.0;
  d->w  = (r[1][1][2] - r[1][1][0] ) / 2.0;
  d->uu = (r[2][1][1] + r[0][1][1] -2*r[1][1][1]);
  d->vv = (r[1][2][1] + r[1][0][1] -2*r[1][1][1]);
  d->ww = (r[1][1][2] + r[1][1][0] -2*r[1][1][1]);
  d->uv = (r[2][2][1] + r[0][0][1] - r[0][2][1] - r[2][0][1]) / 4.0;
  d->uw = (r[2][1][2] + r[0][1][0] - r[0][1][2] - r[2][1][0]) / 4.0;
  d->vw = (r[1][2][2] + r[1][0][0] - r[1][0][2] - r[1][2][0]) / 4.0;
  
} /* estimate_3D_derivatives */



static void estimate_2D_derivatives(VIO_Real r[3][3], 
                                     deriv_2D_struct *d)

{
  VIO_Real        p11, p12, p13;
  VIO_Real        p21, p22, p23;
  VIO_Real        p31, p32, p33;

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
static VIO_BOOL negative_3D_definite(deriv_3D_struct *c)

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

static VIO_BOOL positive_3D_definite(deriv_3D_struct *c)
{
  return ((c->uu > 0.0) &&

          ((c->uu*c->vv - c->uv*c->uv) > 0.0 ) &&

          ((c->uu * (c->vv*c->ww - c->vw*c->vw) -            
            c->uv * (c->uv*c->ww - c->vw*c->uw) +                
            c->uw * (c->uv*c->vw - c->vv*c->uw)) >0.0)
          );
}

static VIO_BOOL positive_3D_semidefinite(deriv_3D_struct *c)
{

  return ((c->uu > -SMALL_EPS) &&
          (c->vv > -SMALL_EPS) &&
          (c->ww > -SMALL_EPS) &&

          ((c->uu*c->vv - c->uv*c->uv) > -SMALL_EPS ) &&
          ((c->uu*c->ww - c->uw*c->uw) > -SMALL_EPS ) &&
          ((c->vv*c->ww - c->vw*c->vw) > -SMALL_EPS ) &&

          ((c->uu * (c->vv*c->ww - c->vw*c->vw) -            
            c->uv * (c->uv*c->ww - c->vw*c->uw) +                
            c->uw * (c->uv*c->vw - c->vv*c->uw)) > -SMALL_EPS) 

          );
}

static VIO_BOOL negative_2D_definite(deriv_2D_struct *c)

{

  return ((c->uu < 0) &&
          ((c->uu*c->vv - c->uv*c->uv) > 0)
          );

}

    
/*
   this procedure is based on the paper of Thirion and Gourdon,
   "Computing the differential characteristics of isointensity surfaces"
   CVIU 61(2) march 190-202, 1995

   The goal is to return the differential charateristics of an
   iso-surface passing through r[1][1][1] contained within a volume
   represented in r[][][].

   The characteristics returned:

       dir_1[], 
       dir_2[] - the vector directions of min amd max curvature
       r_K     - The scalar value of Gaussian curvature
       r_S     - The scalar value of the mean curvature
       r_k1    - The scalar value of the max curvature
       r_k2    - The scalar value of the min curvature
       r_norm[]- the vector normal to the isosurf at r[1][1][1] 
       Lvv     - is actually the mean Lvv value.

   bugs: if the normal to the iso surface has u and v componants
         equal ( du=dv ) then
         the procedure fails to find dir_1 and dir_2.

         Also, the authors state that there is a preferred direction
         with grad(f) = (1,1,1).  In this case, beta=(0,0,0) and
         alpha = (0,0,0) -> so no direction vectors found.

         In either case, the Gaussian, mean, maximum and minimum
         curvature values are found correctly.  The problems concern
         the curvature directions only.  */


VIO_BOOL return_principal_directions(VIO_Real r[3][3][3],
                                           VIO_Real dir_1[3],
                                           VIO_Real dir_2[3],
                                           VIO_Real *r_K,
                                           VIO_Real *r_S,
                                           VIO_Real *r_k1,
                                           VIO_Real *r_k2,
                                           VIO_Real *r_norm,
                                           VIO_Real *r_Lvv,
                                           VIO_Real eps)

{
  deriv_3D_struct 
    d;                        

  VIO_Real
    K,                                /* Gaussian curvature      */
    S,                                /* mean curvature          */
    k1,                                /* value of max curvature  */
    k2,                                /* value of min curvature  */
    Lvv,
    tmp,
    det,sq_det,                        /* determinant             */
    len1,len2,                        /* length of vector        */
    sq_mag_grad,                /* square of gradient mag  */
    beta[3],                        /* the vector beta         */
    alpha[3],                        /* the vector alpha        */
    x,y,z,                        /* first order derivatives */
    xx,yy,zz,xy,xz,yz;                /* second order derivative */
  int 
    i;


  estimate_3D_derivatives(r,&d);

  x  = d.u;                        /* for notational simplicity below */
  y  = d.v;
  z  = d.w;
  xx = d.uu;
  yy = d.vv;
  zz = d.ww;
  xy = d.uv;
  yz = d.vw;
  xz = d.uw;

  print ("in return_principal_directions() covar:\n");
  print ("%12.8f %12.8f %12.8f \n",xx,xy,xz);
  print ("%12.8f %12.8f %12.8f \n",xy,yy,yz);
  print ("%12.8f %12.8f %12.8f \n",xz,yz,zz);

  sq_mag_grad = x*x + y*y + z*z;

  *r_K = *r_S = *r_k1 = *r_k2 = 0.0;

  if (r_norm != NULL) {
    r_norm[VIO_X] = x ;
    r_norm[VIO_Y] = y ;
    r_norm[VIO_Z] = z ;
  }
  
  if ( fabs(sq_mag_grad)<eps ) 
    return(FALSE);
    
                                /* Gaussian curvature: */
  K = (
       x*x*(yy*zz - yz*yz) + 2*y*z*(xz*xy - xx*yz) +
       y*y*(xx*zz - xz*xz) + 2*x*z*(yz*xy - yy*xz) +
       z*z*(xx*yy - xy*xy) + 2*x*y*(xz*yz - zz*xy)
      )
       /(sq_mag_grad*sq_mag_grad);

                                /* Mean curvature */
  S = (
       x*x*(yy + zz) - 2*y*z*yz +
       y*y*(xx + zz) - 2*x*z*xz +
       z*z*(xx + yy) - 2*x*y*xy
      )
       / (2 * sqrt(sq_mag_grad*sq_mag_grad*sq_mag_grad));


  det = S*S - K; 

  if ( fabs(det) < SMALL_EPS ) det = 0.0;

  if (det<0.0) {
    print ("det (S*S - K) is negative, and this shouldn't be...\n");
    det = fabs(det);
  }

  sq_det = sqrt(det);

                                /* min and max curvatures */
  k1 = S + sq_det;
  k2 = S - sq_det;

  Lvv =  sq_mag_grad * S;

  print ("K, S, k1, k2, Lvv, sq_det:   %f %f %f %f %f\n",K, S, k1, k2, Lvv,sq_det);

  len1 = len2 = 0.0;

  if (dir_1  != NULL && dir_2 != NULL) {
    /* calc principal directions */
    beta[VIO_X] = z-y;
    beta[VIO_Y] = x-z;
    beta[VIO_Z] = y-x;
    
    alpha[VIO_X] = (-2.0*z*z*z*xy + y*y*y*zz + 2.0*y*y*y*xz - 2.0*y*y*z*xy
                +2.0*z*z*x*yz  + 2.0*z*z*y*xz - 2.0*y*y*x*yz - 2.0*z*x*y*zz
                +2.0*x*y*z*yy + y*y*z*xx - 2.0*z*z*x*xz + z*x*x*zz
                -x*x*z*yy + 2.0*z*z*y*yz - z*y*y*zz +z*z*z*xx - z*z*z*yy
                -2.0*y*y*x*xz + 2.0*x*x*y*yz - y*y*y*xx + 2.0*x*z*z*xy - y*z*z*xx
                -2.0*z*y*y*yz + y*z*z*yy - 2.0*z*x*x*yz + 2.0*x*y*y*xy
                +x*x*y*zz - x*x*y*yy
                );
    
    alpha[VIO_Y] = (-2.0*x*x*x*yz + z*z*z*xx + 2.0*z*z*z*xy - 2.0*z*z*x*yz
                +2.0*x*x*y*xz  + 2.0*x*x*z*xy - 2.0*z*z*y*xz - 2.0*x*y*z*xx
                +2.0*y*z*x*zz + z*z*x*yy - 2.0*x*x*y*xy + x*y*y*xx
                -y*y*x*zz + 2.0*x*x*z*xz - x*z*z*xx +x*x*x*yy - x*x*x*zz
                -2.0*z*z*y*xy + 2.0*y*y*z*xz - z*z*z*yy + 2.0*y*x*x*yz - z*x*x*yy
                -2.0*x*z*z*xz + z*x*x*zz - 2.0*x*y*y*xz + 2.0*y*z*z*yz
                +y*y*z*xx - y*y*z*zz
                );
    
    alpha[VIO_Z] = (-2.0*y*y*y*xz + x*x*x*yy + 2.0*x*x*x*yz - 2.0*x*x*y*xz
                +2.0*y*y*z*xy  + 2.0*y*y*x*yz - 2.0*x*x*z*xy - 2.0*y*z*x*yy
                +2.0*z*x*y*xx + x*x*y*zz - 2.0*y*y*z*yz + y*z*z*yy
                -z*z*y*xx + 2.0*y*y*x*xy - y*x*x*yy +y*y*y*zz - y*y*y*xx
                -2.0*x*x*z*yz + 2.0*z*z*x*xy - x*x*x*zz + 2.0*z*y*y*xz - x*y*y*zz
                -2.0*y*x*x*xy + x*y*y*xx - 2.0*y*z*z*xy + 2.0*z*x*x*xz
                +z*z*x*yy - z*z*x*xx
                );
    
    for(i=VIO_X; i<=VIO_Z; i++)
      alpha[i] = alpha[i] / (2 * sqrt(sq_mag_grad*sq_mag_grad*sq_mag_grad));
    
    /* t_i = alpha +/- sqrt(det)*beta */
    
    for(i=VIO_X; i<=VIO_Z; i++)
      dir_1[i] = alpha[i] + sq_det*beta[i];
    
    for(i=VIO_X; i<=VIO_Z; i++)
      dir_2[i] = alpha[i] - sq_det*beta[i];
    

    print ("dir1: %12.8f %12.8f %12.8f\n", dir_1[0],dir_1[1],dir_1[2]);
    print ("dir2: %12.8f %12.8f %12.8f\n", dir_2[0],dir_2[1],dir_2[2]);
    print ("alph: %12.8f %12.8f %12.8f\n", alpha[0],alpha[1],alpha[2]);
    print ("beta: %12.8f %12.8f %12.8f\n", beta[0],beta[1],beta[2]);
    
    len1 = sqrt(dir_1[VIO_X]*dir_1[VIO_X] + dir_1[VIO_Y]*dir_1[VIO_Y] + dir_1[VIO_Z]*dir_1[VIO_Z]);
    for(i=VIO_X; i<=VIO_Z; i++)
      dir_1[i] /= len1;
    
    len2 = sqrt(dir_2[VIO_X]*dir_2[VIO_X] + dir_2[VIO_Y]*dir_2[VIO_Y] + dir_2[VIO_Z]*dir_2[VIO_Z]);
    for(i=VIO_X; i<=VIO_Z; i++)
      dir_2[i] /= len2;
  }


  if (fabs(k1)<fabs(k2)) {        /* ensure k1>k2  */
    
    /* swap curvatures */
    tmp= k1;    k1 = k2;    k2 = tmp;

    if (dir_1  != NULL && dir_2 != NULL) {
      /* swap min and max direction vectors */
      tmp = dir_1[VIO_X]; dir_1[VIO_X] = dir_2[VIO_X]; dir_2[VIO_X] = tmp;
      tmp = dir_1[VIO_Y]; dir_1[VIO_Y] = dir_2[VIO_Y]; dir_2[VIO_Y] = tmp;
      tmp = dir_1[VIO_Z]; dir_1[VIO_Z] = dir_2[VIO_Z]; dir_2[VIO_Z] = tmp;
    }
  }

                                /* set return vals */
  *r_k1 = k1;
  *r_k2 = k2;
  *r_K = K;
  *r_S = S;
  *r_Lvv = Lvv;

  if (dir_1  != NULL && dir_2 != NULL) {
    
                                /* ensure that directions found
                                   are perpendicular */

    tmp = dir_1[VIO_X]*dir_2[VIO_X] + dir_1[VIO_Y]*dir_2[VIO_Y] + dir_1[VIO_Z]*dir_2[VIO_Z];
    
    if ( fabs(tmp)>eps  && len1>3.0e-9 && len2>3.0e-9) {
      return(FALSE);
    }
  }


  return(TRUE);
}
        

/*
   The goal of this procedure is to return the differential
   charateristics of an iso-contour passing through r[1][1] 
   contained within an image region represented in r[][].

   The characteristics returned:

       r_tan[] - the tangent vector of the isocontour at r[1][1]
       r_norm[]- the normal vector of the isocontour at r[1][1]
       r_K     - the curvature of the isocontour at r[1][1]
*/


VIO_BOOL return_2D_principal_directions(VIO_Real r[3][3],
                                              VIO_Real norm[3],
                                              VIO_Real tang[3],
                                              VIO_Real *K,
                                              VIO_Real eps)

{
  deriv_2D_struct 
    d;                        

  VIO_Real
    sq_mag_grad,                /* square of gradient mag  */
    mag_grad,                        /* gradient mag            */
    x,y,                        /* first order derivatives */
    xx,yy,xy;                        /* second order derivative */

  *K = 0.0;

  estimate_2D_derivatives(r,&d);

  x  = d.u;                        /* for notational simplicity below */
  y  = d.v;
  xx = d.uu;
  yy = d.vv;
  xy = d.uv;

  sq_mag_grad = x*x + y*y; mag_grad = sqrt(sq_mag_grad);

  if ( fabs(sq_mag_grad)<eps )  {
    return(FALSE);
  }
  else {
    
    norm[VIO_X] = x/mag_grad;        /* vect normal to isocontour */
    norm[VIO_Y] = y/mag_grad; 
    norm[VIO_Z] = 0.0;

    tang[VIO_X] = -y/mag_grad;        /* vect tangent to isocontour */
    tang[VIO_Y] = x/mag_grad; 
    tang[VIO_Z] = 0.0;

                                /* curvature */

    *K = (2*x*y*xy - x*x*yy -xx*y*y) / sqrt(sq_mag_grad * sq_mag_grad * sq_mag_grad);

    return(TRUE);
  }
}
                
                        
VIO_Real return_Lvv(VIO_Real r[3][3][3],
                       VIO_Real eps)
     
{
  deriv_3D_struct 
    d;                        

  VIO_Real
    S,                                /* mean curvature          */
    Lvv,
    sq_mag_grad,                /* square of magnitude of gradient   */
    x,y,z,                        /* first order derivatives */
    xx,yy,zz,xy,xz,yz;                /* second order derivative */

  estimate_3D_derivatives_new(r,&d);

  x  = d.u;  y  = d.v;  z  = d.w;
  xx = d.uu; yy = d.vv; zz = d.ww;
  xy = d.uv; yz = d.vw; xz = d.uw;

  Lvv = 0.0;
  sq_mag_grad = x*x + y*y + z*z;

  if ( fabs(sq_mag_grad) > eps )  {
                                /* Mean curvature */
    S = (
         x*x*(yy + zz) - 2*y*z*yz +
         y*y*(xx + zz) - 2*x*z*xz +
         z*z*(xx + yy) - 2*x*y*xy
         )
          / (2 * sqrt(sq_mag_grad*sq_mag_grad*sq_mag_grad));

    Lvv =  sq_mag_grad * S;
  }

  return(Lvv);
}
                                        




/*
   this procedure will use principal component analysis to extract the
   eigen vectors and eigen values of the 3x3x3 objective function grid
   that represent the maximum, middle and minimum variation of the obj
   func.

   The procedure is based on information taken from "Probabilites et
   analyse des donnees et statistique" by G. Saporta (QA 273 SAP)
   p 159-186

   The charateristics returned are
       dir_1[], 
       dir_2[],
       dir_3[] - the vector directions of max, mid and min curvature,
                 (normalized vectors are returned)
       val[]   - the corresponding eigen values
*/


VIO_BOOL return_local_eigen(VIO_Real r[3][3][3],
                                  VIO_Real dir_1[3],
                                  VIO_Real dir_2[3],
                                  VIO_Real dir_3[3],
                                  VIO_Real val[3])

{
  int 
    eig_flag,iters,cnt,m,n,i,j,k;
  VIO_Real 
    **data, **weighted_data, **covar, **eigvec, *eigval;

  VIO_ALLOC2D(data,27,4);
  VIO_ALLOC2D(weighted_data,27,4);
  VIO_ALLOC2D(eigvec,3,3);
  VIO_ALLOC2D(covar,3,3);
  ALLOC(eigval,3);

  cnt = 0;
                                
  for(i=0; i<3; i++)                /* set up data matrix for easy manipulation */
    for(j=0; j<3; j++)
      for(k=0; k<3; k++){
        data[cnt][0] = i-1.0;
        data[cnt][1] = j-1.0;
        data[cnt][2] = k-1.0;
         data[cnt][3] = r[i][j][k];
        cnt++;
      }
                                /* covar = data[:,0:2]' * w * data[:,0:2],
                                      where w is stored in data[:,3]      */
  for(i=0; i<3; i++)
    for(j=0; j<27; j++) 
      weighted_data[j][i] = data[j][i] * data[j][3];

  for(m=0; m<3; m++) {
    for(n=0; n<3; n++) {
      covar[m][n] = 0.0;
      for(i=0; i<27; i++)                
        covar[m][n] += weighted_data[i][m] * data[i][n];

      if (fabs(covar[m][n]) < SMALL_EPS) covar[m][n] = 0.0;
    }
  }
  
                                  /* is the covariance matrix positive definite? */
  
/*
  flag = (covar[0][0]>0.0 &&
          (covar[0][0]*covar[1][1] - covar[0][1]*covar[1][0])>0 &&
          ((covar[0][0] * (covar[1][1]*covar[2][2] - covar[1][2]*covar[2][1])) -
           (covar[0][1] * (covar[1][0]*covar[2][2] - covar[1][2]*covar[2][0])) +
           (covar[0][2] * (covar[1][0]*covar[2][1] - covar[1][1]*covar[2][0]))
           ) > 0.0
          );

  if (!flag)
    print ("Not positive definite!\n");
*/

                                /* calculate eigen vectors/values, returning the 
                                   eigen vectors in column format within the matrix
                                   eigvec */
  eig_flag = eigen(covar, 3, eigval, eigvec, &iters);

  if (eig_flag) {
    for(i=0; i<3; i++) {
      dir_1[i] = eigvec[i][0];
      dir_2[i] = eigvec[i][1];
      dir_3[i] = eigvec[i][2];
      val[i] = eigval[i];
    }
    
  } 
  else {
    for(i=0; i<3; i++) {
      dir_1[i] = 0.0;
      dir_2[i] = 0.0;
      dir_3[i] = 0.0;
      val[i] = 0.0;
    }
  }
  
  
  VIO_FREE2D(data);
  VIO_FREE2D(weighted_data);
  VIO_FREE2D(covar); 
  VIO_FREE2D(eigvec);
  FREE(eigval);

  return(eig_flag);
}
        

/* 
   The data stored in r[][][] must have a positive semi-definite
   symmetric covariance matrix for eigen values/vectors to be found
   using jacobi (This is tested by positive_3D_semidefinite() in the
   code below).  

   The goal is to have the data fit by a 3D quadratic function,
   with a single (local) minimum.  

   If this is the case, then the routine will return TRUE as well as
   the eigenvectors (dir_1, dir2, dir_3) and the eigenvalues (val[])
   (vectors and values sorted in order of mag(val[]).

*/

VIO_BOOL return_local_eigen_from_hessian(VIO_Real r[3][3][3],
                                               VIO_Real dir_1[3],
                                               VIO_Real dir_2[3],
                                               VIO_Real dir_3[3],
                                               VIO_Real val[3])

{
  int 
    iters,eig_flag,i;
  VIO_Real 
    **covar, **eigvec, *eigval;

  deriv_3D_struct 
    d;                        /* the 1st and second order derivatives */


  VIO_ALLOC2D(covar,3,3);                /* ALLOC the matrices */
  VIO_ALLOC2D(eigvec,3,3);
  ALLOC(eigval,3);

                                /* get the data for the Hessian */
  estimate_3D_derivatives_new(r,&d);

                                /* adjust Hessian matrix, if necessary */

  if ( fabs(d.uu) < SMALL_EPS ) d.uu = 0.0;
  if ( fabs(d.vv) < SMALL_EPS ) d.vv = 0.0;
  if ( fabs(d.ww) < SMALL_EPS ) d.ww = 0.0;
  if ( fabs(d.uv) < SMALL_EPS ) d.uv = 0.0;
  if ( fabs(d.uw) < SMALL_EPS ) d.uw = 0.0;
  if ( fabs(d.vw) < SMALL_EPS ) d.vw = 0.0;


  covar[0][0] = d.uu;
  covar[1][1] = d.vv;
  covar[2][2] = d.ww;
  covar[0][1] = covar[1][0] = d.uv;
  covar[0][2] = covar[2][0] = d.uw;
  covar[1][2] = covar[2][1] = d.vw;

  eig_flag = ( positive_3D_semidefinite(&d) && eigen(covar, 3, eigval, eigvec, &iters));

  if (eig_flag) {
    for(i=0; i<3; i++) {
      val[i] = eigval[i];
      dir_1[i] = eigvec[i][0];
      dir_2[i] = eigvec[i][1];
      dir_3[i] = eigvec[i][2];
    }
  }
  else {
    for(i=0; i<3; i++) {
      val[i] = 0.0;
      dir_1[i] = 0.0;
      dir_2[i] = 0.0;
      dir_3[i] = 0.0;
    }
    dir_3[2] = dir_2[1] = dir_1[0] = 1.0;
  }

  VIO_FREE2D(covar);                /* FREE UP the matrices */
  VIO_FREE2D(eigvec);
  FREE(eigval);


  return(eig_flag);
}
        


