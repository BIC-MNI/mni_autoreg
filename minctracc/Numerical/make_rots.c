#include <def_standard.h>
#include <def_alloc.h>
#include <math.h>

#include "warp_vol.h"
#include "recipes.h"


extern int
  invert_mapping_flag;

extern void
  calc_inv();
extern Boolean
  rotmat_to_ang();


void   make_rots(xmat,data_rot_x,data_rot_y,data_rot_z)
     float
	**xmat,
	data_rot_x,data_rot_y,data_rot_z;
{
   float
      **TRX,
      **TRY,
      **TRZ,
      **T1,
      **T2;
   
   TRX  = matrix(1,4,1,4);
   TRY  = matrix(1,4,1,4);
   TRZ  = matrix(1,4,1,4);
   T1   = matrix(1,4,1,4);
   T2   = matrix(1,4,1,4);

   
   nr_rotxf(TRX,data_rot_x);             /* create the rotate X matrix */
   nr_rotyf(TRY,data_rot_y);             /* create the rotate Y matrix */
   nr_rotzf(TRZ,data_rot_z);             /* create the rotate Z matrix */
   
   nr_multf(TRX,1,4,1,4,  TRY,1,4,1,4,  T1); /* apply rx and ry */
   nr_multf(T1 ,1,4,1,4,  TRZ,1,4,1,4,  xmat); /* apply rz */

   free_matrix(TRX,1,4,1,4);
   free_matrix(TRY,1,4,1,4);
   free_matrix(TRZ,1,4,1,4);
   free_matrix(T1 ,1,4,1,4);
   free_matrix(T2 ,1,4,1,4);

}

public void   
  make_resampling_matrix(xmat,
			 trans_x,trans_y,trans_z,
			 center_x,center_y,center_z,
			 rot_x,rot_y,rot_z,
			 scale_x,scale_y,scale_z)
float
  **xmat,
  trans_x,trans_y,trans_z,
  center_x,center_y,center_z, /* rotation center position */
  rot_x,rot_y,rot_z,   	    /* ang in cntrclkwise radians */
  scale_x,scale_y,scale_z;    /* scale in AP, LR and CC directions */
{

   float
      **TEMP1,
      **TEMP2,
      **R,
      **S,
      **T1,
      **T2;
   int
      i,j;

   TEMP1  = matrix(1,4,1,4);
   TEMP2  = matrix(1,4,1,4);
   R      = matrix(1,4,1,4);
   S      = matrix(1,4,1,4);
   T1     = matrix(1,4,1,4);
   T2     = matrix(1,4,1,4);

  if (!invert_mapping_flag) {
				/* xmat = (-C)*R*S*C*T */

    nr_identf(T1,1,4,1,4);	/* make (-C) */
    T1[4][1] = -center_x;
    T1[4][2] = -center_y;
    T1[4][3] = -center_z;

    make_rots(R,rot_x,rot_y,rot_z); /* make R */

    nr_identf(S,1,4,1,4);	/* make S */
    S[1][1] = scale_x;	
    S[2][2] = scale_y;
    S[3][3] = scale_z;
 
    nr_multf(R,1,4,1,4,  S,1,4,1,4,  T2); /* make   scale*rotation */

    nr_multf(T1,1,4,1,4,  T2,1,4,1,4,  TEMP1); /* apply centering and rotation*scale */
 
    nr_identf(T2,1,4,1,4);
    T2[4][1] = center_x;
    T2[4][2] = center_y;
    T2[4][3] = center_z;

    nr_multf(TEMP1,1,4,1,4, T2,1,4,1,4, TEMP2); /* reposition          */
   
    for (i=1; i<=4; ++i)
      for (j=1; j<=4; ++j)
	xmat[i][j] = TEMP2[i][j];

    xmat[4][1] += trans_x;
    xmat[4][2] += trans_y;
    xmat[4][3] += trans_z;
  }
  else {
				/* xmat = (-T)(-C)(S')(R')(C) where ' is inv*/
    nr_identf(T1,1,4,1,4);

    T1[4][1] = -trans_x-center_x;		
    T1[4][2] = -trans_y-center_y;		
    T1[4][3] = -trans_z-center_z;

    make_rots(R,rot_x,rot_y,rot_z); /* these are inv angles from below */

    nr_identf(S,1,4,1,4);	/* make scaling matrix */
    S[1][1] = 1.0 / scale_x;	
    S[2][2] = 1.0 / scale_y;
    S[3][3] = 1.0 / scale_z;

    nr_multf(S,1,4,1,4,  R,1,4,1,4,  T2); /* make   scale*rotation */

    nr_multf(T1,1,4,1,4,  T2,1,4,1,4,  TEMP1); /* apply centering and rotation*scale */
 
    nr_identf(T2,1,4,1,4);
    T2[4][1] = center_x;
    T2[4][2] = center_y;
    T2[4][3] = center_z;

    nr_multf(TEMP1,1,4,1,4, T2,1,4,1,4, TEMP2); /* reposition          */
   
    for (i=1; i<=4; ++i)
      for (j=1; j<=4; ++j)
	xmat[i][j] = TEMP2[i][j];

  }

   free_matrix(TEMP1,1,4,1,4);
   free_matrix(TEMP2,1,4,1,4);
   free_matrix(R    ,1,4,1,4);
   free_matrix(S    ,1,4,1,4);
   free_matrix(T1   ,1,4,1,4);
   free_matrix(T2   ,1,4,1,4);


}


/* when this routine is called, all calling parameters have valid values
   to map data into target space (forward mapping).
   the object of the routine is to use the information therein, to devise
   the inverse mapping.

   Only the rotation parameters will be modified, since the 
   invert_mapping_flag will control the order of matrix application.
*/

public void
  make_inverted_resampling_matrix(xmat,
				  trans_x,trans_y,trans_z,
				  center_x,center_y,center_z,
				  rot_x,rot_y,rot_z,
				  scale_x,scale_y,scale_z)
float
  **xmat,
  *trans_x,*trans_y,*trans_z,
  *center_x,*center_y,*center_z, /* rotation center position */
  *rot_x,*rot_y,*rot_z,	         /* ang in cntrclkwise degrees! */
  *scale_x,*scale_y,*scale_z;    /* scale in AP, LR and CC directions */
{

  float
    **S,**R,**Ri,
    *ang,
    **TEMP1,
    **TEMP2,
    **T1,
    **T2;
  int
    i,j;

  R     = matrix(1,4,1,4);	/* initial rotation matrix */
  Ri    = matrix(1,4,1,4);	/* inverted matrix */
  TEMP1 = matrix(1,4,1,4);
  ang   = vector(1,3);		/* rot angles for inverted matrix */

  make_rots(R,(*rot_x)*DEG_TO_RAD,(*rot_y)*DEG_TO_RAD,(*rot_z)*DEG_TO_RAD);

  for (i=1; i<=4; ++i)/* copy rot matrix, since calc_inv destroys input */
    for (j=1; j<=4; ++j)
      TEMP1[i][j] = R[i][j];

  calc_inv(TEMP1,Ri,4,4);
    
  if (!rotmat_to_ang(Ri,ang)) {
    print_error("Cannot convert Ri to radians!",__FILE__,__LINE__,0,0,0,0,0);
  }


/*    test by building new new transform matrix, testing it
      against the inv(xmat)                                   */


  S   = matrix(1,4,1,4);
  TEMP2  = matrix(1,4,1,4);
  T1     = matrix(1,4,1,4);
  T2     = matrix(1,4,1,4);

  make_resampling_matrix(TEMP2,
			 *trans_x,*trans_y,*trans_z,
			 *center_x,*center_y,*center_z,
			 ang[1],ang[2],ang[3],
			 *scale_x,*scale_y,*scale_z);

  for (i=1; i<=4; ++i)		/* reset all values of transformation matrix */
    for (j=1; j<=4; ++j)
      T2[i][j] = xmat[i][j];

  calc_inv(T2,TEMP1,4,4);

  nr_multf(xmat,1,4,1,4,  TEMP2,1,4,1,4,  T1); /* T1 should = I */

  for (i=1; i<=4; ++i) {
    for (j=1; j<=4; ++j)
      printf ("%8.4f ",TEMP2[i][j]);
    printf ("|");
    for (j=1; j<=4; ++j)
      printf ("%8.4f ",xmat[i][j]);
    printf ("|");
    for (j=1; j<=4; ++j)
      printf ("%8.4f ",T1[i][j]);
    printf ("\n");
  }

  for (i=1; i<=4; ++i) {
    for (j=1; j<=4; ++j)
      printf ("%8.4f ",TEMP1[i][j]);
    printf ("\n");
  }
  
  for (i=1; i<=4; ++i)		/* reset all values of transformation matrix */
    for (j=1; j<=4; ++j)
      xmat[i][j] = TEMP2[i][j];
  
  *rot_x = ang[1]*RAD_TO_DEG;
  *rot_y = ang[2]*RAD_TO_DEG;
  *rot_z = ang[3]*RAD_TO_DEG;

  free_matrix(TEMP2,1,4,1,4);
  free_matrix(S    ,1,4,1,4);
  free_matrix(T1   ,1,4,1,4);
  free_matrix(T2   ,1,4,1,4);

  free_matrix(R    ,1,4,1,4);
  free_matrix(Ri   ,1,4,1,4);
  free_vector(ang,  1,3);
  free_matrix(TEMP1,1,4,1,4);

}


