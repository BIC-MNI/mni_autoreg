/* ----------------------------- MNI Header -----------------------------------
@NAME       : rotmat_to_ang.c
@INPUT      : rot      - rotation matrix (3x3 in num recipes form) calculated
                         by the calling program.
@OUTPUT     : ang      - vector giving rx,ry and rz rotation angles (in 
                         radians). This vector must be defined by the 
			 calling routine.
@RETURNS    : TRUE if ok, FALSE if error.
@DESCRIPTION: this routine extracts the rotation angles from the rotation
              matrix.  The rotation matrix is assumed to be a 3x3 matrix in
	      numerical recipes form.  It is locally copied into a 
	      4x4 homogeneous matrix for manipulation.

	      we assume that the matrix rotation center is (0,0,0).
	      we assume that the application of this matrix to a vector
	         is done with rot_mat*vec = premultiplication by matrix

	      the resulting angles rx=ang[1],ry=ang[2],rz=ang[3], follow
	      the following assumptions:

	      -all rotations are assumed to be in the range -pi/2..pi/2
	       routine returns FALSE is this is found not to hold
	      -rotations are applied 1 - rx, 2 - ry and 3 - rz
	      -applying these rotations to an identity matrix will
               result in a matrix equal to `rot' (the input matrix)
	      -positive rotations are counter-clockwise when looking
	       down the axis, from the positive end towards the origin.
	      -I assume a coordinate system:
	                  ^ Y
	                  |
	                  |
	                  |
			  |
	                  +---------> X
	                 /
	                /
	               / Z  (towards the viewer).

	      -The problem is posed as:  
	         given a rotation matrix ROT, determine the rotations
		 rx,ry,rz applied in order to give ROT
	       solution:
	         assume the rot matrix is equivalent to a normalized
		 orthogonal local coord sys.  i.e.  row 1 of ROT is
		 the local x direction, row 2 is the local y and row 3
		 is the local z.
	     
		 (note local is lower case, world is UPPER)

	         1- find RZ that brings local x into the XZ plane
		 2- find RY that brings local x*RZ onto X
		 3- find RX that brings local y*RZ*RY onto Y

		 the required rotations are -RX,-RY and -RZ!

@GLOBALS    : 
@CALLS      : mfmult(), rotx(),roty(),rotz(),matrix(),vector()
@CREATED    : Feb 9, 1992 lc
@MODIFIED   : 
Tue Jun  8 08:44:59 EST 1993 LC
   changes all vec*matrix to matrix*vec.  Std is premultiplication by matrix!

---------------------------------------------------------------------------- */
#include <def_mni.h>
#include <recipes.h>

extern char *prog_name;

extern void 
  print_error(char *s, char * d1, int d2, int d3, int d4, int d5, int d6, int d7);


#define EPS  0.00000000001	/* epsilon, should be calculated */

#include "matrix_basics.h"

public Boolean rotmat_to_ang(float **rot, float *ang)
{

   float 
      rx,ry,rz,
      **t,**s,
      **R,
      **Rx,
      **Ry,
      **Rz,
      len,
      i,j,k;

   int
      m,n;

   t  = matrix(1,4,1,1);	/* make two column vectors */
   s  = matrix(1,4,1,1);

   R  = matrix(1,4,1,4);	/* working space matrices */
   Rx = matrix(1,4,1,4);
   Ry = matrix(1,4,1,4);
   Rz = matrix(1,4,1,4);

   nr_identf(R,1,4,1,4);	/* init R homogeneous matrix */

   for (m=1; m<=3; ++m)		/* copy rot matrix into R */
      for (n=1; n<=3; ++n)
	 R[m][n] = rot[m][n];
   
/* ---------------------------------------------------------------
   step one,  find the RZ rotation reqd to bring 
              the local x into the world XZ plane
*/

   for (m=1; m<=3; ++m)		/* extract local x vector, ie the first column */
      t[m][1] = R[m][1];
   t[4][1] = 1.0;

   i = t[1][1];			/* make local vector componants */
   j = t[2][1]; 
   k = t[3][1];

   if (i<EPS) {			/* if i is not already in the positive X range, */
      print_error("step one: rz not in the range -pi/2..pi/2",__FILE__, __LINE__,0,0,0,0,0);
      return(FALSE);
   }

   len = sqrt(i*i + j*j);	/* length of vect x on XY plane */
   if (ABS(len)<EPS) {
      print_error("step one: length of vect x null.",__FILE__, __LINE__,0,0,0,0,0);
      return(FALSE);
   }

   if (ABS(i)>ABS(j)) {
      rz = ABS(fasin(j/len));
   }
   else {
      rz = ABS(facos(i/len));
   }

   if (j>0)			/* what is the counter clockwise angle */
      rz = -rz;                 /* necessary to bring vect x ont XY plane? */
      
  
/*---------------------------------------------------------------
   step two:  find the RY rotation reqd to align 
              the local x on the world X axis 

  (since i was positive above, RY should already by in range -pi/2..pi/2 
  but we'll check it  anyway)                                             */

   for (m=1; m<=3; ++m)		/* extract local x vector */
      t[m][1] = R[m][1];
   t[4][1] = 1.0;

   nr_rotzf(Rz,rz);             /* create the rotate Z matrix */
 
   nr_multf(Rz,1,4,1,4,  t,1,4,1,1,   s);   /* apply RZ, to get x in XZ plane */

   i = s[1][1];			/* make local vector componants */
   j = s[2][1]; 
   k = s[3][1];

   if (i<EPS) {
      print_error("step two: ry not in the range -pi/2..pi/2",__FILE__, __LINE__,0,0,0,0,0);
      return(FALSE);
   }

   len = sqrt(i*i + k*k);		/* length of vect x in XZ plane, after RZ */

   if (ABS(len)<EPS) {
      print_error("step two: length of vect z null.",__FILE__, __LINE__,0,0,0,0,0);
      return(FALSE);
   }

   if (ABS(i)>ABS(k)) {
      ry = ABS(fasin(k/len));
   }
   else {
      ry = ABS(facos(i/len));
   }

   /*    what is the counter clockwise angle necessary to bring  */
   /*    vect x onto X? */
   if (k < 0) { 
      ry = -ry;
   }

   /*--------------------------------------------------------------- */
   /*   step three,rotate around RX to */
   /*              align the local y with Y and z with Z */

   for (m=1; m<=3; ++m)		/* extract local z vector */
      t[m][1] = R[m][3];
   t[4][1] = 1.0;

   nr_rotyf(Ry,ry);             /* create the rotate Y matrix */

				/* t =  roty(ry*180/pi) *(rotz(rz*180/pi) *r(3,:)); */
   nr_multf(Rz,1,4,1,4,  t,1,4,1,1,  s); /* apply RZ, to get x in XZ plane */
   nr_multf(Ry,1,4,1,4,  s,1,4,1,1,  t); /* apply RY, to get x onto X      */

   i = t[1][1];			/* make local vector componants */
   j = t[2][1]; 
   k = t[3][1];

   len = sqrt(j*j + k*k);	/* length of vect x in Y,Z plane */

   if (ABS(len)<EPS) {
      print_error("step three: length of vect z null.",__FILE__, __LINE__,0,0,0,0,0);
      return(FALSE);
   }

   if (ABS(k)>ABS(j)) {
      rx = ABS(fasin(j/len));
   }
   else {
      rx = ABS(facos(k/len));
   }

   if (j< 0) { 
      rx = -rx;
   }
	
   rx = -rx;  /* these are the required rotations */
   ry = -ry;
   rz = -rz;

   ang[1] = rx;
   ang[2] = ry;
   ang[3] = rz;

   free_matrix(t,1,4,1,1);
   free_matrix(s,1,4,1,1);
   free_matrix(R,1,4,1,4);
   free_matrix(Rx,1,4,1,4);
   free_matrix(Ry,1,4,1,4);
   free_matrix(Rz,1,4,1,4);

   return(TRUE);
}
