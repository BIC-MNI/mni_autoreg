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
	      
	      the resulting angles rx=ang[1],ry=ang[2],rz=ang[3], follow
	      the following assumptions:

	      -all rotations are assumed to be in the range -pi/2..pi/2
	       routine returns FALSE is this is found not to hold
	      -rotations are applied 1 - rx, 2 - ry and 3 - rz
	      -applying these rotations to an identity matrix will
               result in a matrix equal to `rot' (the input matrix)
	      -positive rotations are counter-clockwise when looking
	       down the axis, from the positive end towards the origin.
	      -i assume a coordinate system:
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
---------------------------------------------------------------------------- */
#include <def_mni.h>
#include <recipes.h>

extern char *prog_name;

extern void
   print_error();

#define EPS  0.00000000001	/* epsilon, should be calculated */


extern int
   nr_identf(),
   nr_multf(),
   nr_rotxf(),
   nr_rotyf(),
   nr_rotzf();

public Boolean rotmat_to_ang(rot, ang)
float 
   **rot,
   *ang;
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

   t  = matrix(1,1,1,4);
   s  = matrix(1,1,1,4);
   R  = matrix(1,4,1,4);
   Rx = matrix(1,4,1,4);
   Ry = matrix(1,4,1,4);
   Rz = matrix(1,4,1,4);

   nr_identf(R,1,4,1,4);		/* init R homogeneousmatrix */

   for (m=1; m<=3; ++m)		/* copy rot matrix into R */
      for (n=1; n<=3; ++n)
	 R[m][n] = rot[m][n];
   
/* ---------------------------------------------------------------
   step one,  find the RZ rotation reqd to bring 
              the local x into the world XZ plane
*/

   for (m=1; m<=3; ++m)		/* extract local x vector */
      t[1][m] = R[1][m];
   t[1][4] = 1.0;

   i = t[1][1];			/* make local vector componants */
   j = t[1][2]; 
   k = t[1][3];

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
      t[1][m] = R[1][m];
   t[1][4] = 1.0;

   nr_rotzf(Rz,rz);             /* create the rotate Z matrix */
 
   nr_multf(t,1,1,1,4,   Rz,1,4,1,4,   s);   /* apply RZ, to get x in XZ plane */

   i = s[1][1];			/* make local vector componants */
   j = s[1][2]; 
   k = s[1][3];

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
      t[1][m] = R[3][m];
   t[1][4] = 1.0;

   nr_rotyf(Ry,ry);             /* create the rotate Y matrix */

				/* t = (r(3,:)*rotz(rz*180/pi))*roty(ry*180/pi);   */
   nr_multf(t,1,1,1,4,  Rz,1,4,1,4,  s); /* apply RZ, to get x in XZ plane */
   nr_multf(s,1,1,1,4,  Ry,1,4,1,4,  t); /* apply RY, to get x onto X      */

   i = t[1][1];			/* make local vector componants */
   j = t[1][2]; 
   k = t[1][3];

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

   free_matrix(t,1,1,1,4);
   free_matrix(s,1,1,1,4);
   free_matrix(R,1,4,1,4);
   free_matrix(Rx,1,4,1,4);
   free_matrix(Ry,1,4,1,4);
   free_matrix(Rz,1,4,1,4);

   return(TRUE);
}
