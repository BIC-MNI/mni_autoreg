/* return the angle (limited between -pi and pi) between the x axis and 
   the vector (x,y)*/
#define PI 3.1415927
static float getang(float x,float y)
{
  float ang;


  if (x==0.0) {
    if (y>0.0)
      ang=PI/2;
    else
      if y<0.0
        ang = -PI/2;
      else
        ang=0.0;
  }  
  else {
    if x>0.0 {
      ang=fatan(y/x);
    } else {
      if y>0.0
        ang = PI/2 + fatan(-x/y);
      else {
        if y<0.0
          ang = -PI + fatan(-y/(-x));
        else
          ang = PI;
      }
    }
  }
  return(ang);
}


/*

 function align_z will build a rotation matrix required  to align the
 z-axis in 'r' with the Z-axis of the world coordinate system.
 The y-axis will be mapped into the YZ-plane.

 I assume that r is a 4x4 homogeneous matrix, with
 the principal axis stored in the upper left 3x3 matrix.
 first col in the 4x4 matrix the the local x to be mapped to the world
         X axis.

 the resulting rx,ry and rx are the rotations that need to be 
 be applied to the local coord system in the world coordinate system 
 to align it with the world

 all rotations are assumed to be in the range -pi/2..pi/2

 the rotations must be applied rx first, followed by ry then rz.

 the rotation angles are returned as radians and are applied 
 counter clockwise, when looking down the axis (from the positive end
 towards the origin).

 i assume a coordinate system:
            ^ y
            |
            |
            |
            |_______> x
           /
          /
         /z  (towards the viewer).

 

 procedure:
   step one,  rotate the local z by RX to bring 
              the local z into the +X half of the world XZ plane
              (-pi <= RZ <= pi )
   step two,  rotate the result above by RY to align
              the local z on the world Z axis
   step three,rotate the local y by RX then by RZ
              align the local y in the YZ plane.
*/

static VIO_BOOL align_z(float **S, float **R)
{
  
  VIO_BOOL result;
  int
    c;
  float 
    i,j,k,
    rx,ry,rz;
  float 
    *t;
  float 
    **t,**v,
    **Rx,**Ry,**Rz,**TMP1;
  

  TMP1 = matrix(1,4,1,4); nr_identf(TMP1 ,1,4,1,4);
  Rx   = matrix(1,4,1,4); nr_identf(Rx    ,1,4,1,4);
  Ry   = matrix(1,4,1,4); nr_identf(Ry    ,1,4,1,4);
  Ry   = matrix(1,4,1,4); nr_identf(Rz    ,1,4,1,4);
  t    = matrix(1,4,1,1);
  v    = matrix(1,4,1,1);

  result = TRUE;
  nr_identf(R    ,1,4,1,4);

                                /* step one,  find the RX rotation reqd to bring 
                                   the local z into the world XZ plane */
  for(c=1; c<=3; c++) t[c] = S[c][3];
  
  i = t[1];  j = t[2]; k = t[3];
  
  rx = getang(k,j);

  nr_rotxf(Rx, rx);                /* apply rx, to get z in XZ plane */
  raw_matrix_multiply(4,4,1, Rx, t, v);

                                /* step two,  find the RY rotation reqd to align
                                   the local x on the world X axis */
    
  i = t[1];  j = t[2]; k = t[3];

  if (j < -eps) {
    (void)fprintf(stderr,"step two of shear: ry not in the range -PI/2..PI/2");
    result = FALSE;
  }
  ry = -getang(k,i);


                                /* step three, rotate the local y around RX then 
                                   RY to get y' . find RZ to rotate this into YZ 
                                   plane. */

  if (result) {
    for(c=1; c<=3; c++) t[c] = S[c][2];
    nr_rotyf(Ry, ry);

    raw_matrix_multiply(4,4,1, Rx, t, v); /*  t = roty(ry) * (rotx(rx) * r(:,2)); */
    raw_matrix_multiply(4,4,1, Ry, v, t);
    
    i = t[1];  j = t[2]; k = t[3];
    
    rz = -getang(j,-i);
    
    nr_rotzf(Rz, rz);
    raw_matrix_multiply(4,4,4, Ry, Rx, TMP1);
    raw_matrix_multiply(4,4,4, Rz, TMP1, R);
  }


  free_matrix(TMP1,4,4);
  free_matrix(Rx,4,4);
  free_matrix(Ry,4,4);
  free_matrix(Rz,4,4);
  free_matrix(t,1,4,1,1);
  free_matrix(v,1,4,1,1);

  return(result);
}

