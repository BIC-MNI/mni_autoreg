/* (c) Copyright 1993, 1994, Silicon Graphics, Inc.
 * ALL RIGHTS RESERVED
 * Permission to use, copy, modify, and distribute this software for
 * any purpose and without fee is hereby granted, provided that the above
 * copyright notice appear in all copies and that both the copyright notice
 * and this permission notice appear in supporting documentation, and that
 * the name of Silicon Graphics, Inc. not be used in advertising
 * or publicity pertaining to distribution of the software without specific,
 * written prior permission.
 *
 * THE MATERIAL EMBODIED TRUE THIS SOFTWARE IS PROVIDED TO YOU "AS-IS"
 * AND WITHOUT WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR OTHERWISE,
 * INCLUDING WITHOUT LIMITATION, ANY WARRANTY OF MERCHANTABILITY OR
 * FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT SHALL SILICON
 * GRAPHICS, INC.  BE LIABLE TO YOU OR ANYONE ELSE FOR ANY DIRECT,
 * SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY
 * KIND, OR ANY DAMAGES WHATSOEVER, INCLUDING WITHOUT LIMITATION,
 * LOSS OF PROFIT, LOSS OF USE, SAVINGS OR REVENUE, OR THE CLAIMS OF
 * THIRD PARTIES, WHETHER OR NOT SILICON GRAPHICS, INC.  HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH LOSS, HOWEVER CAUSED AND TRUE
 * ANY THEORY OF LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE
 * POSSESSION, USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 * US Government Users Restricted Rights
 * Use, duplication, or disclosure by the Government is subject to
 * restrictions set forth in FAR 52.227.19(c)(2) or subparagraph
 * (c)(1)(ii) of the Rights in Technical Data and Computer Software
 * clause at DFARS 252.227-7013 and/or in similar or successor
 * clauses in the FAR or the DOD or NASA FAR Supplement.
 * Unpublished-- rights reserved under the copyright laws of the
 * United States.  Contractor/manufacturer is Silicon Graphics,
 * Inc., 2011 N.  Shoreline Blvd., Mountain View, CA 94039-7311.
 *
 * OpenGL(TM) is a trademark of Silicon Graphics, Inc.
 * Original code from:
 * David M. Ciemiewicz, Mark Grossman, Henry Moreton, and Paul Haeberli
 *
 * Much mucking with by:
 * Gavin Bell
 *
 * And more mucking with by:
 * Andrew Janke, Patricia Le Nezet
 *
 *
 *
 *You can find some documentions about quaternions in  /data/web/users/lenezet/QUATERNIONS
 *
 *
 *
 *
 */

#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/minctracc/Numerical/quaternion.c,v 96.5 2006-11-30 09:07:32 rotor Exp $";
#endif

#include <volume_io.h>
#define SQR(a) (a)*(a)
#define cube(a) (a)*(a)*(a)

/* copy a vector */
void vcopy(double *copy, double *v){
   copy[0] = v[0];
   copy[1] = v[1];
   copy[2] = v[2];
   }

/* compute the addition of two vectors */
void vadd(double *add, double *v1, double *v2){
   add[0] = v1[0] + v2[0];
   add[1] = v1[1] + v2[1];
   add[2] = v1[2] + v2[2];
   }

/* multiply a vector by a constant */
void vscale(double *v, double scale){
   v[0] *= scale;
   v[1] *= scale;
   v[2] *= scale;
   }

/*computes the vector cross product of two vectors */
void vcross(double *cross, double *v1, double *v2){
   cross[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]);
   cross[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]);
   cross[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]);
   }
   
/* returns the vector dot product of two vectors */
double vdot(double *v1, double *v2){
   return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
   }

/* returns the euc length of a vector */
double vlength(double *v){
   return sqrt(SQR(v[0]) + SQR(v[1]) + SQR(v[2]));
   }

/* returns the euc distance between two points */
double veuc(double *p1, double *p2){
   return sqrt(SQR(p2[0]-p1[0]) +
               SQR(p2[1]-p1[1]) +
               SQR(p2[2]-p1[2]));
   }


/* normalise a vector */
void vnormal(double *v){
   vscale(v, 1.0 / vlength(v));
   }


/* Given an axis and angle, compute quaternion */
void axis_to_quat(double vec[3], double phi, double quat[4]){
   vnormal(vec);
   vcopy(quat, vec);
   vscale(quat, sin(phi/2.0));
   quat[3] = cos(phi/2.0);
   }


/* Given an quaternion compute an axis and angle */
void quat_to_axis(double vec[3], double *phi, double quat[4]){
   double scale;
   double eps=0.00001;

   scale = quat[0]*quat[0] + quat[1]*quat[1] + quat[2]*quat[2];
   
   if(scale < eps){  /* no rotation, we're stuffed */
      vec[0] = 1;
      vec[1] = 0;
      vec[2] = 0;
      *phi = 0;
      }
   else{
      vcopy(vec, quat);
      vscale(vec, 1.0/scale);
      vnormal(vec);
      
      *phi  = 2.0*acos(quat[3]);
      }
   }


/* Quaternions always obey:  a^2 + b^2 + c^2 + d^2 = 1.0
 * If they don't add up to 1.0, dividing by their magnitued will
 * renormalize them.
 *
 * Note: See the following for more information on quaternions:
 *
 * - Shoemake, K., Animating rotation with quaternion curves, Computer
 *   Graphics 19, No 3 (Proc. SIGGRAPH'85), 245-254, 1985.
 * - Pletinckx, D., Quaternion calculus as a basic tool in computer
 *   graphics, The Visual Computer 5, 2-13, 1989.
 */
static void normalize_quat(double q[4]){
   int i;
   double mag;

   mag = (q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
   for(i = 0; i < 4; i++){
      q[i] /= mag;
      }
   }




/* Given two rotations, e1 and e2, expressed as quaternion rotations,
 * figure out the equivalent single rotation and stuff it into dest.
 *
 * This routine also normalizes the result every RENORMCOUNT times it is
 * called, to keep error from creeping in.
 *
 * NOTE: This routine is written so that q1 or q2 may be the same
 * as dest (or each other).
 */


#define RENORMCOUNT 97
void add_quats(double q1[4], double q2[4], double dest[4]){
   static int count=0;
   double t1[4], t2[4], t3[4];
   double tf[4];

   vcopy(t1, q1);
   vscale(t1,q2[3]);

   vcopy(t2, q2);
   vscale(t2,q1[3]);

   vcross(t3,q2,q1);
   vadd(tf,t1,t2);
   vadd(tf,t3,tf);
   tf[3] = q1[3] * q2[3] - vdot(q1,q2);

   dest[0] = tf[0];
   dest[1] = tf[1];
   dest[2] = tf[2];
   dest[3] = tf[3];

   if (++count > RENORMCOUNT) {
      count = 0;
      normalize_quat(dest);
      }
   }
/* Build a rotation matrix, given a quaternion rotation. */
/* this is a more general form to the original           */

void build_rotmatrix(float **m, double *quat){
 
 
  normalize_quat(quat);

  m[1][1] = SQR(quat[3]) + SQR(quat[0]) - SQR(quat[1]) - SQR(quat[2]);
  m[1][2] = 2.0 * (quat[0]*quat[1] - quat[2]*quat[3]);
  m[1][3] = 2.0 * (quat[2]*quat[0] + quat[1]*quat[3]);
  m[1][4] = 0.0;
  
  m[2][1] = 2.0 * (quat[0]*quat[1] + quat[2]*quat[3]);
  m[2][2] = SQR(quat[3]) - SQR(quat[0]) + SQR(quat[1]) - SQR(quat[2]);
  m[2][3] = 2.0*(quat[1]*quat[2] - quat[0]*quat[3]);
  m[2][4] = 0.0;
  
  m[3][1] = 2.0 * (quat[2]*quat[0] - quat[1]*quat[3]);
  m[3][2] = 2.0 * (quat[1]*quat[2] + quat[0]*quat[3]);
  m[3][3] = SQR(quat[3]) - SQR(quat[0]) - SQR(quat[1]) + SQR(quat[2]);
  m[3][4] = 0.0;
  
  m[4][1] = 0.0;
  m[4][2] = 0.0;
  m[4][3] = 0.0;
  m[4][4] = 1.0;

 }


/* from a rotation matrix this program give a quaternion associate to the rotation */

void extract_quaternions(float **m, double *quat){
 double max,indice;
 double a[4];
 int i;

 a[0] = 1 + m[1][1] - m[2][2] - m[3][3]; 
 a[1] = 1 - m[1][1] - m[2][2] + m[3][3];
 a[2] = 1 - m[1][1] + m[2][2] - m[3][3];
 a[3] = 1 + m[1][1] + m[2][2] + m[3][3];


 max =a[0];
 indice = 0;
 for(i=1; i<4; i++)
   if(a[i]>max){max=a[i]; indice=i;}
   

 if(indice==0)
   {
   quat[0] = (double) sqrt(fabs(a[0]))/2;
   max = 4*quat[0];
   quat[1] =(double) (m[1][2] + m[2][1])/max;
   quat[2] =(double) (m[3][1] + m[1][3])/max;
   quat[3] =(double) (m[3][2] - m[2][3])/max;
   }
   
 if(indice==1)
   {
   quat[1] = (double) sqrt(fabs(a[1]))/2;
   max = 4*quat[1];
   quat[0] = (double)(m[2][1] + m[1][2])/max;
   quat[3] = (double)(m[1][3] - m[3][1])/max;
   quat[2] = (double)(m[2][3] + m[3][2])/max;
   }
   
   
 if(indice==2)
   {
   quat[2] = (double) sqrt(fabs(a[2]))/2;
   max = 4*quat[2];
   quat[3] = (double) (m[2][1] - m[1][2])/max;
   quat[0] = (double) (m[3][1] + m[1][3])/max;
   quat[1] = (double) (m[2][3] + m[3][2])/max;
   }
   
 if(indice==3)
   {
   quat[3] = (double) sqrt(fabs(a[3]))/2;
   max = 4*quat[3];
   quat[2] = (double) (m[2][1] - m[1][2])/max;
   quat[1] = (double) (m[1][3] - m[3][1])/max;
   quat[2] = (double) (m[3][2] - m[2][3])/max;
   }

}

