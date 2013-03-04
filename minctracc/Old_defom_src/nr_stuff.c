#include <math.h>


nr_identd( A, m1, m2, n1, n2 )
double **A;
int m1,m2,n1,n2;
{

   int i,j;

   for (i=m1; i<=m2; ++i)
      for (j=n1; j<=n2; ++j) {
         if (i==j) 
            A[i][j] = 1.0;
         else
            A[i][j] = 0.0;
      }
   
}

nr_identf( A, m1, m2, n1, n2 )
float **A;
int m1,m2,n1,n2;
{

   int i,j;

   for (i=m1; i<=m2; ++i)
      for (j=n1; j<=n2; ++j) {
         if (i==j) 
            A[i][j] = 1.0;
         else
            A[i][j] = 0.0;
      }
   
}

nr_copyd( A, m1,m2, n1,n2, B )
/* space must already be allocated for B */
double **A, **B;
int m1,m2,n1,n2;
{
   int i,j;

   for (i=m1; i<=m2; ++i)
      for (j=n1; j<=n2; ++j)
         B[i][j] = A[i][j];
}

nr_copyf( A, m1,m2, n1,n2, B )
/* space must already be allocated for B */
float **A, **B;
int m1,m2,n1,n2;
{
   int i,j;

   for (i=m1; i<=m2; ++i)
      for (j=n1; j<=n2; ++j)
         B[i][j] = A[i][j];
}



nr_identd( A, m1, m2, n1, n2 )
double **A;
int m1,m2,n1,n2;
{

   int i,j;

   for (i=m1; i<=m2; ++i)
      for (j=n1; j<=n2; ++j) {
         if (i==j) 
            A[i][j] = 1.0;
         else
            A[i][j] = 0.0;
      }
   
}

nr_identf( A, m1, m2, n1, n2 )
float **A;
int m1,m2,n1,n2;
{

   int i,j;

   for (i=m1; i<=m2; ++i)
      for (j=n1; j<=n2; ++j) {
         if (i==j) 
            A[i][j] = 1.0;
         else
            A[i][j] = 0.0;
      }
   
}


/* --------------------------------------------------------
   rx =[1   0      0      0 
        0  cos(a)  sin(a) 0
        0 -sin(a)  cos(a) 0
        0   0      0      1];
*/


nr_rotxd(M,a)
double 
   **M;
double
   a;
{
   nr_identd(M,1,4,1,4);

   M[2][2] = cos(a);    M[2][3] = sin(a);
   M[3][2] = -sin(a);   M[3][3] = cos(a);
}


nr_rotxf(M,a)
float
   **M;
float
   a;
{
   nr_identf(M,1,4,1,4);

   M[2][2] = fcos(a);    M[2][3] = fsin(a);
   M[3][2] = -fsin(a);   M[3][3] = fcos(a);
}


/* --------------------------------------------------------
ry = [  cos(a)   0 -sin(a)  0 
        0       1   0       0
        sin(a)  0  cos(a)   0
        0   0      0      1];
*/

nr_rotyd(M,a)
double 
   **M;
double
   a;
{

   nr_identd(M,1,4,1,4);

   M[1][1] = cos(a);    M[1][3] = -sin(a);
   M[3][1] = sin(a);   M[3][3] = cos(a);
}

nr_rotyf(M,a)
float
   **M;
float
   a;
{

   nr_identf(M,1,4,1,4);

   M[1][1] = fcos(a);    M[1][3] = -fsin(a);
   M[3][1] = fsin(a);   M[3][3] = fcos(a);
}


/* --------------------------------------------------------
rz = [cos(a)  sin(a) 0  0
      -sin(a) cos(a) 0  0
        0     0      1  0
        0     0      0  1];
*/

nr_rotzd(M,a)
double 
   **M;
double
   a;
{

   nr_identd(M,1,4,1,4);

   M[1][1] = cos(a);   M[1][2] = sin(a);
   M[2][1] = -sin(a);  M[2][2] = cos(a);
}

nr_rotzf(M,a)
float
   **M;
float
   a;
{

   nr_identf(M,1,4,1,4);

   M[1][1] = fcos(a);   M[1][2] = fsin(a);
   M[2][1] = -fsin(a);  M[2][2] = fcos(a);
}
