#include <volume_io.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include "cov_to_praxes.h"

time_t time(time_t *tloc);

char *prog_name;

VIO_BOOL eigen(double **inputMat, int ndim, 
                     double *eigen_val, double **eigen_vec, 
                     int    *iters);

main(int argc, char *argv[]) {

  VIO_Real
    **mat, **mat_saved,**mat_input, *eig_val, **eig_vec, **diag_val,
    temp;
  int
    counter,flag,iters,i,j,k,dim;
  float
     **axes;
  union
    {
      long   l;
      char   c[4];
    } seedval;
  time_t t;
  char tmp;
  

  prog_name = argv[0];

  dim = 3;
  
  ALLOC2D(mat_input,dim+1,dim+1);
  ALLOC2D(mat_saved,dim+1,dim+1);
  ALLOC2D(mat,      dim+1,dim+1);
  ALLOC2D(eig_vec,  dim+1,dim+1);
  ALLOC(eig_val,    dim+1);
  ALLOC2D(axes,     dim+1,dim+1);
  ALLOC2D(diag_val, dim+1,dim+1);
  
  t = 2*time(NULL);
  seedval.l = t; 
  
  tmp = seedval.c[0]; seedval.c[0] = seedval.c[3]; seedval.c[3] = tmp; 
  tmp = seedval.c[1]; seedval.c[1] = seedval.c[2]; seedval.c[2] = tmp;
  
  srand48(seedval.l);
  
  
/*
   test cov_to_praxes

  mat_input[1][1] = drand48();
  mat_input[2][2] = drand48();
  mat_input[3][3] = drand48();
  mat_input[1][2] = mat_input[2][1] = drand48();
  mat_input[1][3] = mat_input[3][1] = drand48();
  mat_input[2][3] = mat_input[3][2] = drand48();
  

  for(i=1; i<=3; i++) {
    for(j=1; j<=3; j++)
      print ("%10.6f ",mat_input[i][j]);
    print ("\n");
  }

  cov_to_praxes( dim, mat_input, axes);


  print("-> %f %f %f\n", axes[1][1], axes[1][2], axes[1][3]);
  print("-> %f %f %f\n", axes[2][1], axes[2][2], axes[2][3]);
  print("-> %f %f %f\n", axes[3][1], axes[3][2], axes[3][3]);

*/

/*
   test eigen 
*/

  counter = 0;
  
  do {
    
    mat_input[0][0] = drand48();
    mat_input[1][1] = drand48();
    mat_input[2][2] = drand48();
    mat_input[0][1] = mat_input[1][0] = drand48();
    mat_input[0][2] = mat_input[2][0] = drand48();
    mat_input[1][2] = mat_input[2][1] = drand48();
    
    
    for(i=0; i<=2; i++) 
      for(j=0; j<=2; j++) 
        mat_saved[i][j]=mat_input[i][j];

    flag = (mat_input[0][0]>0.0 &&
            (mat_input[0][0]*mat_input[1][1] - mat_input[0][1]*mat_input[1][0])>0 &&
            ((mat_input[0][0] * (mat_input[1][1]*mat_input[2][2] - mat_input[1][2]*mat_input[2][1])) -
             (mat_input[0][1] * (mat_input[1][0]*mat_input[2][2] - mat_input[1][2]*mat_input[2][0])) +
             (mat_input[0][2] * (mat_input[1][0]*mat_input[2][1] - mat_input[1][1]*mat_input[2][0]))
             ) > 0.0
            );
    
    (void)eigen(mat_saved, 3, 
                eig_val, eig_vec, &iters); 
    
    if (flag) {
      counter++;
      print ("%d \n",counter);
    }

    if ( !flag && eig_val[0]>0.0 &&eig_val[1]>0.0 &&eig_val[2]>0.0 ) {
      
      for(i=0; i<=2; i++) {
        for(j=0; j<=2; j++) {
          print ("%10.6f ",mat_input[i][j]);
          mat[i][j] = mat_input[i][j];
        }
        print ("\n");
      }
      
      if (flag) 
        print ("positive def\n");
      else
        print ("not positive def\n");

      
      print ("iters: %d\n",iters);
      for(i=0; i<=2; i++) {
        for(j=0; j<=2; j++) {
          if (i==j)
            diag_val[i][j] = eig_val[i];
          else
            diag_val[i][j] = 0.0;
        }
      }
      
      
      
      for(i=0; i<=2; i++) {
        for(j=0; j<=2; j++) {
          print ("%10.6f ",eig_vec[i][j]);
        }
        print ("   ");
        for(j=0; j<=2; j++) {
          print ("%10.6f ",diag_val[i][j]);
        }
        print ("\n");
      }
      
      
      for(i=0; i<=2; i++) {
        for(j=0; j<=2; j++) {
          mat_input[i][j] = 0.0;
          for(k=0; k<=2; k++) {
            mat_input[i][j] += mat[i][k] * eig_vec[k][j];
          }
        }
      }
      
      print ("input_mat * eig_vec:\n");
      for(i=0; i<=2; i++) {
        for(j=0; j<=2; j++) {
          print ("%10.6f ",mat_input[i][j]);
        }
        print ("\n");
      }
    
      for(i=0; i<=2; i++) {
        for(j=0; j<=2; j++) {
          mat_input[i][j] = 0.0;
          for(k=0; k<=2; k++) {
            mat_input[i][j] += eig_vec[i][k] * diag_val[k][j];
          }
        }
      }
      
      print ("eig_vec * diag(eig_val):\n");
      
      for(i=0; i<=2; i++) {
        for(j=0; j<=2; j++) {
          print ("%10.6f ",mat_input[i][j]);
        }
        print ("\n");
      }
    }
  } while (TRUE);

  FREE2D(mat_input);
  FREE2D(mat);
  FREE2D(eig_vec);
  FREE(eig_val);
  FREE2D(axes);
  FREE2D(diag_val);




}

