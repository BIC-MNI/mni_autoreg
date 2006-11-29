/* ----------------------------- MNI Header -----------------------------------
@NAME       : test_eigen.c
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : status
@DESCRIPTION: Program to test calculation of eigen vectors and quadratic fitting.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Mon Sep 25 08:45:43 MET 1995
@MODIFIED   : $Log: test_eigen.c,v $
@MODIFIED   : Revision 1.5  2006-11-29 09:09:31  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 1.4  2005/07/20 20:45:47  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 1.3  2004/02/12 05:54:06  rotor
@MODIFIED   :  * removed public/private defs
@MODIFIED   :
@MODIFIED   : Revision 1.2  2002/03/26 14:15:33  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 1.1  1999/10/25 19:52:11  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :

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
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Extra_progs/test_eigen.c,v 1.5 2006-11-29 09:09:31 rotor Exp $";
#endif

#include <volume_io.h>

/* Constants */
#ifndef TRUE
#  define TRUE 1
#  define FALSE 0
#endif

#define VERY_SMALL_EPS 0.0001        /* this is data dependent! */

VIO_BOOL return_local_eigen(VIO_Real r[3][3][3],
                                  VIO_Real dir_1[3],
                                  VIO_Real dir_2[3],
                                  VIO_Real dir_3[3],
                                  VIO_Real val[3]);

VIO_BOOL return_local_eigen_from_hessian(VIO_Real r[3][3][3],
                                  VIO_Real dir_1[3],
                                  VIO_Real dir_2[3],
                                  VIO_Real dir_3[3],
                                  VIO_Real val[3]);

VIO_BOOL return_3D_disp_from_quad_fit(VIO_Real r[3][3][3], 
                                            VIO_Real *dispu, 
                                            VIO_Real *dispv, 
                                            VIO_Real *dispw);        

VIO_BOOL return_3D_disp_from_min_quad_fit(VIO_Real r[3][3][3], 
                                                VIO_Real *dispu, 
                                                VIO_Real *dispv, 
                                                VIO_Real *dispw);        


void setup_val(VIO_Real val[3][3][3]) {
  int i,j,k;

/*
  val[0][0][0] = 0.8;  val[1][0][0] = 0.8;  val[2][0][0] = 0.8;
  val[0][1][0] = 0.8;  val[1][1][0] = 0.8;  val[2][1][0] = 0.8;
  val[0][2][0] = 0.8;  val[1][2][0] = 0.8;  val[2][2][0] = 0.8;

  val[0][0][1] = 0.90;  val[1][0][1] = 0.90;  val[2][0][1] = 0.90;
  val[0][1][1] = 0.90;  val[1][1][1] = 1.0 ;  val[2][1][1] = 0.90;
  val[0][2][1] = 0.90;  val[1][2][1] = 0.90;  val[2][2][1] = 0.90;

  val[0][0][2] = 0.8;  val[1][0][2] = 0.8;  val[2][0][2] = 0.8;
  val[0][1][2] = 0.8;  val[1][1][2] = 0.8;  val[2][1][2] = 0.8;
  val[0][2][2] = 0.8;  val[1][2][2] = 0.8;  val[2][2][2] = 0.8;
*/

  val[0][0][0] = 0.4;  val[1][0][0] = 0.38;  val[2][0][0] = 0.7;
  val[0][1][0] = 0.38;  val[1][1][0] = 0.55;  val[2][1][0] = 0.6;
  val[0][2][0] = 0.7;  val[1][2][0] = 0.6;  val[2][2][0] = 0.7;

  val[0][0][1] = 0.38;  val[1][0][1] = 0.45; val[2][0][1] = 0.6;
  val[0][1][1] = 0.45;  val[1][1][1] = 0.5;  val[2][1][1] = 0.55;
  val[0][2][1] = 0.6;  val[1][2][1] = 0.55;  val[2][2][1] = 0.6;

  val[0][0][2] = 0.7;  val[1][0][2] = 0.6;  val[2][0][2] = 0.7;
  val[0][1][2] = 0.6;  val[1][1][2] = 0.55;  val[2][1][2] = 0.6;
  val[0][2][2] = 0.7;  val[1][2][2] = 0.6;  val[2][2][2] = 0.7;


  for(i=0; i<=2; i++)
    for(j=0; j<=2; j++)
      for(k=0; k<=2; k++) 
        val[i][j][k] = 1.0;

  for(i=0; i<=2; i++) {
    val[0][0][i] *= 1.21;
    val[0][1][i] *= 1.2;
    val[0][2][i] *= 1.21;
    val[2][0][i] *= 1.1;
    val[2][1][i] *= 1.;
    val[2][2][i] *= 1.3;
  }

  val[1][1][1] = 1.0;

}

print_val(VIO_Real r[3][3][3]) 
{
  int u,v,w;

  print ("input data:\n");
  for(v=0; v<3; v++) {
    for(w=0; w<3; w++) {
      for(u=0; u<3; u++)
        print ("%5.3f ",r[u][v][w]);
      print ("    ");
    } 
    print ("\n");
  }
  print ("\n");

}

VIO_BOOL return_principal_directions(VIO_Real r[3][3][3],
                                           VIO_Real dir_1[3],
                                           VIO_Real dir_2[3],
                                           VIO_Real *r_K,
                                           VIO_Real *r_S,
                                           VIO_Real *r_k1,
                                           VIO_Real *r_k2,
                                           VIO_Real *r_norm,
                                           VIO_Real *r_Lvv,
                                           VIO_Real eps);


char *prog_name;

int main(int argc, char *argv[])
{

  VIO_Real 
    tmp, max_val, min_val, intensity_threshold,
    max_val_x, max_val_y, max_val_z,
    min_val_x, min_val_y, min_val_z,
    du, dv, dw,
    K, S, k1, k2, Lvv,
    dir_max[3],
    dir_mid[3],
    dir_min[3],
    dir_vals[3],
    val[3][3][3];
  
  int
    flag,i,j,k;


  prog_name = argv[0];

  if (argc!=1) {
    print("usage: %s \n", prog_name);
    exit(EXIT_FAILURE);
  }


                /* eigen value analysis on Hessian */

  print ("\n***************** Eigen value analysis on Hessian matrix\n");
  for(i=0; i<3; i++)  {
    dir_min[i] = dir_mid[i] = dir_max[i] = dir_vals[i] = -1000000.0;
  }
  setup_val(val);
  print_val(val);
  flag = return_local_eigen_from_hessian(val, dir_min, dir_mid, dir_max, dir_vals);

  if (flag) 
    print ("flag is True (covar matrix is pos semidef)\n"); 
  else 
    print ("flag is False (covar is not pos semidef, or no eigen vals found)\n");

  print ("val:");
  for(i=0; i<3; i++)
    print ("%12.8f ",dir_vals[i]);
  print ("\n");

  print ("vecs:\n");
  for(i=0; i<3; i++) 
    print ("     %12.8f %12.8f %12.8f\n",dir_min[i], dir_mid[i],dir_max[i]);

  if (flag) {
    flag = return_3D_disp_from_min_quad_fit(val, &du, &dv, &dw);

    if (flag) 
      print ("flag is True (covar matrix is pos def)\n"); 
    else 
      print ("flag is False (covar is not pos def, or no eigen vals found)\n");

    print ("disp: %f %f %f\n", du, dv, dw);
  }



                /* eigen value analysis on intensities */

  print ("\n***************** Eigen value analysis on intensities\n");

  for(i=0; i<3; i++)  {
    dir_min[i] = dir_mid[i] = dir_max[i] = dir_vals[i] = -1000000.0;
  }

  setup_val(val);
  flag = return_local_eigen(val, dir_min, dir_mid, dir_max, dir_vals);

  if (flag) print ("flag is True\n"); else print ("flag is False\n");

  for(i=0; i<3; i++)
    print ("%12.8f ",dir_vals[i]);
  print ("\n");

  for(i=0; i<3; i++) 
    print ("%12.8f %12.8f %12.8f\n",dir_min[i], dir_mid[i],dir_max[i]);


  
    
  print ("\n***************** Principal direction analysis\n");

  for(i=0; i<3; i++)  {
    dir_min[i] = dir_mid[i] = dir_max[i] = dir_vals[i] = -1000000;
  }
  setup_val(val);

  flag = return_principal_directions(val,dir_mid,dir_min, &K, &S, 
                                     &k1, &k2, dir_max, &Lvv, 0.0000001);
  
  dir_vals[0] = k2;
  dir_vals[1] = k1;
  dir_vals[2] = sqrt(dir_max[0]*dir_max[0] + 
                     dir_max[1]*dir_max[1] + 
                     dir_max[2]*dir_max[2]);
  
  if (dir_vals[2] > 0.0) {
    for(i=0; i<3; i++)
      dir_max[i] /= dir_vals[2];
  }
  
  if (flag) print ("flag is True\n"); else print ("flag is False\n");
  for(i=0; i<3; i++)
    print ("%12.8f ",dir_vals[i]);
  print ("\n");

  for(i=0; i<3; i++) 
    print ("%12.8f %12.8f %12.8f\n",dir_min[i], dir_mid[i],dir_max[i]);


  exit(EXIT_SUCCESS);
}

