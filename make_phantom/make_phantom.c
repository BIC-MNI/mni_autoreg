/* ----------------------------- MNI Header -----------------------------------
   @NAME       : make_phantom.c
   @INPUT      : data used to define an object in a minc volume
                 
       usage:   make_phantom [options] outputfile.mnc
   
   @OUTPUT     : volume data containing either a voxelated ellipse or rectangle.

   @RETURNS    : TRUE if ok, VIO_ERROR if error.

   @COPYRIGHT  :
              Copyright 1993 Louis Collins, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
   @CREATED    : Wed Mar 16 20:20:50 EST 1994  Louis Collins
   @MODIFIED   : $Log: make_phantom.c,v $
   @MODIFIED   : Revision 1.8  2011-12-08 01:53:04  rotor
   @MODIFIED   :  * added a fix from Vlad to handle cases where an object is on the boundary
   @MODIFIED   :
   @MODIFIED   : Revision 1.7  2006/11/28 08:57:13  rotor
   @MODIFIED   :  * many changes to allow clean build against minc 2.0
   @MODIFIED   :
   @MODIFIED   : Revision 1.6  2005/07/20 20:45:32  rotor
   @MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
   @MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
   @MODIFIED   :     * Removed old VOLUME_IO cruft #defines
   @MODIFIED   :     * Fixed up all Makefile.am's in subdirs
   @MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
   @MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
   @MODIFIED   :
   @MODIFIED   : Revision 1.5  2004/02/12 05:53:40  rotor
   @MODIFIED   :  * removed /static defs
   @MODIFIED   :
   @MODIFIED   : Revision 1.4  2000/02/20 04:00:58  stever
   @MODIFIED   : * use new history_string() function to generate history strings
   @MODIFIED   :   when outputting MNI files (.mnc, .xfm)
   @MODIFIED   : * removed unused vax routines from Proglib
   @MODIFIED   : * tuned configure script; CPPFLAGS and LDFLAGS are now left alone,
   @MODIFIED   :   for the installer to use
   @MODIFIED   :
   @MODIFIED   : Revision 1.3  2000/02/15 19:01:59  stever
   @MODIFIED   : Add tests for param2xfm, minctracc -linear.
   @MODIFIED   :
   @MODIFIED   : Revision 1.2  2000/02/02 20:10:13  stever
   @MODIFIED   : * minctracc/Testing was a copy of Testing with one extra test only;
   @MODIFIED   :   folded the extra test into Testing, and removed minctracc/Testing
   @MODIFIED   : * minor source changes to placate GCC's -Wall option
   @MODIFIED   :
   @MODIFIED   : Revision 1.1  1996/08/21 15:25:05  louis
   @MODIFIED   : Initial revision
   @MODIFIED   :
 * Revision 1.1  94/03/16  22:58:40  louis
 * Initial revision
 * 
---------------------------------------------------------------------------- */

#include <config.h>
#include <float.h>
#include <limits.h>
#include <volume_io.h>
#include <Proglib.h>
#include <ParseArgv.h>
#include "make_phantom.h"

static char *default_dim_names[VIO_N_DIMENSIONS] = { MIxspace, MIyspace, MIzspace };

#define SUBSTEPS  9

#ifndef ROUND
#define  ROUND( x )     (int)floor( (double) (x) + 0.5 )
#endif

void  set_min_max(VIO_Real *voxel, 
                          int *minx, int *maxx, 
                          int *miny, int *maxy, 
                          int *minz, int *maxz)
{
  int x,y,z;

  /* If voxel[VIO_X] = 0.6, then *minx will be set to 1.
     Is that correct, or do we want *minx = 0 ? */
  x = VIO_ROUND(voxel[VIO_X]);
  y = VIO_ROUND(voxel[VIO_Y]);
  z = VIO_ROUND(voxel[VIO_Z]);
  
  if (x > *maxx) *maxx = x;
  if (x < *minx) *minx = x;

  if (y > *maxy) *maxy = y;
  if (y < *miny) *miny = y;

  if (z > *maxz) *maxz = z;
  if (z < *minz) *minz = z;
}

void  set_min_max_r(VIO_Real *voxel, 
                            VIO_Real *minx, VIO_Real *maxx, 
                            VIO_Real *miny, VIO_Real *maxy, 
                            VIO_Real *minz, VIO_Real *maxz)
{
  if (voxel[VIO_X] > *maxx) *maxx = voxel[VIO_X];
  if (voxel[VIO_X] < *minx) *minx = voxel[VIO_X];

  if (voxel[VIO_Y] > *maxy) *maxy = voxel[VIO_Y];
  if (voxel[VIO_Y] < *miny) *miny = voxel[VIO_Y];

  if (voxel[VIO_Z] > *maxz) *maxz = voxel[VIO_Z];
  if (voxel[VIO_Z] < *minz) *minz = voxel[VIO_Z];
}


VIO_Real partial_elliptical(VIO_Volume data, VIO_Real center[], VIO_Real r2[], VIO_Real step[], 
                               int vx, int vy, int vz)
{
  VIO_Real istep,jstep,kstep,coord[3];
  int total,total_in,i,j,k,m;
  VIO_Real frac;

  total = total_in = 0;

  for(i = 0; i <= SUBSTEPS; ++i) {

    istep = -0.5 + ((VIO_Real)i/(VIO_Real)SUBSTEPS);

      for(j = 0; j <= SUBSTEPS; ++j) {

      jstep = -0.5 + ((VIO_Real)j/(VIO_Real)SUBSTEPS);


      for(k = 0; k <= SUBSTEPS; ++k) {
  
        kstep = -0.5 + ((VIO_Real)k/(VIO_Real)SUBSTEPS);

        convert_3D_voxel_to_world(data, 
                                  ((VIO_Real)vx+istep), ((VIO_Real)vy+jstep), ((VIO_Real)vz+kstep), 
                                  &coord[VIO_X], &coord[VIO_Y], &coord[VIO_Z]);
        
   for(m=0; m<3; ++m){
           coord[m] = center[m] - coord[m];
           coord[m] = coord[m]*coord[m] / r2[m];
      }
        
        if ( (coord[VIO_X] + coord[VIO_Y] + coord[VIO_Z]) <= 1.0000001 ) 
          total_in++;

        total++;
      }
    }
  }

  if (total > 0)
    frac = (VIO_Real)total_in/(VIO_Real)total;
  else
    frac = 0.0;


/*
if (frac>0.0);
  print("tot_in = %3d, tot = %3d, frac = %f\n",total_in, total, frac);
*/


  return(frac);
}


VIO_BOOL is_inside_ellipse(VIO_Volume data, VIO_Real center[], VIO_Real r2[], int i, int j, int k)
{
  VIO_Real coord[3];
  int m;

  convert_3D_voxel_to_world(data, 
                            (VIO_Real)i, (VIO_Real)j, (VIO_Real)k, 
                            &coord[VIO_X], &coord[VIO_Y], &coord[VIO_Z]);

  for(m=0; m<3; ++m){
     coord[m] = center[m] - coord[m];
     coord[m] = coord[m]*coord[m] / r2[m];
  }
  
  return (coord[VIO_X] + coord[VIO_Y] + coord[VIO_Z]) <= 1.0000001;
}


int main (int argc, char *argv[] )
{   
  FILE 
    *ofd;
  char 
    *outfilename;  
  VIO_Status 
    status;
  
  VIO_Volume
    data;
  VIO_Real
    fraction,zero, one,
    r2[3],
    coord[3],
    voxel[3];
  char *history;
  VIO_Real
    w_minz,w_maxz,
    w_miny,w_maxy,
    w_minx,w_maxx;
  int 
    in1,in2,in3,in4,in5,in6,in7,
    minz,maxz,
    miny,maxy,
    minx,maxx,
    i,j,k,m;

  /* set default values */

  prog_name = argv[0];
  outfilename = NULL;

  history = history_string(argc, argv);

  /* Call ParseArgv to interpret all command line args */
  if (ParseArgv(&argc, argv, argTable, 0) || (argc!=2)) {
    (void) fprintf(stderr, 
                   "\nUsage: %s [<options>] output.mnc\n", 
                   prog_name);
    (void) fprintf(stderr,"       %s [-help]\n\n", prog_name);
    exit(EXIT_FAILURE);
  }
  
  /* set up necessary file names */
  outfilename  = argv[1];
  
  if (!clobber_flag && file_exists(outfilename)) {
    print ("File %s exists.\n",outfilename);
    print ("Use -clobber to overwrite.\n");
    return VIO_ERROR;
  }

  /* check to see if the output file can be written */
  status = open_file( outfilename , WRITE_FILE, BINARY_FORMAT,  &ofd );
  if ( status != VIO_OK ) {
    print ("filename `%s' cannot be opened.", outfilename);
    return VIO_ERROR;
  }
  status = close_file(ofd);
  (void)remove(outfilename);   

  /******************************************************************************/
  /*             create volume data                                             */
  /******************************************************************************/
  
  data = create_volume(3, default_dim_names, datatype, is_signed, 0.0, 0.0);
  set_volume_voxel_range(data, voxel_range[0], voxel_range[1]);
  set_volume_real_range( data, real_range[0],  real_range[1]);
  set_volume_sizes(data, count);
  set_volume_separations(data, step);

  voxel[0] = 0.0;
  voxel[1] = 0.0;
  voxel[2] = 0.0;

  set_volume_translation(data, voxel, start);
  alloc_volume_data(data);

  zero = CONVERT_VALUE_TO_VOXEL(data, background);
  one = CONVERT_VALUE_TO_VOXEL(data, fill_value);

  if (debug) print ("zero: real = %f, voxel = %f\n", background, zero);
  if (debug) print ("one: real = %f, voxel = %f\n", fill_value, one);

  /******************************************************************************/
  /*             write out background value                                     */
  /******************************************************************************/

  for(i=count[VIO_X]; i--; )
    for(j=count[VIO_Y]; j--; )
      for(k=count[VIO_Z]; k--; ){
        SET_VOXEL_3D(data, i,j,k, zero);
      }

  /******************************************************************************/
  /*             build object data                                              */
  /******************************************************************************/

  /* begin by finding bounding box for data */
  maxx = maxy = maxz = -INT_MAX;
  minx = miny = minz = INT_MAX;
  w_maxx = w_maxy = w_maxz = -DBL_MAX;
  w_minx = w_miny = w_minz = DBL_MAX;

/* sometimes we may want the object to be outside the volume
  for(i=3; i<3; i++){
     coord[i] = center[i];
     }
  convert_3D_world_to_voxel(data, 
                            coord[VIO_X], coord[VIO_Y], coord[VIO_Z], 
                            &voxel[VIO_X], &voxel[VIO_Y], &voxel[VIO_Z]);
  set_min_max(voxel, &minx, &maxx, &miny, &maxy, &minz, &maxz);
  set_min_max_r(coord, &w_minx, &w_maxx, &w_miny, &w_maxy, &w_minz, &w_maxz);
  if (debug) print ("%f %f %f\n",voxel[VIO_X], voxel[VIO_Y], voxel[VIO_Z]);
*/

  for(i=0; i<3; i++){
     coord[i] = center[i];
     }
  coord[VIO_Z] += width[VIO_Z]/2.0;
  
  convert_3D_world_to_voxel(data, 
                            coord[VIO_X], coord[VIO_Y], coord[VIO_Z], 
                            &voxel[VIO_X], &voxel[VIO_Y], &voxel[VIO_Z]);  
  if (debug) print ("%f %f %f\n",voxel[VIO_X], voxel[VIO_Y], voxel[VIO_Z]);
  set_min_max(voxel, &minx, &maxx, &miny, &maxy, &minz, &maxz);
  set_min_max_r(coord, &w_minx, &w_maxx, &w_miny, &w_maxy, &w_minz, &w_maxz);

  for(i=0; i<3; i++){
     coord[i] = center[i];
     }
  coord[VIO_Z] -= width[VIO_Z]/2.0;
  convert_3D_world_to_voxel(data, 
                            coord[VIO_X], coord[VIO_Y], coord[VIO_Z], 
                            &voxel[VIO_X], &voxel[VIO_Y], &voxel[VIO_Z]);  
  if (debug) print ("%f %f %f\n",voxel[VIO_X], voxel[VIO_Y], voxel[VIO_Z]);
  set_min_max(voxel, &minx, &maxx, &miny, &maxy, &minz, &maxz);
  set_min_max_r(coord, &w_minx, &w_maxx, &w_miny, &w_maxy, &w_minz, &w_maxz);
  
  for(i=0; i<3; i++){
     coord[i] = center[i];
     }
  coord[VIO_Y] += width[VIO_Y]/2.0;
  convert_3D_world_to_voxel(data, 
                            coord[VIO_X], coord[VIO_Y], coord[VIO_Z], 
                            &voxel[VIO_X], &voxel[VIO_Y], &voxel[VIO_Z]);  
  if (debug) print ("%f %f %f\n",voxel[VIO_X], voxel[VIO_Y], voxel[VIO_Z]);
  set_min_max(voxel, &minx, &maxx, &miny, &maxy, &minz, &maxz);
  set_min_max_r(coord, &w_minx, &w_maxx, &w_miny, &w_maxy, &w_minz, &w_maxz);

  for(i=0; i<3; i++){
     coord[i] = center[i];
     }
  coord[VIO_Y] -= width[VIO_Y]/2.0;
  convert_3D_world_to_voxel(data, 
                            coord[VIO_X], coord[VIO_Y], coord[VIO_Z], 
                            &voxel[VIO_X], &voxel[VIO_Y], &voxel[VIO_Z]);  
  if (debug) print ("%f %f %f\n",voxel[VIO_X], voxel[VIO_Y], voxel[VIO_Z]);
  set_min_max(voxel, &minx, &maxx, &miny, &maxy, &minz, &maxz);
  set_min_max_r(coord, &w_minx, &w_maxx, &w_miny, &w_maxy, &w_minz, &w_maxz);

  for(i=0; i<3; i++){
     coord[i] = center[i];
     }
  coord[VIO_X] += width[VIO_X]/2.0;
  convert_3D_world_to_voxel(data, 
                            coord[VIO_X], coord[VIO_Y], coord[VIO_Z], 
                            &voxel[VIO_X], &voxel[VIO_Y], &voxel[VIO_Z]);  
  if (debug) print ("%f %f %f\n",voxel[VIO_X], voxel[VIO_Y], voxel[VIO_Z]);
  set_min_max(voxel, &minx, &maxx, &miny, &maxy, &minz, &maxz);
  set_min_max_r(coord, &w_minx, &w_maxx, &w_miny, &w_maxy, &w_minz, &w_maxz);

  for(i=0; i<3; i++){
     coord[i] = center[i];
     }
  coord[VIO_X] -= width[VIO_X]/2.0;
  convert_3D_world_to_voxel(data, 
                            coord[VIO_X], coord[VIO_Y], coord[VIO_Z], 
                            &voxel[VIO_X], &voxel[VIO_Y], &voxel[VIO_Z]);  
  if (debug) print ("%f %f %f\n",voxel[VIO_X], voxel[VIO_Y], voxel[VIO_Z]);
  set_min_max(voxel, &minx, &maxx, &miny, &maxy, &minz, &maxz);
  set_min_max_r(coord, &w_minx, &w_maxx, &w_miny, &w_maxy, &w_minz, &w_maxz);

  if (debug) print ("mins, maxs: %3d %3d %3d %3d %3d %3d\n", minx, maxx, miny, maxy, minz, maxz);


  /* check to see if limits of object fit in data volume */
  if ( maxx < 0  || maxy < 0 || maxz < 0 ||
       minx > count[VIO_X] || miny > count[VIO_Y] || minz > count[VIO_Z]) {
    print ("There is no voxel overlap between object and volume as defined.\n");
    return VIO_ERROR;
  }

  /* add one voxel all around, for partial volume calculations. */
  if (object==ELLIPSE) {
    --minx;  --miny;  --minz;
    ++maxx;  ++maxy;  ++maxz;
  }
  
  /* readjust for edge of volume */
  if (minx < 0) minx = 0;
  if (maxx > count[VIO_X]-1) maxx = count[VIO_X]-1;
  if (miny < 0) miny = 0;
  if (maxy > count[VIO_Y]-1) maxy = count[VIO_Y]-1;
  if (minz < 0) minz = 0;
  if (maxz > count[VIO_Z]-1) maxz = count[VIO_Z]-1;


  if (debug) print ("mins, maxs: %3d %3d %3d %3d %3d %3d\n", minx, maxx, miny, maxy, minz, maxz);

  if (object==RECTANGLE) {
    for(i=minx; i<=maxx; i++){
      for(j=miny; j<=maxy; j++){
         for(k=minz; k<=maxz; k++){
            SET_VOXEL_3D(data, i,j,k, one);
            }
         }
      }
  } else {                        /* ellipsoid */

    for(m=0; m<3; m++){
      r2[m] = width[m]*width[m]/4.0;
      }

    for(i=minx; i<=maxx; i++){
      for(j=miny; j<=maxy; j++){
         for(k=minz; k<=maxz; k++){

          in1 = is_inside_ellipse(data, center, r2, i+1,j,  k );
          in2 = is_inside_ellipse(data, center, r2, i-1,j,  k);
          in3 = is_inside_ellipse(data, center, r2, i,  j+1,k);
          in4 = is_inside_ellipse(data, center, r2, i,  j-1,k);
          in5 = is_inside_ellipse(data, center, r2, i,  j,  k+1);
          in6 = is_inside_ellipse(data, center, r2, i,  j,  k-1);
          in7 = is_inside_ellipse(data, center, r2, i,  j,  k);
          
          if (in1 && in2 && in3 && in4 && in5 && in6 && in7) {
             SET_VOXEL_3D(data, i,j,k, one);
             }
          else {
            if (in1 || in2 || in3 || in4 || in5 || in6 || in7) {

              if (partial_flag) {
                VIO_Real value;
                fraction = partial_elliptical(data, center, r2, step, i,  j,  k);
                fraction = fraction*edge_value;
                value = CONVERT_VALUE_TO_VOXEL(data, fraction);
                SET_VOXEL_3D(data, i,j,k, value);
              }
              else {
                  if (in7) {
                      SET_VOXEL_3D(data, i,j,k, one);
                  }
              }
              
            }
          }
         }
      }          
    }
    
  }
  
  
  status = output_volume(outfilename, NC_UNSPECIFIED, FALSE, 0.0, 0.0, data, 
                         history, (minc_output_options *)NULL);  
  if (status != VIO_OK)
    print("problems writing volume data for %s.", outfilename);
  
  
  return(status);
   
}

