/* ----------------------------- MNI Header -----------------------------------
@NAME       : interpolation.c
@DESCRIPTION: File containing routines to interpolate voxel values
              from minc volumes using different interpolation kernels.
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

@CREATED    : Wed May 26 13:05:44 EST 1993 LC using routines from NEELIN's
              mincresample.
@MODIFIED   :  $Log: interpolation.c,v $
@MODIFIED   :  Revision 96.7  2006-11-30 09:07:33  rotor
@MODIFIED   :   * many more changes for clean minc 2.0 build
@MODIFIED   :
@MODIFIED   :  Revision 96.5  2005/07/20 20:45:52  rotor
@MODIFIED   :      * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :      * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :      * Removed old VOLUME_IO cruft #defines
@MODIFIED   :      * Fixed up all Makefile.am's in subdirs
@MODIFIED   :      * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :      * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   :  Revision 96.4  2004/02/12 06:04:53  rotor
@MODIFIED   :   * removed /static defs
@MODIFIED   :
@MODIFIED   :  Revision 96.3  2002/03/26 14:15:47  stever
@MODIFIED   :  Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   :  Revision 96.2  2000/03/17 01:11:31  stever
@MODIFIED   :  code simplification
@MODIFIED   :
@MODIFIED   :  Revision 96.1  1999/10/25 19:59:17  louis
@MODIFIED   :  final checkin before switch to CVS
@MODIFIED   :
 * Revision 96.0  1996/08/21  18:22:15  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.5  1996/08/12  14:16:15  louis
 * Release of MNI_AutoReg version 1.0
 *
 * Revision 1.10  1996/08/12  14:16:13  louis
 * Pre-release
 *
 * Revision 1.9  1995/02/22  08:56:06  collins
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.8  94/04/06  11:48:38  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.7  94/02/21  16:35:39  louis
 * version before feb 22 changes
 * 
 * Revision 1.6  93/11/15  16:26:46  louis
 * working version, with new library, with RCS revision stuff,
 * before deformations included
 * 
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/minctracc/Volume/interpolation.c,v 96.7 2006-11-30 09:07:33 rotor Exp $";
#endif

#include <volume_io.h>
#include "minctracc_point_vector.h"

#define VOL_NDIMS 3

/* ----------------------------- MNI Header -----------------------------------
@NAME       : nearest_neighbour_interpolant
@INPUT      : volume - pointer to volume data
              coord - point at which volume should be interpolated in voxel 
                 units (with 0 being first point of the volume).
@OUTPUT     : result - interpolated value.
@RETURNS    : TRUE if coord is within the volume, FALSE otherwise.
@DESCRIPTION: Routine to interpolate a volume at a point with nearest
              neighbour interpolation.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : February 12, 1993 (Peter Neelin)
@MODIFIED   : Fri May 28 09:06:12 EST 1993 Louis Collins
               mod to use david's volume_struct
---------------------------------------------------------------------------- */
int nearest_neighbour_interpolant(VIO_Volume volume, 
                                         PointR *coord, double *result)
{
   long 
     ind0, ind1, ind2, 
     max[3];
   int sizes[3];

   /* Check that the coordinate is inside the volume */
   
  get_volume_sizes(volume, sizes);
  max[0]=sizes[0];
  max[1]=sizes[1];
  max[2]=sizes[2];
   
   if ((Point_x( *coord ) < -0.5) || (Point_x( *coord ) >= max[0]-0.5) ||
       (Point_y( *coord ) < -0.5) || (Point_y( *coord ) >= max[1]-0.5) ||
       (Point_z( *coord ) < -0.5) || (Point_z( *coord ) >= max[2]-0.5)) {

      *result = CONVERT_VOXEL_TO_VALUE( volume, get_volume_voxel_min(volume));

      return FALSE;
   }

   /* Get the whole part of the coordinate */
   ind0 = (long) floor(Point_x( *coord ) + 0.5);
   ind1 = (long) floor(Point_y( *coord ) + 0.5);
   ind2 = (long) floor(Point_z( *coord ) + 0.5);

   /* Get the value */
   GET_VALUE_3D( *result ,  volume, ind0  , ind1  , ind2  );

   return TRUE;

}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : trilinear_interpolant
@INPUT      : volume - pointer to volume data
              coord - point at which volume should be interpolated in voxel 
                 units (with 0 being first point of the volume).
@OUTPUT     : result - interpolated TRUE value.
@RETURNS    : TRUE if coord is within the volume, FALSE otherwise.
@DESCRIPTION: Routine to interpolate a volume at a point with tri-linear
              interpolation.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : February 10, 1993 (Peter Neelin)
@MODIFIED   : Fri May 28 09:06:12 EST 1993 Louis Collins
               mod to use david's volume_struct
---------------------------------------------------------------------------- */
int trilinear_interpolant(VIO_Volume volume, 
                                 PointR *coord, double *result)
{
  long ind0, ind1, ind2, max[3];
  int sizes[3];
  int flag;
  double temp_result;
  static double f0, f1, f2, r0, r1, r2, r1r2, r1f2, f1r2, f1f2;
  static double v000, v001, v010, v011, v100, v101, v110, v111;
  
  /* Check that the coordinate is inside the volume */
  
  get_volume_sizes(volume, sizes);
  max[0]=sizes[0];
  max[1]=sizes[1];
  max[2]=sizes[2];
  
  if ((Point_x( *coord ) < 0) || (Point_x( *coord ) >= max[0]-1) ||
      (Point_y( *coord ) < 0) || (Point_y( *coord ) >= max[1]-1) ||
      (Point_z( *coord ) < 0) || (Point_z( *coord ) >= max[2]-1)) {
    
    flag = nearest_neighbour_interpolant(volume, coord, &temp_result) ;
    *result = temp_result;
    return(flag);
  }
    
  /* Get the whole part of the coordinate */ 
  ind0 = (long) floor(Point_x( *coord ));
  ind1 = (long) floor(Point_y( *coord ));
  ind2 = (long) floor(Point_z( *coord ));
  if (ind0 >= max[0]-1) ind0 = max[0]-1;
  if (ind1 >= max[1]-1) ind1 = max[1]-1;
  if (ind2 >= max[2]-1) ind2 = max[2]-1;
  
  /* Get the relevant voxels */
  GET_VALUE_3D( v000 ,  volume, ind0  , ind1  , ind2   ); 
  GET_VALUE_3D( v001 ,  volume, ind0  , ind1  , ind2+1 ); 
  GET_VALUE_3D( v010 ,  volume, ind0  , ind1+1, ind2   ); 
  GET_VALUE_3D( v011 ,  volume, ind0  , ind1+1, ind2+1 ); 
  GET_VALUE_3D( v100 ,  volume, ind0+1, ind1  , ind2   ); 
  GET_VALUE_3D( v101 ,  volume, ind0+1, ind1  , ind2+1 ); 
  GET_VALUE_3D( v110 ,  volume, ind0+1, ind1+1, ind2   ); 
  GET_VALUE_3D( v111 ,  volume, ind0+1, ind1+1, ind2+1 ); 

  /* Get the fraction parts */
  f0 = Point_x( *coord ) - ind0;
  f1 = Point_y( *coord ) - ind1;
  f2 = Point_z( *coord ) - ind2;
  r0 = 1.0 - f0;
  r1 = 1.0 - f1;
  r2 = 1.0 - f2;
  
  /* Do the interpolation */
  r1r2 = r1 * r2;
  r1f2 = r1 * f2;
  f1r2 = f1 * r2;
  f1f2 = f1 * f2;
  
  *result =
    r0 *  (r1r2 * v000 +
           r1f2 * v001 +
           f1r2 * v010 +
           f1f2 * v011);
  *result +=
    f0 *  (r1r2 * v100 +
           r1f2 * v101 +
           f1r2 * v110 +
           f1f2 * v111);

  
  return TRUE;
  
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : do_Ncubic_interpolation
@INPUT      : volume - pointer to volume data
              index - indices to start point in volume
                 (bottom of 4x4x4 cube for interpolation)
              cur_dim - dimension to be interpolated (0 = volume, 1 = slice,
                 2 = line)
@OUTPUT     : result - interpolated value.
@RETURNS    : (nothing)
@DESCRIPTION: Routine to interpolate a volume, slice or line (specified by
              cur_dim).
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : February 12, 1993 (Peter Neelin)
@MODIFIED   : Fri May 28 09:06:12 EST 1993 Louis Collins
               mod to use david's volume_struct
---------------------------------------------------------------------------- */
void do_Ncubic_interpolation(VIO_Volume volume, 
                             long index[], int cur_dim, 
                             double frac[], double *result)
{
   long base_index;
   double v0, v1, v2, v3, u;

   /* Save index that we will change */
   base_index = index[cur_dim];

   /* If last dimension, then just get the values */
   if (cur_dim == VOL_NDIMS-1) {
     GET_VOXEL_3D( v0 ,  volume, index[0] ,index[1], index[2] );
     index[cur_dim]++;
     GET_VOXEL_3D( v1 ,  volume, index[0] ,index[1], index[2] );
     index[cur_dim]++;
     GET_VOXEL_3D( v2 ,  volume, index[0] ,index[1], index[2] );
     index[cur_dim]++;
     GET_VOXEL_3D( v3 ,  volume, index[0] ,index[1], index[2] );
   }

   /* Otherwise, recurse */
   else {
      do_Ncubic_interpolation(volume, index, cur_dim+1, frac, &v0);
      index[cur_dim]++;
      do_Ncubic_interpolation(volume, index, cur_dim+1, frac, &v1);
      index[cur_dim]++;
      do_Ncubic_interpolation(volume, index, cur_dim+1, frac, &v2);
      index[cur_dim]++;
      do_Ncubic_interpolation(volume, index, cur_dim+1, frac, &v3);
   }

   /* Restore index */
   index[cur_dim] = base_index;

   /* Scale values for slices */
   if (cur_dim==0) {
      v0 = CONVERT_VOXEL_TO_VALUE(volume,v0);
      v1 = CONVERT_VOXEL_TO_VALUE(volume,v1);
      v2 = CONVERT_VOXEL_TO_VALUE(volume,v2);
      v3 = CONVERT_VOXEL_TO_VALUE(volume,v3);
   }

   /* Get fraction */
   u = frac[cur_dim];

   /* Do tricubic interpolation (code from Dave MacDonald).
      Gives v1 and v2 at u = 0 and 1 and gives continuity of intensity
      and first derivative. */
   *result =
     ( (v1) + (u) * (
       0.5 * ((v2)-(v0)) + (u) * (
       (v0) - 2.5 * (v1) + 2.0 * (v2) - 0.5 * (v3) + (u) * (
       -0.5 * (v0) + 1.5 * (v1) - 1.5 * (v2) + 0.5 * (v3)  )
                                 )
                    )
     );

   return;
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : tricubic_interpolant
@INPUT      : volume - pointer to volume data
              coord - point at which volume should be interpolated in voxel 
                 units (with 0 being first point of the volume).
@OUTPUT     : result - interpolated value.
@RETURNS    : TRUE if coord is within the volume, FALSE otherwise.
@DESCRIPTION: Routine to interpolate a volume at a point with tri-cubic
              interpolation.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : February 12, 1993 (Peter Neelin)
@MODIFIED   : Fri May 28 09:06:12 EST 1993 Louis Collins
               mod to use david's volume_struct
---------------------------------------------------------------------------- */
int tricubic_interpolant(VIO_Volume volume, 
                                PointR *coord, double *result)
{
   long ind0, ind1, ind2, max[3], index[VOL_NDIMS];
   double frac[VOL_NDIMS];
   int sizes[3];
   int flag;
   double temp_result;

   /* Check that the coordinate is inside the volume */

   get_volume_sizes(volume, sizes);
   max[0] = sizes[0];
   max[1] = sizes[1];
   max[2] = sizes[2];
   
   if ((Point_x( *coord ) < 0) || (Point_x( *coord ) >= max[0]-1) ||
       (Point_y( *coord ) < 0) || (Point_y( *coord ) >= max[1]-1) ||
       (Point_z( *coord ) < 0) || (Point_z( *coord ) >= max[2]-1)) {

     flag = nearest_neighbour_interpolant(volume, coord, &temp_result) ;
     *result = temp_result;
     return(flag);
   }


   /* Get the whole and fractional part of the coordinate */
   ind0 = (long) floor(Point_x( *coord ));
   ind1 = (long) floor(Point_y( *coord ));
   ind2 = (long) floor(Point_z( *coord ));
   frac[0] = Point_x( *coord ) - ind0;
   frac[1] = Point_y( *coord ) - ind1;
   frac[2] = Point_z( *coord ) - ind2;
   ind0--;
   ind1--;
   ind2--;

   /* Check for edges - do linear interpolation at edges */
   if ((ind0 >= max[0]-3) || (ind0 < 0) ||
       (ind1 >= max[1]-3) || (ind1 < 0) ||
       (ind2 >= max[2]-3) || (ind2 < 0)) {
      return trilinear_interpolant(volume, coord, result);
   }
   index[0]=ind0; index[1]=ind1; index[2]=ind2;

   /* Do the interpolation */
   do_Ncubic_interpolation(volume, index, 0, frac, result);

   return TRUE;

}


/* A point is not masked if it is a point we should consider.
   If the mask volume is NULL, we consider all points.
   Otherwise, consider a point if the mask volume value is > 0.
*/
int point_not_masked( VIO_Volume volume, 
                             VIO_Real wx, VIO_Real wy, VIO_Real wz)
{
    double result;
    PointR coord;

    if ( volume == NULL )
        return TRUE;

    convert_3D_world_to_voxel( volume, wx, wy, wz, 
                               &Point_x(coord), &Point_y(coord), &Point_z(coord) );
    
    /* interpolation returns TRUE iff coordinate is inside volume */
    if ( nearest_neighbour_interpolant(volume,&coord,&result) ) {
        return (result > 0.0);
    }

    return(FALSE) ;
}


int voxel_point_not_masked( VIO_Volume volume, 
                                   VIO_Real vx, VIO_Real vy, VIO_Real vz)
{
    double result;
    PointR coord;
  
    if ( volume == NULL )
        return TRUE;

    fill_Point(coord, vx, vy, vz);
    
    if ( nearest_neighbour_interpolant(volume,&coord,&result) ) {
        return (result > 0.0);
    }

    return(FALSE) ;
}

