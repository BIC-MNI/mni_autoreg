/* ----------------------------- MNI Header -----------------------------------
@NAME       : interpolation.c
@DESCRIPTION: File containing routines to interpolate voxel values
              from minc volumes using different interpolation kernels.
@METHOD     : 
@GLOBALS    : 
@CREATED    : Wed May 26 13:05:44 EST 1993 LC using routines from NEELIN's
              mincresample.
@MODIFIED   : 


---------------------------------------------------------------------------- */


#include <def_mni.h>
#include "def_point_vector.h"

#define VOL_NDIMS 3


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
	      Wed Sep 15 09:01:17 EDT 1993 LC
               swapped x & z references.  x varies fastest, z, slowest
---------------------------------------------------------------------------- */
public int trilinear_interpolant(Volume volume, 
                                 PointR *coord, double *result)
{
  long ind0, ind1, ind2, max[3];
  int sizes[3];
  static double f0, f1, f2, r0, r1, r2, r1r2, r1f2, f1r2, f1f2;
  static double v000, v001, v010, v011, v100, v101, v110, v111;
  
  /* Check that the coordinate is inside the volume */
  
  get_volume_sizes(volume, sizes);
  max[0]=sizes[0];
  max[1]=sizes[1];
  max[2]=sizes[2];
  
  if ((Point_z( *coord ) < 0) || (Point_z( *coord ) > max[0]-1) ||
      (Point_y( *coord ) < 0) || (Point_y( *coord ) > max[1]-1) ||
      (Point_x( *coord ) < 0) || (Point_x( *coord ) > max[2]-1)) {
    
    *result = 0.0;
    
    return FALSE;
  }
  
  
  /* Get the whole part of the coordinate */ 
  ind0 = (long) Point_z( *coord );
  ind1 = (long) Point_y( *coord );
  ind2 = (long) Point_x( *coord );
  if (ind0 >= max[0]-1) ind0 = max[0]-1;
  if (ind1 >= max[1]-1) ind1 = max[1]-1;
  if (ind2 >= max[2]-1) ind2 = max[2]-1;
  
  /* Get the relevant voxels */
  GET_VOXEL_3D( v000 ,  volume, ind0  , ind1  , ind2   ); 
  v000=CONVERT_VOXEL_TO_VALUE(volume,v000);
  GET_VOXEL_3D( v001 ,  volume, ind0  , ind1  , ind2+1 ); 
  v001=CONVERT_VOXEL_TO_VALUE(volume,v001);
  GET_VOXEL_3D( v010 ,  volume, ind0  , ind1+1, ind2   ); 
  v010=CONVERT_VOXEL_TO_VALUE(volume,v010);
  GET_VOXEL_3D( v011 ,  volume, ind0  , ind1+1, ind2+1 ); 
  v011=CONVERT_VOXEL_TO_VALUE(volume,v011);
  GET_VOXEL_3D( v100 ,  volume, ind0+1, ind1  , ind2   ); 
  v101=CONVERT_VOXEL_TO_VALUE(volume,v100);
  GET_VOXEL_3D( v101 ,  volume, ind0+1, ind1  , ind2+1 ); 
  v101=CONVERT_VOXEL_TO_VALUE(volume,v101);
  GET_VOXEL_3D( v110 ,  volume, ind0+1, ind1+1, ind2   ); 
  v110=CONVERT_VOXEL_TO_VALUE(volume,v110);
  GET_VOXEL_3D( v111 ,  volume, ind0+1, ind1+1, ind2+1 ); 
  v111=CONVERT_VOXEL_TO_VALUE(volume,v111);

   /* Get the fraction parts */
   f0 = Point_z( *coord ) - ind0;
   f1 = Point_y( *coord ) - ind1;
   f2 = Point_x( *coord ) - ind2;
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
public void do_Ncubic_interpolation(Volume volume, 
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
	      Wed Sep 15 09:01:17 EDT 1993 LC
               swapped x & z references.  x varies fastest, z, slowest
---------------------------------------------------------------------------- */
public int tricubic_interpolant(Volume volume, 
                                PointR *coord, double *result)
{
   long ind0, ind1, ind2, max[3], index[VOL_NDIMS];
   double frac[VOL_NDIMS];
   int sizes[3];

   /* Check that the coordinate is inside the volume */

   get_volume_sizes(volume, sizes);
   max[0] = sizes[0];
   max[1] = sizes[1];
   max[2] = sizes[2];
   
   if ((Point_z( *coord ) < 0) || (Point_z( *coord ) > max[0]-1) ||
       (Point_y( *coord ) < 0) || (Point_y( *coord ) > max[1]-1) ||
       (Point_x( *coord ) < 0) || (Point_x( *coord ) > max[2]-1)) {

      *result = 0.0;

      return FALSE;
   }


   /* Get the whole and fractional part of the coordinate */
   ind0 = (long) Point_z( *coord );
   ind1 = (long) Point_y( *coord );
   ind2 = (long) Point_x( *coord );
   frac[0] = Point_z( *coord ) - ind0;
   frac[1] = Point_y( *coord ) - ind1;
   frac[2] = Point_x( *coord ) - ind2;
   ind0--;
   ind1--;
   ind2--;

   /* Check for edges - do linear interpolation at edges */
   if ((ind0 > max[0]-3) || (ind0 < 0) ||
       (ind1 > max[1]-3) || (ind1 < 0) ||
       (ind2 > max[2]-3) || (ind2 < 0)) {
      return trilinear_interpolant(volume, coord, result);
   }
   index[0]=ind0; index[1]=ind1; index[2]=ind2;

   /* Do the interpolation */
   do_Ncubic_interpolation(volume, index, 0, frac, result);

   return TRUE;

}

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
	      Wed Sep 15 09:01:17 EDT 1993 LC
               swapped x & z references.  x varies fastest, z, slowest
---------------------------------------------------------------------------- */
public int nearest_neighbour_interpolant(Volume volume, 
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
   
   if ((Point_z( *coord ) < -0.5) || (Point_z( *coord ) > max[0]-0.5) ||
       (Point_y( *coord ) < -0.5) || (Point_y( *coord ) > max[1]-0.5) ||
       (Point_x( *coord ) < -0.5) || (Point_x( *coord ) > max[2]-0.5)) {

      *result = 0.0;

      return FALSE;
   }

   /* Get the whole part of the coordinate */
   ind0 = (long) (Point_z( *coord ) + 0.5);
   ind1 = (long) (Point_y( *coord ) + 0.5);
   ind2 = (long) (Point_x( *coord ) + 0.5);

   /* Get the value */
   GET_VOXEL_3D( *result ,  volume, ind0  , ind1  , ind2  );

   *result = CONVERT_VOXEL_TO_VALUE( volume, *result);

   return TRUE;

}


