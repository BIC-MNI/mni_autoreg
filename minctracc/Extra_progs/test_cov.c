#include <volume_io.h>

#include "constants.h"
#include "matrix_basics.h"
#include "make_rots.h"
#include "point_vector.h"

#define DO_TRANSFORM(result, transformation, coord) \
   general_transform_point(transformation, \
      Point_x(coord), Point_y(coord), Point_z(coord), \
      &(Point_x(result)), &(Point_y(result)), &(Point_z(result)) )

#define INTERPOLATE_TRUE_VALUE(volume, coord, result) \
   trilinear_interpolant(volume, coord, result)

static char *default_dim_names[N_DIMENSIONS] =
   { MIzspace, MIyspace, MIxspace };


BOOLEAN vol_to_cov(Volume d1, Volume m1, float *centroid, float **covar, double *step)
{

  VectorR
    vector_step,
    slice_step,
    row_step,
    col_step;

  PointR 
    starting_offset,
    starting_origin,
    starting_position,
    slice,
    row,
    col,
    voxel;

  Real
    tx,ty,tz,
    voxel_value;
  int
    i,r,c,s,
    limits[VOL_NDIMS];

  float
    t,
    sxx,syy,szz,
    sxy,syz,sxz,
    sx,sy,sz,si; 
  Real
    thickness[3];
  int
    sizes[3];

  Real
    true_value,
    position[3];

  get_volume_separations(d1, thickness);
  get_volume_sizes(d1, sizes);
  
				/* build sampling lattice info */
  for_less( i, 0, 3) {	
    step[i] *= thickness[i] / ABS( thickness[i]);
  }

  fill_Vector( col_step,   step[COL_IND], 0.0,     0.0 );
  fill_Vector( row_step,   0.0,     step[ROW_IND], 0.0 );
  fill_Vector( slice_step, 0.0,     0.0,     step[SLICE_IND] );

  convert_3D_voxel_to_world(d1, 0.0, 0.0, 0.0, &tx, &ty, &tz); 

  fill_Point( starting_origin, tx, ty, tz);

  for_less( i, 0, 3) {		/* for each dim, get # of steps in that direction,
				   and set starting offset */
    t = sizes[i] * thickness[i] / step[i];
    limits[i] = (int)( ABS( t ) );
    
    Point_coord( starting_offset, (i) ) = 
      ( (sizes[i]-1)*thickness[i] - (limits[i] * step[i] ) ) / 2.0;
  }
  
  ADD_POINTS( starting_position, starting_origin, starting_offset ); /*  */

				/* calculate centroids first */

  sx = 0.0;
  sy = 0.0;
  sz = 0.0;
  si = 0.0;

  for_inclusive(s,0,limits[SLICE_IND]) {

    SCALE_VECTOR( vector_step, slice_step, s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    for_inclusive(r,0,limits[ROW_IND]) {

      SCALE_VECTOR( vector_step, row_step, r);
      ADD_POINT_VECTOR( row, slice, vector_step );

      SCALE_POINT( col, row, 1.0); /* init first col position */
      for_inclusive(c,0,limits[COL_IND]) {

	convert_3D_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);

	fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
	
	if (point_not_masked(m1, Point_x(col), Point_y(col), Point_z(col))) {
	  
	  if (INTERPOLATE_TRUE_VALUE( d1, &voxel, &true_value )) {
	    
	    sx +=  Point_x(col) * true_value;
	    sy +=  Point_y(col) * true_value;
	    sz +=  Point_z(col) * true_value;
	    
	    si += true_value;
	  }
	  /* else requested voxel is just outside volume., so ignore it */

	} 
	
	ADD_POINT_VECTOR( col, col, col_step );
	
      }
    }
  }

  if (si!=0.0) {
    centroid[1] = sx/ si;
    centroid[2] = sy/ si;
    centroid[3] = sz/ si;
    
    sxx = syy = szz = 0.0;
    sxy = syz = sxz = 0.0;
    
				/* now calculate variances and co-variances */

    for_inclusive(s,0,limits[Z]) {
      
      SCALE_VECTOR( vector_step, slice_step, s);
      ADD_POINT_VECTOR( slice, starting_position, vector_step );
      
      for_inclusive(r,0,limits[Y]) {
	
	SCALE_VECTOR( vector_step, row_step, r);
	ADD_POINT_VECTOR( row, slice, vector_step );
	
	SCALE_POINT( col, row, 1.0); /* init first col position */
	for_inclusive(c,0,limits[X]) {
	  
	  convert_3D_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);

	  fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
	
	  if (point_not_masked(m1, Point_x(col), Point_y(col), Point_z(col))) {
	    
	    if (INTERPOLATE_TRUE_VALUE( d1, &voxel, &true_value )) {
	      
      	      sxx += (Point_x( col )-centroid[1]) * (Point_x( col )-centroid[1]) * true_value;
	      syy += (Point_y( col )-centroid[2]) * (Point_y( col )-centroid[2]) * true_value;
	      szz += (Point_z( col )-centroid[3]) * (Point_z( col )-centroid[3]) * true_value;
	      sxy += (Point_x( col )-centroid[1]) * (Point_y( col )-centroid[2]) * true_value;
	      syz += (Point_y( col )-centroid[2]) * (Point_z( col )-centroid[3]) * true_value;
	      sxz += (Point_x( col )-centroid[1]) * (Point_z( col )-centroid[3]) * true_value;

	    }
	    /* else requested voxel is just outside volume., so ignore it */

	  } 
	  
	  ADD_POINT_VECTOR( col, col, col_step );
	
	}
      }
    }
    
    covar[1][1] = sxx/si; covar[1][2] = sxy/si; covar[1][3] = sxz/si;
    covar[2][1] = sxy/si; covar[2][2] = syy/si; covar[2][3] = syz/si;
    covar[3][1] = sxz/si; covar[3][2] = syz/si; covar[3][3] = szz/si;
    
    return(TRUE);
    
  }
  else {
    return(FALSE);
  }
}



char *prog_name;
main(int argc, char *argv[])
{
  Volume vol;
  float **cov, *cog;
  double step[3];
  Real x,y,z,r,s,c;

  if (argc!=2) {
    print ("usage: %s filename.mnc\n",argv[0]);
    exit(EXIT_FAILURE);
  }

  prog_name = argv[0];
  input_volume(argv[1],3,default_dim_names /*(char **)NULL*/, NC_UNSPECIFIED, FALSE, 0.0,0.0,
	       TRUE, &vol, (minc_input_options *)NULL);


  step[0] = 4.0;
  step[1] = 4.0;
  step[2] = 4.0;

  ALLOC2D(cov,4,4);
  ALLOC(cog,4);

  vol_to_cov(vol, NULL, cog, cov, step );

  (void)print ("%f %f %f\n",cog[1],cog[2],cog[3]);

  FREE2D(cov);
  FREE(cog);

}















