#include <volume_io.h>




#include <volume_io.h>
#include <recipes.h>

#include "constants.h"
#include "matrix_basics.h"
#include "cov_to_praxes.h"
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


char *prog_name;


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



BOOLEAN get_cog(char *file, double *c1)
{
  Volume vol;
  float **cov, *cog;
  double step[3];
  Real x,y,z,r,s,c;

  input_volume(file,3,default_dim_names /*(char **)NULL*/, NC_UNSPECIFIED, FALSE, 0.0,0.0,
	       TRUE, &vol, (minc_input_options *)NULL);

  step[0] = 4.0;
  step[1] = 4.0;
  step[2] = 4.0;

  cov = matrix(1,3,1,3); 
  cog = vector(1,3);

  if ( vol_to_cov(vol, NULL, cog, cov, step ) ) {
    c1[0] = cog[1];
    c1[1] = cog[1];
    c1[2] = cog[3];
    return(TRUE);
  }
  else
    return(FALSE);


}


main(int argc, char *argv[])
{
  float
    c1[4];
  double 
    cent[3],
    rots[3],
    tran[3],
    scal[3],
    sher[6];
  Transform 
    *mat;
  General_transform
    gt;
  char *xfmfile;
  int i;
  Volume data;
  
  if (argc<2) {
    print ("usage:  xfm_to_param  file.xfm  [file.mnc]\n");
    exit(EXIT_FAILURE);
  }

  prog_name = argv[0];
  xfmfile   = argv[1];

  for_less(i,0,3) {
    cent[i] = rots[i] = tran[i] = scal[i] = 0.0;
  }
  for_less(i,0,6) sher[i] = 0.0;

  if (argc>2) {
    print ("mnc = %s\n",argv[2]);
    if (! get_cog(argv[2], cent) ) {
      print("Cannot calculate the COG of volume %s\n.", argv[2] );
      return(FALSE);
    }
  }

  if (input_transform_file(xfmfile, &gt)!=OK) {
    (void)fprintf(stderr, "Error reading transformation file.\n");
    exit(EXIT_FAILURE);
  }
  
  mat = get_linear_transform_ptr(&gt);

  extract2_parameters_from_matrix(mat,    cent, tran, scal, sher, rots);

  print("after parameter extraction\n");
  print("-center      %10.5f %10.5f %10.5f\n", cent[0], cent[1], cent[2]);
  print("-translation %10.5f %10.5f %10.5f\n", tran[0], tran[1], tran[2]);
  print("-rotation    %10.5f %10.5f %10.5f\n", 
	rots[0]*180.0/3.1415927, rots[1]*180.0/3.1415927, rots[2]*180.0/3.1415927);
  print("-scale       %10.5f %10.5f %10.5f\n", scal[0], scal[1], scal[2]);
  print("-shear       %10.5f %10.5f %10.5f\n", sher[0], sher[1], sher[2]);
  
  
}
