
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <volume_io.h>
#include <config.h>
#include <Proglib.h>
#include <ParseArgv.h>

#include "constants.h"
#include "matrix_basics.h"
#include "make_rots.h"
#include "point_vector.h"
#include "interpolation.h"



#define DO_TRANSFORM(result, transformation, coord) \
   general_transform_point(transformation, \
      Point_x(coord), Point_y(coord), Point_z(coord), \
      &(Point_x(result)), &(Point_y(result)), &(Point_z(result)) )

#define INTERPOLATE_TRUE_VALUE(volume, coord, result) \
   trilinear_interpolant(volume, coord, result)

static char *default_dim_names[VIO_N_DIMENSIONS] =
   { MIzspace, MIyspace, MIxspace };


char *prog_name;


VIO_BOOL vol_to_cov(VIO_Volume d1, VIO_Volume m1, float centroid[4], float covar[4][4], double *step)
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

  VIO_Real
    tx,ty,tz;
  int
    i,r,c,s,
    limits[VOL_NDIMS];

  float
    t,
    sxx,syy,szz,
    sxy,syz,sxz,
    sx,sy,sz,si; 
  VIO_Real
    thickness[3];
  int
    sizes[3];

  VIO_Real
    true_value;

  get_volume_separations(d1, thickness);
  get_volume_sizes(d1, sizes);
  
                                /* build sampling lattice info */
  for(i=0; i<3; i++) {        
    step[i] *= thickness[i] / fabs( thickness[i]);
  }

  fill_Vector( col_step,   step[COL_IND], 0.0,     0.0 );
  fill_Vector( row_step,   0.0,     step[ROW_IND], 0.0 );
  fill_Vector( slice_step, 0.0,     0.0,     step[SLICE_IND] );

  convert_3D_voxel_to_world(d1, 0.0, 0.0, 0.0, &tx, &ty, &tz); 

  fill_Point( starting_origin, tx, ty, tz);

  for(i=0; i<3; i++) {                /* for each dim, get # of steps in that direction,
                                   and set starting offset */
    t = sizes[i] * thickness[i] / step[i];
    limits[i] = abs( t );
    
    Point_coord( starting_offset, (i) ) = 
      ( (sizes[i]-1)*thickness[i] - (limits[i] * step[i] ) ) / 2.0;
  }
  
  ADD_POINTS( starting_position, starting_origin, starting_offset ); /*  */

                                /* calculate centroids first */

  sx = 0.0;
  sy = 0.0;
  sz = 0.0;
  si = 0.0;

  for(s=0; s<=limits[SLICE_IND]; s++) {

    SCALE_VECTOR( vector_step, slice_step, s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    for(r=0; r<=limits[ROW_IND]; r++) {

      SCALE_VECTOR( vector_step, row_step, r);
      ADD_POINT_VECTOR( row, slice, vector_step );

      SCALE_POINT( col, row, 1.0); /* init first col position */
      for(c=0; c<=limits[COL_IND]; c++) {

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

    for(s=0; s<=limits[VIO_Z]; s++) {
      
      SCALE_VECTOR( vector_step, slice_step, s);
      ADD_POINT_VECTOR( slice, starting_position, vector_step );
      
      for(r=0; r<=limits[VIO_Y]; r++) {
        
        SCALE_VECTOR( vector_step, row_step, r);
        ADD_POINT_VECTOR( row, slice, vector_step );
        
        SCALE_POINT( col, row, 1.0); /* init first col position */
        for(c=0; c<=limits[VIO_X]; c++) {
          
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



VIO_BOOL get_cog(char *file, double *c1)
{
  VIO_Volume vol;
  float cov[4][4], cog[4];
  double step[3];

  input_volume(file,3,default_dim_names /*(char **)NULL*/, NC_UNSPECIFIED, FALSE, 0.0,0.0,
               TRUE, &vol, (minc_input_options *)NULL);

  step[0] = 4.0;
  step[1] = 4.0;
  step[2] = 4.0;

  if ( vol_to_cov(vol, NULL, cog, cov, step ) ) {
    c1[0] = cog[1];
    c1[1] = cog[2];
    c1[2] = cog[3];
    return(TRUE);
  }
  else
    return(FALSE);


}

/* Centre of rotation and scale */
double cent[3];


void process_linear_transform( Transform* mat )
{
  int i;
  double 
    rots[3],
    tran[3],
    scal[3],
    sher[6];

  for(i=0; i<3; i++) {
    rots[i] = tran[i] = scal[i] = 0.0;
  }
  for(i=0; i<6; i++) sher[i] = 0.0;

  extract2_parameters_from_matrix(mat,    cent, tran, scal, sher, rots);

  print("after parameter extraction\n");
  print("-center      %10.5f %10.5f %10.5f\n", cent[0], cent[1], cent[2]);
  print("-translation %10.5f %10.5f %10.5f\n", tran[0], tran[1], tran[2]);
  print("-rotation    %10.5f %10.5f %10.5f\n", 
        rots[0]*180.0/3.1415927, rots[1]*180.0/3.1415927, rots[2]*180.0/3.1415927);
  print("-scale       %10.5f %10.5f %10.5f\n", scal[0], scal[1], scal[2]);
  print("-shear       %10.5f %10.5f %10.5f\n", sher[0], sher[1], sher[2]);
  
}


void process_general_transform( General_transform* t);

void process_concatenated_transform( General_transform* t )
{
    int n = get_n_concated_transforms( t );
    int i;

    print("[CONCATENATED TRANSFORM]\n");
    for( i = 0; i < n; ++i ) {
        print("Transform %d: ", i+1);
        process_general_transform( get_nth_general_transform(t,i) );
        print("\n\n");
    }
}
            

void process_general_transform( General_transform* t )
{
    switch (t->type) {
    case LINEAR:
        process_linear_transform( get_linear_transform_ptr(t) );
        break;
    case THIN_PLATE_SPLINE:
        print("[THIN PLATE SPLINE]\n");
        break;
    case USER_TRANSFORM:
        print("[USER SPECIFIED TRANSFORM]\n");
        break;
    case CONCATENATED_TRANSFORM:
        process_concatenated_transform( t );
        break;
    case GRID_TRANSFORM:
        print("[GRID TRANSFORM]\n");
        break;
    }
}


int main(int argc, char *argv[])
{
  VIO_General_transform gt;
  char *xfmfile;
  int i;

    
  static ArgvInfo argTable[] = {
    {"-center",      ARGV_FLOAT, (char *) 3, (char *)cent,
         "Force center of rotation and scale."},
    {"-version", ARGV_FUNC, (char *) print_version_info, (char *)MNI_AUTOREG_LONG_VERSION,
         "Print out version info and exit."},
    {NULL, ARGV_END, NULL, NULL, NULL}
  };
   
  for(i=0; i<3; i++) {
    cent[i] = 0.0;
  }

  if (ParseArgv(&argc, argv, argTable, 0) || (argc<2)) {
    (void) fprintf(stderr, "Usage: %s file.xfm [file.mnc] [options] \n",
                   argv[0]);
    exit(EXIT_FAILURE);
  }
  
  prog_name = argv[0];
  xfmfile   = argv[1];

  if (argc>2) {
    print ("mnc = %s\n",argv[2]);
    if (! get_cog(argv[2], cent) ) {
      print("Cannot calculate the COG of volume %s\n.", argv[2] );
      exit(EXIT_FAILURE);
    }
  }

  if (input_transform_file(xfmfile, &gt)!=OK) {
    (void)fprintf(stderr, "Error reading transformation file.\n");
    exit(EXIT_FAILURE);
  }

  process_general_transform( &gt );
  exit(EXIT_SUCCESS);
}


