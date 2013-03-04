/* ----------------------------- MNI Header -----------------------------------
@NAME       : apodize data.c
@INPUT      : volume of data + depth of apodization on each surface.
@CREATED    : Mon Jul  5 13:35:04 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */

#include <config.h>
#include <volume_io.h>
#include <Proglib.h>

#define  FLOOR( x )     ((int) floor(x))
#define  ROUND( x )     FLOOR( (double) (x) + 0.5 )
#define  ABS( x )   ( ((x) > 0) ? (x) : (-(x)) )
#define  INTERPOLATE( alpha, a, b ) ((a) + (alpha) * ((b) - (a)))

#define GWID1 0.75
#define GWID2 0.9

/************************************************************/
/* return the value of the normal distribution (that has a  */
/* normalized height of 1.0,  given sigma, x and mu         */
/*  all in mm                                               */
/************************************************************/
static float normal_height(float fwhm,float mu,float x)
{
   float sigma,t2,f;

   sigma = fwhm/2.35482;

   if (sigma==0) {
      if (x==mu)
         f = 1.0;
      else
         f = 0;
   }
   else {
      t2 = (x-mu)*(x-mu)/(2*sigma*sigma);
      f = exp(-t2);
   }


   return(f);
}

void apodize_data(VIO_Volume data, int xyzv[VIO_MAX_DIMENSIONS],
                         double xramp1,double xramp2,
                         double yramp1,double yramp2,
                         double zramp1,double zramp2)
{
  int
    num_steps,
    row,col,slice;
  float
    scale1,scale2, scale;
  
  VIO_Real
    valid_min, valid_max,
    real_value,
    voxel_value;

  int
    sizes[3], pos[3];
  VIO_Real
    step[3], ramp1[3], ramp2[3];

  get_volume_sizes(data, sizes);
  get_volume_separations(data, step);
  get_volume_voxel_range(data, &valid_min, &valid_max);

  ramp1[VIO_X] = xramp1;
  ramp1[VIO_Y] = yramp1;
  ramp1[VIO_Z] = zramp1;
  ramp2[VIO_X] = xramp2;
  ramp2[VIO_Y] = yramp2;
  ramp2[VIO_Z] = zramp2;

  int dim0, dim1, dim2;

  for( dim0 = 0; dim0 < 3; dim0++ ) {

    dim1 = ( dim0 + 1 ) % 3;
    dim2 = ( dim0 + 2 ) % 3;

    if (ramp1[dim0] > ABS(step[xyzv[dim0]])/2) {

      /* number of slices to be apodized */
      num_steps = ROUND( 1.25*ramp1[dim0]/ABS(step[xyzv[dim0]]) + 0.5);      

      if (num_steps>ABS(sizes[xyzv[dim0]])) {
        print_error_and_line_num("FWHM is greater than slice dimension\n",__FILE__, __LINE__);
      }
    
      for(slice=0; slice<num_steps; slice++) {
      
        scale1 = normal_height( ramp1[dim0]*GWID1, 1.25*ramp1[dim0], (float)slice*ABS(step[xyzv[dim0]]));
        scale2 = normal_height( ramp1[dim0]*GWID2, 0.0, (float)(num_steps - 1 - slice)*ABS(step[xyzv[dim0]]));
        scale = INTERPOLATE( slice/(num_steps-1.0) , scale1, scale2);

        for(row=0; row<sizes[xyzv[dim1]]; row++) {
          for(col=0; col<sizes[xyzv[dim2]]; col++) {
            pos[xyzv[dim0]] = slice;
            pos[xyzv[dim1]] = row;
            pos[xyzv[dim2]] = col;
            GET_VOXEL_3D(voxel_value, data, pos[0], pos[1], pos[2] );
            if (voxel_value >= valid_min && voxel_value <= valid_max) {
              real_value = CONVERT_VOXEL_TO_VALUE(data, voxel_value);
              real_value *= scale;
              voxel_value = CONVERT_VALUE_TO_VOXEL(data, real_value);
              SET_VOXEL_3D( data, pos[0], pos[1], pos[2], voxel_value);
            }
          }
        }
      }
    }
    
    if (ramp2[dim0] > ABS(step[xyzv[dim0]])/2 ) {

      /* number of slices to be apodized */
      num_steps = ROUND( 1.25*ramp2[dim0]/ABS(step[xyzv[dim0]]) + 0.5);      

      if (num_steps>sizes[xyzv[dim0]]) {
        print_error_and_line_num("FWHM is greater than slice dimension\n",__FILE__, __LINE__);
      }
    
      for(slice=0; slice<num_steps; slice++) {
      
        scale1 = normal_height( ramp2[dim0]*GWID1, 1.25*ramp2[dim0], (float)slice*ABS(step[xyzv[dim0]]));
        scale2 = normal_height( ramp2[dim0]*GWID2, 0.0, (float)(num_steps - 1 - slice)*ABS(step[xyzv[dim0]]));
        scale = INTERPOLATE( slice/(num_steps-1.0) , scale1, scale2);

        for(row=0; row<sizes[xyzv[dim1]]; row++) {
          for(col=0; col<sizes[xyzv[dim2]]; col++) {
            pos[xyzv[dim0]] = sizes[xyzv[dim0]]-1-slice;
            pos[xyzv[dim1]] = row;
            pos[xyzv[dim2]] = col;
            GET_VOXEL_3D(voxel_value, data, pos[0], pos[1], pos[2] );
            if (voxel_value >= valid_min && voxel_value <= valid_max) {
              real_value = CONVERT_VOXEL_TO_VALUE(data, voxel_value);
              real_value *= scale;
              voxel_value = CONVERT_VALUE_TO_VOXEL(data, real_value);
              SET_VOXEL_3D( data, pos[0], pos[1], pos[2], voxel_value );
            }
          }
        }
      }
    }
  }
  
}

