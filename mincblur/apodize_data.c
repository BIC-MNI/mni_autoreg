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

void apodize_data(VIO_Volume data, 
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
    sizes[3];
  VIO_Real
    step[3];

  get_volume_sizes(data, sizes);
  get_volume_separations(data, step);
  get_volume_voxel_range(data, &valid_min, &valid_max);

  if (zramp1 > ABS(step[VIO_Z])/2) {

    /* number of slices to be apodized */
    num_steps = ROUND( 1.25*zramp1/ABS(step[VIO_Z]) + 0.5);      

    if (num_steps>ABS(sizes[VIO_Z])) {
      print_error_and_line_num("FWHM is greater than slice dimension\n",__FILE__, __LINE__);
    }
    
    for(slice=0; slice<num_steps; slice++) {
      
      scale1 = normal_height( zramp1*GWID1, 1.25*zramp1, (float)slice*ABS(step[VIO_Z]));
      scale2 = normal_height( zramp1*GWID2, 0.0, (float)(num_steps - 1 - slice)*ABS(step[VIO_Z]));

      scale = INTERPOLATE( slice/(num_steps-1.0) , scale1, scale2);

      for(row=0; row<sizes[VIO_Y]; row++) {
        
        for(col=0; col<sizes[VIO_X]; col++) {

          GET_VOXEL_3D(voxel_value, data, col, row, slice);
          if (voxel_value >= valid_min && voxel_value <= valid_max) {

            real_value = CONVERT_VOXEL_TO_VALUE(data, voxel_value);

            real_value *= scale;

            voxel_value = CONVERT_VALUE_TO_VOXEL(data, real_value);

            SET_VOXEL_3D( data, col, row, slice, voxel_value);


          }

        }
        
      }
    }
    
  }
  
  if (zramp2 > ABS(step[VIO_Z])/2 ) {

    /* number of slices to be apodized */
    num_steps = ROUND( 1.25*zramp2/ABS(step[VIO_Z]) + 0.5);      

    if (num_steps>sizes[VIO_Z]) {
      print_error_and_line_num("FWHM is greater than slice dimension\n",__FILE__, __LINE__);
    }
    
    for(slice=0; slice<num_steps; slice++) {
      
      scale1 = normal_height( zramp2*GWID1, 1.25*zramp2, (float)slice*ABS(step[VIO_Z]));
      scale2 = normal_height( zramp2*GWID2, 0.0, (float)(num_steps - 1 - slice)*ABS(step[VIO_Z]));
      scale = INTERPOLATE( slice/(num_steps-1.0) , scale1, scale2);

      for(row=0; row<sizes[VIO_Y]; row++) {
        
        for(col=0; col<sizes[VIO_X]; col++) {

          GET_VOXEL_3D(voxel_value, data, col, row, sizes[VIO_Z]-1-slice);
          if (voxel_value >= valid_min && voxel_value <= valid_max) {
            real_value = CONVERT_VOXEL_TO_VALUE(data, voxel_value);
            real_value *= scale;
            voxel_value = CONVERT_VALUE_TO_VOXEL(data, real_value);
            SET_VOXEL_3D( data, col, row, sizes[VIO_Z]-1-slice, voxel_value);
          }
        }
        
      }
    }
    
  }
  
  if (yramp1 > ABS(step[VIO_Y])/2 ) {
                                /* number of rows to be apodized */
    num_steps = ROUND( 1.25*yramp1/ABS(step[VIO_Y]) + 0.5);      

    if (num_steps>sizes[VIO_Y]) {
      print_error_and_line_num("FWHM is greater than col dimension\n",__FILE__, __LINE__);
    }
    
    for(row=0; row<num_steps; row++) {
      
      scale1 = normal_height( yramp1*GWID1, 1.25*yramp1, (float)row*ABS(step[VIO_Y]));
      scale2 = normal_height( yramp1*GWID2, 0.0, (float)(num_steps - 1 - row)*ABS(step[VIO_Y]));
      scale = INTERPOLATE( row/(num_steps-1.0) , scale1, scale2);

      for(slice=0; slice<sizes[VIO_Z]; slice++) {
        
        for(col=0; col<sizes[VIO_X]; col++) {
          
          GET_VOXEL_3D(voxel_value, data, col, row, slice);
          if (voxel_value >= valid_min && voxel_value <= valid_max) {
            real_value = CONVERT_VOXEL_TO_VALUE(data, voxel_value);
            real_value *= scale;
            voxel_value = CONVERT_VALUE_TO_VOXEL(data, real_value);
            SET_VOXEL_3D( data, col, row, slice, voxel_value);
          }
        }
        
      }
    }
    
  }

  if (yramp2 > ABS(step[VIO_Y])/2) {
                                /* number of rows to be apodized */
    num_steps = ROUND( 1.25*yramp2/ABS(step[VIO_Y]) + 0.5);      
    
    if (num_steps>sizes[VIO_Y]) {
      print_error_and_line_num("FWHM is greater than col dimension\n",__FILE__, __LINE__);
    }
    
    for(row=0; row<num_steps; row++) {
      
      scale1 = normal_height( yramp2*GWID1, 1.25*yramp2, (float)row*ABS(step[VIO_Y]));
      scale2 = normal_height( yramp2*GWID2, 0.0, (float)(num_steps - 1 - row)*ABS(step[VIO_Y]));
      scale = INTERPOLATE( row/(num_steps-1.0) , scale1, scale2);
      
      for(slice=0; slice<sizes[VIO_Z]; slice++) {
        
        for(col=0; col<sizes[VIO_X]; col++) {
          
          GET_VOXEL_3D(voxel_value, data, col, sizes[VIO_Y]-1-row, slice);
          if (voxel_value >= valid_min && voxel_value <= valid_max) {
            real_value = CONVERT_VOXEL_TO_VALUE(data, voxel_value);
            real_value *= scale;
            voxel_value = CONVERT_VALUE_TO_VOXEL(data, real_value);
            SET_VOXEL_3D( data, col, sizes[VIO_Y]-1-row, slice, voxel_value);
          }
        }
        
      }
    }
    
  }
  
  if (xramp1 > ABS(step[VIO_X])/2) {
                                /* number of slices to be apodized */
    num_steps = ROUND( 1.25*xramp1/ABS(step[VIO_X]) + 0.5);      
    
    if (num_steps>sizes[VIO_X]) {
      print_error_and_line_num("FWHM is greater than col dimension\n",__FILE__, __LINE__);
    }
    
    for(col=0; col<num_steps; col++) {


      scale1 = normal_height( xramp1*GWID1, 1.25*xramp1, (float)col*ABS(step[VIO_X]));
      scale2 = normal_height( xramp1*GWID2, 0.0, (float)(num_steps - 1 - col)*ABS(step[VIO_X]));
      scale = INTERPOLATE( col/(num_steps-1.0) , scale1, scale2);

      for(slice=0; slice<sizes[VIO_Z]; slice++) {
        
        for(row=0; row<sizes[VIO_Y]; row++) {
          
          GET_VOXEL_3D( voxel_value, data, col, row, slice);
          if (voxel_value >= valid_min && voxel_value <= valid_max) {
            real_value = CONVERT_VOXEL_TO_VALUE(data, voxel_value);
            real_value *= scale;
            voxel_value = CONVERT_VALUE_TO_VOXEL(data, real_value);
            SET_VOXEL_3D( data, col, row, slice, voxel_value);
          }
          
        }
        
      }
    }
    
  }
  
  if (xramp2 > ABS(step[VIO_X])/2) {
                                /* number of slices to be apodized */
    num_steps = ROUND( 1.25*xramp2/step[VIO_X] + 0.5);      

    if (num_steps>sizes[VIO_X]) {
      print_error_and_line_num("FWHM is greater than col dimension\n",__FILE__, __LINE__);
    }
    
    for(col=0; col<num_steps; col++) {
      
      scale1 = normal_height( xramp2*GWID1, 1.25*xramp2, (float)col*ABS(step[VIO_X]));
      scale2 = normal_height( xramp2*GWID2, 0.0, (float)(num_steps - 1 - col)*ABS(step[VIO_X]));
      scale = INTERPOLATE( col/(num_steps-1.0) , scale1, scale2);

      for(slice=0; slice<sizes[VIO_Z]; slice++) {
        
        for(row=0; row<sizes[VIO_Y]; row++) {
          
          GET_VOXEL_3D( voxel_value, data, sizes[VIO_X]-1-col, row, slice);
          if (voxel_value >= valid_min && voxel_value <= valid_max) {
            real_value = CONVERT_VOXEL_TO_VALUE(data, voxel_value);
            real_value *= scale;
            voxel_value = CONVERT_VALUE_TO_VOXEL(data, real_value);
            SET_VOXEL_3D( data, sizes[VIO_X]-1-col, row, slice, voxel_value);
          }
          
        }
        
      }
    }
    
  }
  
  
  
  
}

