/* ----------------------------- MNI Header -----------------------------------
@NAME       : apodize data.c
@INPUT      : volume of data + depth of apodization on each surface.
@CREATED    : Mon Jul  5 13:35:04 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */

#include <config.h>
#include <volume_io.h>
#include <print_error.h>


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

void apodize_data(Volume data, 
			 double xramp1,double xramp2,
			 double yramp1,double yramp2,
			 double zramp1,double zramp2)
{
  int
    num_steps,
    row,col,slice;
  float
    scale1,scale2, scale;
  
  Real
    valid_min, valid_max,
    real_value,
    voxel_value;

  int
    sizes[3];
  Real
    step[3];

  get_volume_sizes(data, sizes);
  get_volume_separations(data, step);
  get_volume_voxel_range(data, &valid_min, &valid_max);

  if (zramp1 > ABS(step[Z])/2) {

    /* number of slices to be apodized */
    num_steps = ROUND( 1.25*zramp1/ABS(step[Z]) + 0.5);      

    if (num_steps>ABS(sizes[Z])) {
      print_error_and_line_num("FWHM is greater than slice dimension\n",__FILE__, __LINE__);
    }
    
    for_less(slice, 0, num_steps) {
      
      scale1 = normal_height( zramp1*GWID1, 1.25*zramp1, (float)slice*ABS(step[Z]));
      scale2 = normal_height( zramp1*GWID2, 0.0, (float)(num_steps - 1 - slice)*ABS(step[Z]));

      scale = INTERPOLATE( slice/(num_steps-1.0) , scale1, scale2);

      for_less(row, 0, sizes[Y]) {
	
	for_less(col, 0, sizes[X]) {

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
  
  if (zramp2 > ABS(step[Z])/2 ) {

    /* number of slices to be apodized */
    num_steps = ROUND( 1.25*zramp2/ABS(step[Z]) + 0.5);      

    if (num_steps>sizes[Z]) {
      print_error_and_line_num("FWHM is greater than slice dimension\n",__FILE__, __LINE__);
    }
    
    for_less(slice, 0, num_steps) {
      
      scale1 = normal_height( zramp2*GWID1, 1.25*zramp2, (float)slice*ABS(step[Z]));
      scale2 = normal_height( zramp2*GWID2, 0.0, (float)(num_steps - 1 - slice)*ABS(step[Z]));
      scale = INTERPOLATE( slice/(num_steps-1.0) , scale1, scale2);

      for_less(row, 0, sizes[Y]) {
	
	for_less(col, 0, sizes[X]) {

	  GET_VOXEL_3D(voxel_value, data, col, row, sizes[Z]-1-slice);
	  if (voxel_value >= valid_min && voxel_value <= valid_max) {
	    real_value = CONVERT_VOXEL_TO_VALUE(data, voxel_value);
	    real_value *= scale;
	    voxel_value = CONVERT_VALUE_TO_VOXEL(data, real_value);
	    SET_VOXEL_3D( data, col, row, sizes[Z]-1-slice, voxel_value);
	  }
	}
	
      }
    }
    
  }
  
  if (yramp1 > ABS(step[Y])/2 ) {
				/* number of rows to be apodized */
    num_steps = ROUND( 1.25*yramp1/ABS(step[Y]) + 0.5);      

    if (num_steps>sizes[Y]) {
      print_error_and_line_num("FWHM is greater than col dimension\n",__FILE__, __LINE__);
    }
    
    for_less(row, 0, num_steps) {
      
      scale1 = normal_height( yramp1*GWID1, 1.25*yramp1, (float)row*ABS(step[Y]));
      scale2 = normal_height( yramp1*GWID2, 0.0, (float)(num_steps - 1 - row)*ABS(step[Y]));
      scale = INTERPOLATE( row/(num_steps-1.0) , scale1, scale2);

      for_less(slice, 0, sizes[Z]) {
	
	for_less(col, 0, sizes[X]) {
	  
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

  if (yramp2 > ABS(step[Y])/2) {
				/* number of rows to be apodized */
    num_steps = ROUND( 1.25*yramp2/ABS(step[Y]) + 0.5);      
    
    if (num_steps>sizes[Y]) {
      print_error_and_line_num("FWHM is greater than col dimension\n",__FILE__, __LINE__);
    }
    
    for_less(row, 0, num_steps) {
      
      scale1 = normal_height( yramp2*GWID1, 1.25*yramp2, (float)row*ABS(step[Y]));
      scale2 = normal_height( yramp2*GWID2, 0.0, (float)(num_steps - 1 - row)*ABS(step[Y]));
      scale = INTERPOLATE( row/(num_steps-1.0) , scale1, scale2);
      
      for_less(slice, 0, sizes[Z]) {
	
	for_less(col, 0, sizes[X]) {
	  
	  GET_VOXEL_3D(voxel_value, data, col, sizes[Y]-1-row, slice);
	  if (voxel_value >= valid_min && voxel_value <= valid_max) {
	    real_value = CONVERT_VOXEL_TO_VALUE(data, voxel_value);
	    real_value *= scale;
	    voxel_value = CONVERT_VALUE_TO_VOXEL(data, real_value);
	    SET_VOXEL_3D( data, col, sizes[Y]-1-row, slice, voxel_value);
	  }
	}
	
      }
    }
    
  }
  
  if (xramp1 > ABS(step[X])/2) {
				/* number of slices to be apodized */
    num_steps = ROUND( 1.25*xramp1/ABS(step[X]) + 0.5);      
    
    if (num_steps>sizes[X]) {
      print_error_and_line_num("FWHM is greater than col dimension\n",__FILE__, __LINE__);
    }
    
    for_less(col, 0, num_steps) {


      scale1 = normal_height( xramp1*GWID1, 1.25*xramp1, (float)col*ABS(step[X]));
      scale2 = normal_height( xramp1*GWID2, 0.0, (float)(num_steps - 1 - col)*ABS(step[X]));
      scale = INTERPOLATE( col/(num_steps-1.0) , scale1, scale2);

      for_less(slice, 0, sizes[Z]) {
	
	for_less(row, 0, sizes[Y]) {
	  
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
  
  if (xramp2 > ABS(step[X])/2) {
				/* number of slices to be apodized */
    num_steps = ROUND( 1.25*xramp2/step[X] + 0.5);      

    if (num_steps>sizes[X]) {
      print_error_and_line_num("FWHM is greater than col dimension\n",__FILE__, __LINE__);
    }
    
    for_less(col, 0, num_steps) {
      
      scale1 = normal_height( xramp2*GWID1, 1.25*xramp2, (float)col*ABS(step[X]));
      scale2 = normal_height( xramp2*GWID2, 0.0, (float)(num_steps - 1 - col)*ABS(step[X]));
      scale = INTERPOLATE( col/(num_steps-1.0) , scale1, scale2);

      for_less(slice, 0, sizes[Z]) {
	
	for_less(row, 0, sizes[Y]) {
	  
	  GET_VOXEL_3D( voxel_value, data, sizes[X]-1-col, row, slice);
	  if (voxel_value >= valid_min && voxel_value <= valid_max) {
	    real_value = CONVERT_VOXEL_TO_VALUE(data, voxel_value);
	    real_value *= scale;
	    voxel_value = CONVERT_VALUE_TO_VOXEL(data, real_value);
	    SET_VOXEL_3D( data, sizes[X]-1-col, row, slice, voxel_value);
	  }
	  
	}
	
      }
    }
    
  }
  
  
  
  
}

