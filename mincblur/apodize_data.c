#include <standard.h>
#include <ms_iffheader.h>

/************************************************************/
/* return the value of the normal dist at x, given c,sigma  */
/* and mu ----   all in mm                                  */
/************************************************************/
private float normal_dist(c,fwhm,mu,x)
float c,fwhm,mu,x;
{
   float sigma,t1,t2,t3,f;

   sigma = fwhm/2;

   if (sigma==0) {
      if (x==mu)
	 f = c;
      else
	 f = 0;
   }
   else {
      t1 = c / (sqrt(2*PI) * sigma);
      t2 = (x-mu)*(x-mu)/(2*sigma*sigma);
      t3 = exp(-t2);
      f = t1*t3;
   }


   return(f);
}

apodize_data(data, xramp1,xramp2,yramp1,yramp2,zramp1,zramp2)

DATA 
   *data;
float 
   xramp1,xramp2,yramp1,yramp2,zramp1,zramp2;
     
{

   int
      num_steps,
      i,j,k;
   char
      *c1, *c2, ct;
   float
      scale;

   if (zramp1 > 0.0001) {
				/* number of slices to be apodized */
      num_steps = ROUND( zramp1/data->slice_width + 0.5);      


      for (i=0; i < num_steps; ++i) {

	 scale = normal_dist(sqrt(2*PI)*0.5*zramp1*0.75,zramp1*0.75,zramp1,(float)i*data->slice_width);

	 for (j = 0; j < data->rows; ++j) {

	    c1 = data->voxels + i*data->slice_size + j*data->cols*data->bytes_per_voxel;

	    for (k=0; k<data->cols; k++) {
	       *c1 = *c1 * scale;
	       c1++;
	    }

	 }
      }

   }

   if (zramp2 > 0.0001) {
				/* number of slices to be apodized */
      num_steps = ROUND( zramp2/data->slice_width + 0.5);      

      for (i=0; i < num_steps; ++i) {

	 scale = normal_dist(sqrt(2*PI)*0.5*zramp2*0.75,zramp2*0.75,zramp2,(float)i*data->slice_width);

	 for (j = 0; j < data->rows; ++j) {

	    c1 = data->voxels + (data->slices-1-i)*data->slice_size + j*data->cols*data->bytes_per_voxel;

	    for (k=0; k<data->cols; k++) {
	       *c1 = *c1 * scale;
	       c1++;
	    }

	 }
      }

   }

   if (yramp1 > 0.0001) {
				/* number of rows to be apodized */
      num_steps = ROUND( yramp1/data->pixel_size_row + 0.5);      

      for (i=0; i < data->slices; ++i) {

	 for (j = 0; j < num_steps; ++j) {

	 scale = normal_dist(sqrt(2*PI)*0.5*yramp1*0.75,yramp1*0.75,yramp1,(float)j*data->pixel_size_row);

	    c1 = data->voxels + i*data->slice_size + j*data->cols*data->bytes_per_voxel;

	    for (k=0; k<data->cols; k++) {
	       *c1 = *c1 * scale;
	       c1++;
	    }

	 }
      }

   }

   if (yramp2 > 0.0001) {
				/* number of rows to be apodized */
      num_steps = ROUND( yramp2/data->pixel_size_row + 0.5);      

      for (i=0; i < data->slices; ++i) {

	 for (j = 0; j < num_steps; ++j) {

	 scale = normal_dist(sqrt(2*PI)*0.5*yramp2*0.75,yramp2*0.75,yramp2,(float)j*data->pixel_size_row);

	    c1 = data->voxels + i*data->slice_size + (data->rows-1-j)*data->cols*data->bytes_per_voxel;

	    for (k=0; k<data->cols; k++) {
	       *c1 = *c1 * scale;
	       c1++;
	    }

	 }
      }

   }

   if (xramp1 > 0.0001) {
				/* number of slices to be apodized */
      num_steps = ROUND( xramp1/data->pixel_size_col + 0.5);      

      for (i=0; i < data->slices; ++i) {

	 for (j = 0; j < data->rows; ++j) {
	 
	    c1 = data->voxels + i*data->slice_size + j*data->cols*data->bytes_per_voxel;

	    for (k=0; k<num_steps; k++) {

	       scale = normal_dist(sqrt(2*PI)*0.5*xramp1*0.75,xramp1*0.75,xramp1,
				   (float)k*data->pixel_size_col);

	       *c1 = *c1 * scale;
	       c1++;
	    }

	 }
      }

   }


   if (xramp2 > 0.0001) {
				/* number of slices to be apodized */
      num_steps = ROUND( xramp2/data->pixel_size_col + 0.5);      

      for (i=0; i < data->slices; ++i) {

	 for (j = 0; j < data->rows; ++j) {
	 
	    c1 = data->voxels + i*data->slice_size + ((j+1)*data->cols-1)*data->bytes_per_voxel;

	    for (k=0; k<num_steps; k++) {

	       scale = normal_dist(sqrt(2*PI)*0.5*xramp2*0.75,xramp2*0.75,xramp2,
				   (float)k*data->pixel_size_col);

	       *c1 = *c1 * scale;
	       c1--;
	    }

	 }
      }

   }



}

