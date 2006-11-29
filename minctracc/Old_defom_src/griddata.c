#include <def_standard.h>
#include <ms_iffheader.h>

make_grid(data)

     DATA *data;

{
   int
      i,j,k;
   char
      *c1, *c2, ct;



                                /* make lines on every 16th row */
   for (i=0; i < data->slices; ++i)
      for (j = 8; j < data->rows; j+=16) {
         
         c1 = data->voxels + i*data->slice_size + j*data->cols*data->bytes_per_voxel;
         
         for (k=0; k<data->cols; k++) {
            *c1 = 255;
            c1++;
         }
      }
   
   
   for (i=0; i < data->slices; ++i)
      for (j = 0; j < data->rows; ++j) {
         for (k=8; k<data->cols; k+=16) {
            
            c1 = data->voxels + i*data->slice_size + (j*data->cols + k)*data->bytes_per_voxel;
            *c1 = 255;
            
         }
      }
   
   for (i=0; i < data->slices; ++i) 
      for (j = 8; j < data->rows; j+=16) {
         for (k=8; k<data->cols; k+=16) {
            
            c1 = data->voxels + i*data->slice_size + (j*data->cols + k)*data->bytes_per_voxel;
            *c1 = 255;
            
         }
      }
   

      


}

