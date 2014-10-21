/* ----------------------------- MNI Header -----------------------------------
@NAME       : blur_volume.c
@INPUT      : data - a pointer to a volume_struct of data
              fwhm - full-width-half-maximum of the gaussian blurring kernel
              outfile - name of the base filename to store the <name>_blur.mnc
              ndim - =1, do blurring in the z direction only,
                     =2, do blurring in the x and y directions only,
                     =3, blur in all three directions.
@OUTPUT     : creates and stores the blurred volume in an output file.
@RETURNS    : status variable - VIO_OK or ERROR.
@DESCRIPTION: This routine convolves each row, column and slice with a
              gaussian kernel.  The convolution is accomplished in the fourier
              domain by multiplying the fourier transformations of both the
              data and the kernel.
@METHOD     : 
@GLOBALS    : 
@CALLS      : stuff from volume_support.c and libmni.a
@CREATED    : Wed Jun 23 09:04:34 EST 1993  Louis Collins 
                 from code originally written for blur_grad working on .iff files.
@COPYRIGHT  :
              Copyright 1995 Louis Collins, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.

@MODIFIED   : $Log: blur_volume.c,v $
@MODIFIED   : Revision 96.5  2009-07-23 22:34:00  claude
@MODIFIED   : cleanup in mincblur for VIO_X, VIO_Y, VIO_Z
@MODIFIED   :
@MODIFIED   : Revision 96.4  2006/11/28 09:12:21  rotor
@MODIFIED   :  * fixes to allow clean compile against minc 2.0
@MODIFIED   :
@MODIFIED   : Revision 96.3  2005/07/20 20:45:39  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 96.2  2004/02/12 05:53:48  rotor
@MODIFIED   :  * removed /static defs
@MODIFIED   :
@MODIFIED   : Revision 96.1  2000/01/27 16:51:27  louis
@MODIFIED   : working version
@MODIFIED   :
 * Revision 96.0  1996/08/21  18:22:24  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:22:17  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:16:28  louis
 * Never released version 0.95
 *
 * Revision 1.13  1996/08/12  14:16:19  louis
 * Pre-release
 *
 * Revision 1.12  1995/09/18  09:02:42  collins
 * new functional version of mincblur with new and improved default behavior.
 * By default, only the blurred volume is created. If you want the gradient
 * magnitude volumes, you have to ask for it (-gradient).
 *
 * The temporary files (corresponding to the REAL data, and partial derivatives)
 * are removed by mincblur, so now it runs cleanly.  Unfortunately, these files
 * are _not_ deleted if mincblur fails, or is stopped.  I will have to use
 * unlink for this, but its a bit to much to do right now, since I would have
 * to change the way files are dealt with in gradient_volume.c, blur_volume.c
 * and gradmag_volume.c.
 *
 * Also, it is possible to keep the partial derivative volumes (-partial).
 *
 * this version is in mni_reg-0.1j
 *
 * Revision 1.11  1995/09/18  06:45:42  collins
 * this file is a working version of mincblur.  All references to numerical
 * recipes routines have been removed.  This version is included in the
 * package mni_reg-0.1i
 *
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/mincblur/blur_volume.c,v 96.5 2009-07-23 22:34:00 claude Exp $";
#endif

#include <float.h>
#include <volume_io.h>
#include <config.h>
#include <Proglib.h>
#include "blur_support.h"

extern int debug;

int ms_volume_reals_flag;

void fft1(float *signal, int numpoints, int direction);

VIO_Status blur3D_volume(VIO_Volume data, int xyzv[VIO_MAX_DIMENSIONS],
                            double fwhmx, double fwhmy, double fwhmz, 
                            char *infile,
                            char *outfile, 
                            FILE *reals_fp,
                            int kernel_type, char *history)
{ 
  float 
    *fdata,                        /* floating point storage for blurred volume */
    *f_ptr,                        /* pointer to fdata */
    tmp,
    *dat_vector,                /* temp storage of original row, col or slice vect. */
    *dat_vecto2,                /* storage of result of dat_vector*kern             */
    *kern;                        /* convolution kernel                               */
                                /*  place it back into data->voxels                 */

  VIO_Real
    lowest_val,
    max_val, 
    min_val;
    
  int                                
    total_voxels,                
    vector_size_data,                /* original size of row, col or slice vector        */
    kernel_size_data,                /* original size of kernel vector                   */
    array_size_pow2,                /* actual size of vector/kernel data used in FFT    */
                                /* routines - needs to be a power of two            */
    array_size;
  int   
    data_offset;                /* offset required to place original data (size n)  */
                                /*  into array (size m=2^n) so that data is centered*/

  
  register int 
    slice_limit,
    row,col,slice,                /* counters to access original data                 */
    vindex;                        /* counter to access vector and vecto2              */

  int 
    slice_size,                        /* size of each data step - in bytes              */
    row_size, col_size;
                
  char
    full_outfilename[1024];        /* name of output file */

  VIO_progress_struct 
    progress;                        /* used to monitor progress of calculations         */

  VIO_Status 
    status;
  
  int
    sizes[3],                        /* number of rows, cols and slices */
    pos[3];                          /* Input order of rows, cols, slices */
  VIO_Real
    steps[3];                        /* size of voxel step from center to center in x,y,z */


  /*---------------------------------------------------------------------------------*/
  /*             start by setting up the raw data.                                   */
  /*---------------------------------------------------------------------------------*/

  get_volume_sizes(data, sizes);          /* rows,cols,slices */
  get_volume_separations(data, steps);
  
  slice_size = sizes[xyzv[VIO_X]] * sizes[xyzv[VIO_Y]];    /* sizeof one slice  */
  col_size   = sizes[xyzv[VIO_Y]];                         /* sizeof one column */
  row_size   = sizes[xyzv[VIO_X]];                         /* sizeof one row    */

  total_voxels = sizes[xyzv[VIO_X]]*sizes[xyzv[VIO_Y]]*sizes[xyzv[VIO_Z]];

  ALLOC(fdata, total_voxels);

  lowest_val = get_volume_voxel_min(data);
    

  get_volume_real_range(data, &min_val, &max_val);

  if (debug) 
    print("Volume def min and max: = %f %f\n", min_val, max_val);

  max_val = -FLT_MAX;
  min_val = FLT_MAX;

  f_ptr = fdata;
  for(slice=0; slice<sizes[xyzv[VIO_Z]]; slice++) {
    pos[xyzv[VIO_Z]] = slice;
    for(row=0; row<sizes[xyzv[VIO_Y]]; row++) {
      pos[xyzv[VIO_Y]] = row;
      for(col=0; col<sizes[xyzv[VIO_X]]; col++) {
        pos[xyzv[VIO_X]] = col;

        GET_VOXEL_3D( tmp, data, pos[0], pos[1], pos[2] );

        if (tmp <= lowest_val)
          tmp = lowest_val;

        *f_ptr = CONVERT_VOXEL_TO_VALUE(data, tmp);
        if (max_val < *f_ptr) max_val = *f_ptr;
        if (min_val > *f_ptr) min_val = *f_ptr;
        f_ptr++;
      }
    }
  }

  if (debug) print("before blur min/max = %f %f\n", min_val, max_val);
  

  /* note data is stored by rows (along x), then by cols (along y) then slices (along z) */
  
  /*-----------------------------------------------------------------------------*/
  /*             determine   size of data structures needed                      */
  
  vector_size_data = sizes[xyzv[VIO_X]]; 
  kernel_size_data = (int)(((4*fwhmx)/VIO_ABS(steps[xyzv[VIO_X]])) + 0.5);
  
  if (kernel_size_data > MAX(vector_size_data,256))
    kernel_size_data =  MAX(vector_size_data,256);
  
  /*             array_size_pow2 will hold the size of the arrays for FFT convolution,
                 remember that ffts require arrays 2^n in length                          */
  
  array_size_pow2  = next_power_of_two(vector_size_data+kernel_size_data+1);
  array_size = 2*array_size_pow2+1;  /* allocate 2*, since each point is a    */
                                     /* complex number for FFT, and the plus 1*/
                                     /* is for the zero offset FFT routine    */
  ALLOC(dat_vector, array_size);
  ALLOC(dat_vecto2, array_size);
  ALLOC(kern,       array_size);

  /*--------------------------------------------------------------------------------------*/
  /*                get ready to start up the transformation.                             */
  /*--------------------------------------------------------------------------------------*/
  
  initialize_progress_report( &progress, FALSE, sizes[xyzv[VIO_Z]] + sizes[xyzv[VIO_Z]] + sizes[xyzv[VIO_X]] + 1,
                             "Blurring volume" );
  
  /*--------------------------------------------------------------------------------------*/
  /*                start with rows - i.e. the d/dx volume                                */
  /*--------------------------------------------------------------------------------------*/
  
  /*    1st calculate kern array for gaussian kernel*/
  
  make_kernel(kern,(float)(VIO_ABS(steps[xyzv[VIO_X]])),fwhmx,array_size_pow2,kernel_type);
  fft1(kern,array_size_pow2,1);
  
  /*    calculate offset for original data to be placed in vector            */
  
  data_offset = (array_size_pow2-sizes[xyzv[VIO_X]])/2;
  
  /*    2nd now convolve this kernel with the rows of the dataset            */
  
  slice_limit = fwhmx > 0 ? sizes[xyzv[VIO_Z]] : 0;

  for (slice = 0; slice < slice_limit; slice++) {      /* for each slice */
    
    for (row = 0; row < sizes[xyzv[VIO_Y]]; row++) {           /* for each row   */
      
      f_ptr = fdata + slice*slice_size + row*sizes[xyzv[VIO_X]];
      memset(dat_vector,0,(2*array_size_pow2+1)*sizeof(float));
      
      for (col=0; col< sizes[xyzv[VIO_X]]; col++) {        /* extract the row */
        dat_vector[1 +2*(col+data_offset)  ] = *f_ptr++;
      }
      
      fft1(dat_vector,array_size_pow2,1);
      muli_vects(dat_vecto2,dat_vector,kern,array_size_pow2);
      fft1(dat_vecto2,array_size_pow2,-1);
      
      f_ptr = fdata + slice*slice_size + row*sizes[xyzv[VIO_X]];
      for (col=0; col< sizes[xyzv[VIO_X]]; col++) {        /* put the row back */
        
        vindex = 1 + 2*(col+data_offset);
        
        *f_ptr++ = dat_vecto2[vindex]/array_size_pow2;
        

        
      }
      
      
    }
    update_progress_report( &progress, slice+1 );
  }
  
  FREE(dat_vector);
  FREE(dat_vecto2);
  FREE(kern);


  
  /*--------------------------------------------------------------------------------------*/
  /*                 now do cols - i.e. the d/dy volume                                   */
  /*--------------------------------------------------------------------------------------*/
  
  /*-----------------------------------------------------------------------------*/
  /*             determine   size of data structures needed                      */
  
  f_ptr = fdata;
  
  vector_size_data = sizes[xyzv[VIO_Y]];
  kernel_size_data = (int)(((4*fwhmy)/(VIO_ABS(steps[xyzv[VIO_Y]]))) + 0.5);
  
  if (kernel_size_data > MAX(vector_size_data,256))
    kernel_size_data =  MAX(vector_size_data,256);
  
  
  /*             array_size_pow2 will hold the size of the arrays for FFT convolution,
                 remember that ffts require arrays 2^n in length                          */
  
  array_size_pow2  = next_power_of_two(vector_size_data+kernel_size_data+1);
  array_size = 2*array_size_pow2+1;  /* allocate 2*, since each point is a    */
                                     /* complex number for FFT, and the plus 1*/
                                     /* is for the zero offset FFT routine    */
  
  ALLOC(dat_vector, array_size);
  ALLOC(dat_vecto2, array_size);
  ALLOC(kern,       array_size);
  
  /*    1st calculate kern array for gaussian kernel*/
  
  make_kernel(kern,(float)(VIO_ABS(steps[xyzv[VIO_Y]])),fwhmy,array_size_pow2,kernel_type);
  fft1(kern,array_size_pow2,1);
  
  /*    calculate offset for original data to be placed in vector            */
  
  data_offset = (array_size_pow2-sizes[xyzv[VIO_Y]])/2;
  
  /*    2nd now convolve this kernel with the rows of the dataset            */
  
  slice_limit = fwhmy > 0 ? sizes[xyzv[VIO_Z]] : 0;


  for (slice = 0; slice < slice_limit; slice++) {      /* for each slice */
    
    for (col = 0; col < sizes[xyzv[VIO_X]]; col++) {           /* for each col   */
      
      /*         f_ptr = fdata + slice*slice_size + row*sizeof(float); */
      
      f_ptr = fdata + slice*slice_size + col;
      
      
      memset(dat_vector,0,(2*array_size_pow2+1)*sizeof(float));
      
      for (row=0; row< sizes[xyzv[VIO_Y]]; row++) {        /* extract the col */
        dat_vector[1 +2*(row+data_offset) ] = *f_ptr;
        f_ptr += row_size;
      }
      
      
      fft1(dat_vector,array_size_pow2,1);
      muli_vects(dat_vecto2,dat_vector,kern,array_size_pow2);
      fft1(dat_vecto2,array_size_pow2,-1);
      
      f_ptr = fdata + slice*slice_size + col;
      for (row=0; row< sizes[xyzv[VIO_Y]]; row++) {        /* put the col back */
        
        vindex = 1 + 2*(row+data_offset);
        
        *f_ptr = dat_vecto2[vindex]/array_size_pow2;
        
        f_ptr += row_size;
        
        
      }
      
    }
    update_progress_report( &progress, slice+sizes[xyzv[VIO_Z]]+1 );
    
  }
  
  FREE(dat_vector);
  FREE(dat_vecto2);
  FREE(kern);
  
  
  
  /*--------------------------------------------------------------------------------------*/
  /*                 now do slices - i.e. the d/dz volume                                 */
  /*--------------------------------------------------------------------------------------*/
  
  /*-----------------------------------------------------------------------------*/
  /*             determine   size of data structures needed                      */
  
  
  f_ptr = fdata;
  
  vector_size_data = sizes[xyzv[VIO_Z]];
  kernel_size_data = (int)(((4*fwhmz)/(VIO_ABS(steps[xyzv[VIO_Z]]))) + 0.5);
  
  if (kernel_size_data > MAX(vector_size_data,256))
    kernel_size_data =  MAX(vector_size_data,256);
  
  /*             array_size_pow2 will hold the size of the arrays for FFT convolution,
                 remember that ffts require arrays 2^n in length                          */
  
  array_size_pow2  = next_power_of_two(vector_size_data+kernel_size_data+1);
  array_size = 2*array_size_pow2+1;  /* allocate 2*, since each point is a    */
                                     /* complex number for FFT, and the plus 1*/
                                     /* is for the zero offset FFT routine    */

  ALLOC(dat_vector, array_size); 
  ALLOC(dat_vecto2, array_size); 
  ALLOC(kern,       array_size); 
  
  max_val = -FLT_MAX;
  min_val = FLT_MAX;
    
  if ( fwhmz > 0 ){
    
    /*    1st calculate kern array for gaussian kernel*/
    
    make_kernel(kern,(float)(VIO_ABS(steps[xyzv[VIO_Z]])),fwhmz,array_size_pow2,kernel_type);
    fft1(kern,array_size_pow2,1);
    
    /*    calculate offset for original data to be placed in vector            */
    
    data_offset = (array_size_pow2-sizes[xyzv[VIO_Z]])/2;
    
    
    /*    2nd now convolve this kernel with the slices of the dataset            */
    
    for (col = 0; col < sizes[xyzv[VIO_X]]; col++) {      /* for each column */
      
      for (row = 0; row < sizes[xyzv[VIO_Y]]; row++) {           /* for each row   */
        
        f_ptr = fdata + col*col_size + row;
        
        memset(dat_vector,0,(2*array_size_pow2+1)*sizeof(float));
        
        for (slice=0; slice< sizes[xyzv[VIO_Z]]; slice++) {        /* extract the slice vector */
          dat_vector[1 +2*(slice+data_offset) ] = *f_ptr;
          f_ptr += slice_size;
        }
        
        fft1(dat_vector,array_size_pow2,1);
        muli_vects(dat_vecto2,dat_vector,kern,array_size_pow2);

        fft1(dat_vecto2,array_size_pow2,-1);
        
        f_ptr = fdata + col*col_size + row;
        
        for (slice=0; slice< sizes[xyzv[VIO_Z]]; slice++) {        /* put the vector back */
          
          vindex = 1 + 2*(slice+data_offset);
          
          *f_ptr = dat_vecto2[vindex]/array_size_pow2;
          
          if (max_val<*f_ptr) max_val = *f_ptr;
          if (min_val>*f_ptr) min_val = *f_ptr;
          
          f_ptr += slice_size;
        }
        
        
      }
      update_progress_report( &progress, col + 2*sizes[xyzv[VIO_Z]] + 1 );
      
    }
    
  }  /* if ndim */
  else {

    for (slice = 0; slice < sizes[xyzv[VIO_Z]]; slice++) {      /* for each slice */
      for (col = 0; col < sizes[xyzv[VIO_X]]; col++) {             /* for each column */
        for (row = 0; row < sizes[xyzv[VIO_Y]]; row++) {           /* for each row   */
          if (max_val<*f_ptr) max_val = *f_ptr;
          if (min_val>*f_ptr) min_val = *f_ptr;
          f_ptr++;
        }
      }
    }
  }



  terminate_progress_report( &progress );
  
  if (debug) print("after  blur min/max = %f %f\n", min_val, max_val);
  
  FREE(dat_vector);
  FREE(dat_vecto2);
  FREE(kern);
  
  if (reals_fp != (FILE *)NULL) {
    status = io_binary_data(reals_fp,WRITE_FILE, fdata, sizeof(float), total_voxels);
    if (status != VIO_OK) 
      print_error_and_line_num("problems writing blurred reals data...",__FILE__, __LINE__);
  }
  
/* set up the correct info to copy the data back out in mnc */

  f_ptr = fdata;
  
  set_volume_real_range(data, min_val, max_val);

  printf("Making byte volume...\n" );
  for(slice=0; slice<sizes[xyzv[VIO_Z]]; slice++) {
    pos[xyzv[VIO_Z]] = slice;
    for(row=0; row<sizes[xyzv[VIO_Y]]; row++) {
      pos[xyzv[VIO_Y]] = row;
      for(col=0; col<sizes[xyzv[VIO_X]]; col++) {
        pos[xyzv[VIO_X]] = col;
        tmp = CONVERT_VALUE_TO_VOXEL(data, *f_ptr);
         SET_VOXEL_3D( data, pos[0], pos[1], pos[2], tmp);
        f_ptr++;
      }
    }
  }

  FREE(fdata);
  
  snprintf(full_outfilename, sizeof(full_outfilename), "%s_blur.mnc",outfile);

  status = output_modified_volume(full_outfilename, NC_UNSPECIFIED, FALSE, 
                                  min_val, max_val, data, infile, history, 
                                  (minc_output_options *)NULL);


  if (status != VIO_OK)
    print_error_and_line_num("problems writing blurred data...",__FILE__, __LINE__);

  return(status);

  
}

