/* ----------------------------- MNI Header -----------------------------------
@NAME       : blur_volume.c
@INPUT      : data - a pointer to a volume_struct of data
              fwhm - full-width-half-maximum of the gaussian blurring kernel
	      outfile - name of the base filename to store the <name>_blur.mnc
	      ndim - =2, do blurring in the x and y directions only,
                     =3, blur in all three directions.
@OUTPUT     : creates and stores the blurred volume in an output file.
@RETURNS    : status variable - OK or ERROR.
@DESCRIPTION: This routine convolves each row, column and slice with a
              gaussian kernel.  The convolution is accomplished in the fourier
	      domain by multiplying the fourier transformations of both the
	      data and the kernel.
@METHOD     : 
@GLOBALS    : 
@CALLS      : stuff from volume_support.c and libmni.a
@CREATED    : Wed Jun 23 09:04:34 EST 1993  Louis Collins 
                 from code originally written for blur_grad working on .iff files.
@MODIFIED   : 
---------------------------------------------------------------------------- */
#include <def_mni.h>
#include <volume_support.h>
#include <recipes.h>


public Status blur3D_volume(Volume data,
			    float fwhm, 
			    char *outfile, 
			    int ndim)
{ 
  float 
    *fdata,			/* floating point storage for blurred volume */
    *f_ptr,			/* pointer to fdata */
    tmp,
    max_val, 
    min_val,
    
    *dat_vector,		/* temp storage of original row, col or slice vect. */
    *dat_vecto2,		/* storage of result of dat_vector*kern                 */
    *kern,			/* convolution kernel                               */
    k1,				/* width (in pixels) of gaussian                    */
    val;			/* temp storage, to take value from dat_vector, and     */
				/*  place it back into data->voxels                 */
  int				
    i, total_voxels,		
    vector_size_data,		/* original size of row, col or slice vector        */
    kernel_size_data,		/* original size of kernel vector                   */
    array_size_pow2;		/* actual size of vector/kernel data used in FFT    */
				/* routines - needs to be a power of two            */
  int   
    data_offset;		/* offset required to place original data (size n)  */
				/*  into array (size m=2^n) so that data is centered*/

  unsigned char 
    *c_ptr;			/* pointer to access original data                  */
  unsigned short 
    *s_ptr;			/* pointer to access original data                  */
  
  register int 
    slice_limit,
    row,col,slice,		/* counters to access original data                 */
    k,kindex,			/* counters to access kernel                        */
    vindex;			/* counter to access vector and vecto2              */

  int 
    voxel_size, slice_size,	/* size of each data step - in bytes              */
    row_size, col_size;
                
  float
    sum,			/* sum for normalization of kernel                  */
    t1,t2,t3;			/* temp values for calculation of kernel            */

  FILE 
    *f,				/* file used to print out debug info about kernel(s)*/
    *ofp;			/* file used tp write out dx,dy or dz volume        */
  char
    full_outfilename[256];	/* name of output file */

  progress_struct 
    progress;			/* used to monitor progress of calculations         */

  Status 
    status;
  
  int
    sizes[3];			/* number of rows, cols and slices */
  Real
    steps[3];			/* size of voxel step from center to center in x,y,z */

  /*---------------------------------------------------------------------------------*/
  /*             start by setting up the raw data.                                   */
  /*---------------------------------------------------------------------------------*/

  get_volume_sizes(data, sizes);          /* rows,cols,slices */
  get_volume_separations(data, steps);
  
  slice_size = sizes[Y] * sizes[X];    /* sizeof one slice  */
  col_size   = sizes[Y];               /* sizeof one column */
  row_size   = sizes[X];               /* sizeof one row    */
  voxel_size = sizeof(float);                           
  
  total_voxels = sizes[Y]*sizes[X]*sizes[Z];
  
  ALLOC(fdata, total_voxels);
  
  printf("Making float volume..." );
  f_ptr = fdata;
  for_less( slice, 0, sizes[Z])
    for_less( row, 0, sizes[Y])
      for_less( col, 0, sizes[X]) {
	GET_VOXEL_3D( tmp, data, col, row, slice);
	*f_ptr++ = CONVERT_VOXEL_TO_VALUE(data, tmp);
      }
  printf("done\n");
  
  
  /* note data is stored by rows (along x), then by cols (along y) then slices (along z) */
  
  /*-----------------------------------------------------------------------------*/
  /*             determine   size of data structures needed                      */
  
  vector_size_data = sizes[X]; 
  kernel_size_data = (int)(((4*fwhm)/steps[X]) + 0.5);
  
  if (kernel_size_data > MAX(vector_size_data,256))
    kernel_size_data =  MAX(vector_size_data,256);
  
  /*             array_size_pow2 will hold the size of the arrays for FFT convolution,
		 remember that ffts require arrays 2^n in length                          */
  
  array_size_pow2  = next_power_of_two(vector_size_data+kernel_size_data);

  dat_vector = vector(0,2*array_size_pow2+1); /* allocate 2*, since each point is a    */
  dat_vecto2 = vector(0,2*array_size_pow2+1); /* complex number for FFT, and the plus 1*/
  kern       = vector(0,2*array_size_pow2+1); /* is for Num. Rec. FFT routines         */

  /*--------------------------------------------------------------------------------------*/
  /*                get ready to start up the transformation.                             */
  /*--------------------------------------------------------------------------------------*/
  
  initialize_progress_report( &progress, FALSE, sizes[Z] + sizes[Z] + sizes[X] + 1,
			     "Blurring volume" );
  
  /*--------------------------------------------------------------------------------------*/
  /*                start with rows - i.e. the d/dx volume                                */
  /*--------------------------------------------------------------------------------------*/
  
  /*    1st calculate kern array for gaussian kernel*/
  
  make_kernel(kern,steps[X],fwhm,array_size_pow2);
  four1(kern,array_size_pow2,1);
  
  /*    calculate offset for original data to be placed in vector            */
  
  data_offset = (array_size_pow2-sizes[X])/2;
  max_val = -100000.0;
  min_val = 100000.0;
  
  
  /*    2nd now convolve this kernel with the rows of the dataset            */
  
  
  if (ndim<3)
    slice_limit = 1;
  else
    slice_limit = sizes[Z];
  
  for (slice = 0; slice < slice_limit; slice++) {      /* for each slice */
    
    for (row = 0; row < sizes[Y]; row++) {           /* for each row   */
      
      f_ptr = fdata + slice*slice_size + row*sizes[X];
      memset(dat_vector,0,(2*array_size_pow2+1)*sizeof(float));
      
      for (col=0; col< sizes[X]; col++) {        /* extract the row */
	dat_vector[1 +2*(col+data_offset)  ] = *f_ptr++;
      }
      
      four1(dat_vector,array_size_pow2,1);
      muli_vects(dat_vecto2,dat_vector,kern,array_size_pow2);
      four1(dat_vecto2,array_size_pow2,-1);
      
      f_ptr = fdata + slice*slice_size + row*sizes[X];
      for (col=0; col< sizes[X]; col++) {        /* put the row back */
	
	vindex = 1 + 2*(col+data_offset);
	
	*f_ptr++ = dat_vecto2[vindex]/array_size_pow2;
	
	
      }
      
      
    }
    update_progress_report( &progress, slice+1 );
  }
  
  free_vector(dat_vector, 0,2*array_size_pow2+1); 
  free_vector(dat_vecto2, 0,2*array_size_pow2+1); 
  free_vector(kern      , 0,2*array_size_pow2+1); 
    
  
  /*--------------------------------------------------------------------------------------*/
  /*                 now do cols - i.e. the d/dy volume                                   */
  /*--------------------------------------------------------------------------------------*/
  
  /*-----------------------------------------------------------------------------*/
  /*             determine   size of data structures needed                      */
  
  c_ptr = data->voxels;
  s_ptr = (unsigned short *)data->voxels;
  f_ptr = fdata;
  
  vector_size_data = sizes[Y];
  kernel_size_data = (int)(((4*fwhm)/steps[Y]) + 0.5);
  
  if (kernel_size_data > MAX(vector_size_data,256))
    kernel_size_data =  MAX(vector_size_data,256);
  
  
  /*             array_size_pow2 will hold the size of the arrays for FFT convolution,
		 remember that ffts require arrays 2^n in length                          */
  
  array_size_pow2  = next_power_of_two(vector_size_data+kernel_size_data);
  
  dat_vector = vector(0,2*array_size_pow2+1); /* allocate 2*, since each point is a    */
  dat_vecto2 = vector(0,2*array_size_pow2+1); /* complex number for FFT, and the plus 1*/
  kern       = vector(0,2*array_size_pow2+1); /* is for Num. Rec. FFT routines         */
  
  /*    1st calculate kern array for gaussian kernel*/
  
  make_kernel(kern,steps[Y],fwhm,array_size_pow2);
  four1(kern,array_size_pow2,1);
  
  /*    calculate offset for original data to be placed in vector            */
  
  data_offset = (array_size_pow2-sizes[Y])/2;
  
  /*    2nd now convolve this kernel with the rows of the dataset            */
  
  max_val = -100000.0;
  min_val = 100000.0;
  
  if (ndim<3)
    slice_limit = 1;
  else
    slice_limit = sizes[Z];
  
  for (slice = 0; slice < slice_limit; slice++) {      /* for each slice */
    
    for (col = 0; col < sizes[X]; col++) {           /* for each col   */
      
      /*	 f_ptr = fdata + slice*slice_size + row*sizeof(float); */
      
      f_ptr = fdata + slice*slice_size + col;
      
      
      memset(dat_vector,0,(2*array_size_pow2+1)*sizeof(float));
      
      for (row=0; row< sizes[Y]; row++) {        /* extract the col */
	dat_vector[1 +2*(row+data_offset) ] = *f_ptr;
	f_ptr += row_size;
      }
      
      
      four1(dat_vector,array_size_pow2,1);
      muli_vects(dat_vecto2,dat_vector,kern,array_size_pow2);
      four1(dat_vecto2,array_size_pow2,-1);
      
      f_ptr = fdata + slice*slice_size + col;
      for (row=0; row< sizes[Y]; row++) {        /* put the col back */
	
	vindex = 1 + 2*(row+data_offset);
	
	*f_ptr = dat_vecto2[vindex]/array_size_pow2;
	
	f_ptr += row_size;
	
	
      }
      
    }
    update_progress_report( &progress, slice+sizes[Z]+1 );
    
  }
  
  free_vector(dat_vector, 0,2*array_size_pow2+1); 
  free_vector(dat_vecto2, 0,2*array_size_pow2+1); 
  free_vector(kern      , 0,2*array_size_pow2+1); 
  
  
  
  /*--------------------------------------------------------------------------------------*/
  /*                 now do slices - i.e. the d/dz volume                                 */
  /*--------------------------------------------------------------------------------------*/
  
  /*-----------------------------------------------------------------------------*/
  /*             determine   size of data structures needed                      */
  
  
  c_ptr = data->voxels;
  s_ptr = (unsigned short *)data->voxels;
  f_ptr = fdata;
  
  vector_size_data = sizes[Z];
  kernel_size_data = (int)(((4*fwhm)/steps[Z]) + 0.5);
  
  if (kernel_size_data > MAX(vector_size_data,256))
    kernel_size_data =  MAX(vector_size_data,256);
  
  /*             array_size_pow2 will hold the size of the arrays for FFT convolution,
		 remember that ffts require arrays 2^n in length                          */
  
  array_size_pow2  = next_power_of_two(vector_size_data+kernel_size_data);
  
  dat_vector = vector(0,2*array_size_pow2+1); /* allocate 2*, since each point is a    */
  dat_vecto2 = vector(0,2*array_size_pow2+1); /* complex number for FFT, and the plus 1*/
  kern       = vector(0,2*array_size_pow2+1); /* is for Num. Rec. FFT routines         */
  
  if (ndim>2) {
    
    /*    1st calculate kern array for gaussian kernel*/
    
    make_kernel(kern,steps[Z],fwhm,array_size_pow2);
    four1(kern,array_size_pow2,1);
    
    /*    calculate offset for original data to be placed in vector            */
    
    data_offset = (array_size_pow2-sizes[Z])/2;
    
    
    /*    2nd now convolve this kernel with the slices of the dataset            */
    
    max_val = -100000.0;
    min_val = 100000.0;
    
    
    
    for (col = 0; col < sizes[X]; col++) {      /* for each column */
      
      for (row = 0; row < sizes[Y]; row++) {           /* for each row   */
	
	if (data->bytes_per_voxel==1) {
	  /*	    f_ptr = fdata + col*col_size + row*sizeof(float); */
	  f_ptr = fdata + col*col_size + row;
	  
	  memset(dat_vector,0,(2*array_size_pow2+1)*sizeof(float));
	  
	  for (slice=0; slice< sizes[Z]; slice++) {        /* extract the slice vector */
	    dat_vector[1 +2*(slice+data_offset) ] = *f_ptr;
	    f_ptr += slice_size;
	  }
	  
	  four1(dat_vector,array_size_pow2,1);
	  muli_vects(dat_vecto2,dat_vector,kern,array_size_pow2);
	  four1(dat_vecto2,array_size_pow2,-1);
	  
	  f_ptr = fdata + col*col_size + row;
	  
	  for (slice=0; slice< sizes[Z]; slice++) {        /* put the vector back */
	    
	    vindex = 1 + 2*(slice+data_offset);
	    
	    *f_ptr = dat_vecto2[vindex]/array_size_pow2;
	    
	    if (max_val<*f_ptr) max_val = *f_ptr;
	    if (min_val>*f_ptr) min_val = *f_ptr;
	    
	    f_ptr += slice_size;
	  }
	}
	
	
      }
      update_progress_report( &progress, col + 2*sizes[Z] + 1 );
      
    }
    
  }  /* if ndim */
  else {
    max_val = -100000.0;
    min_val = 100000.0;
    
    for (col = 0; col < sizes[X]; col++) {      /* for each column */
      for (row = 0; row < sizes[Y]; row++) {           /* for each row   */
	if (max_val<*f_ptr) max_val = *f_ptr;
	if (min_val>*f_ptr) min_val = *f_ptr;
	f_ptr++;
      }
    }
  }
  terminate_progress_report( &progress );
  
  printf ("blurred min = %f, max = %f\n",min_val, max_val);
  
  free_vector(dat_vector, 0,2*array_size_pow2+1); 
  free_vector(dat_vecto2, 0,2*array_size_pow2+1); 
  free_vector(kern      , 0,2*array_size_pow2+1); 
  
  
  if (ms_volume_reals_flag) {
    sprintf(full_outfilename,"%s_reals",outfile);
    status = open_file(full_outfilename, WRITE_FILE, BINARY_FORMAT, &ofp);
    if (status != OK) 
      PRINT("problems opening blurred reals data...");
    else {
      status = io_binary_data(ofp,WRITE_FILE, fdata, sizeof(float), total_voxels);
      if (status != OK) 
	PRINT("problems writing blurred reals data...");
    }
    close_file(ofp);
  }
  
/* set up the correct info to copy the data back out in mnc */


/*
  *
  *
  *   FINISH THESE WRITE ROUTINES:
  *
  *
  *
*/


  f_ptr = fdata;
  
  data->fp_max = max_val;
  data->fp_min = min_val;
  
  printf("Making byte volume..." );
  for_less( slice, 0, sizes[Z])
    for_less( row, 0, sizes[Y])
      for_less( col, 0, sizes[X]) {
	tmp = CONVERT_VALUE_TO_VOXEL(data, *f_ptr);
 	PUT_VOXEL_3D( data, col, row, slice, tmp);
	*f_ptr++;
      }


  printf("done\n");
  
  sprintf(full_outfilename,"%s_blur.iff",outfile);
  status = write_new_iff_file(data,full_outfilename);
  if (status != OK) 
    PRINT("problems writing blurred data...");
  FREE1(status, data->voxels); FREE1(status, data); 
  
  return(status);


  FREE(fdata);

  
}

