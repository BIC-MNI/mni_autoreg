/* ----------------------------- MNI Header -----------------------------------
   @NAME       : gradmag_volume.c
   @INPUT      : 
   @OUTPUT     : 
   @RETURNS    : 
   @DESCRIPTION: routines used to make blurred and graient volumes
   @METHOD     : 
   @GLOBALS    : 
   @CALLS      : 
   @CREATED    : October 10, 1991 louis collins
   @MODIFIED   : 
	Wed Jun 16 13:21:18 EST 1993 LC
	  rewrite for minc files using libmni.a and libminc.a
   ---------------------------------------------------------------------------- */

#include <def_mni.h>

#define EQUAL3(a,b,c) (((a)==(b)) && ((b)==(c)) && ((a)==(c)))


public calc_gradient_magnitude(char *infilename)
{   
   int argn;
   
   char 
      *usage = "filename",
      infilename1[256],
      infilename2[256],
      infilename3[256],
      outfilename[256],
      time_str[250];
   FILE 
      *ifd1,
      *ifd2,
      *ifd3,
      *ofd;

   unsigned char 
      *in_slice1,
      *in_slice2,
      *in_slice3,
      *out_data,
      *c_ptr1,
      *c_ptr2,
      *c_ptr3,
      *c_out;

   unsigned short 
      *s_ptr1,
      *s_ptr2,
      *s_ptr3,
      *s_out;

   float 
      *f_ptr,
      tx,ty,tz,tv;

   int 
      total_voxels,
      header_size1,
      header_size2,
      header_size3,
      x1,x2,y1,y2,z1,z2,
      value,
      slice_size;

   progress_struct 
      progress;

   register 
      int col,row,slice;

   Status 
      status;

   DATA 
      *data1,*data2,*data3, *data4;
   
   /* set default values */
   
   ifd1 = ifd2 = ifd3 = ofd = NULL;
   
   sprintf (infilename1,"%s_dx.iff",infilename);
   sprintf (infilename2,"%s_dy.iff",infilename);
   sprintf (infilename3,"%s_dz.iff",infilename);
   sprintf (outfilename,"%s_dxyz.iff",infilename);

   status = open_file( infilename1 , READ_FILE, BINARY_FORMAT, &ifd1 );
   if ( status != OK ) 
      print_error ("filename `%s' not found.", __FILE__, __LINE__, infilename1, 0,0,0,0);

   status = open_file( infilename2 , READ_FILE, BINARY_FORMAT, &ifd2 );
   if ( status != OK ) 
      print_error ("filename `%s' not found.", __FILE__, __LINE__, infilename2, 0,0,0,0);

   status = open_file( infilename3 ,READ_FILE, BINARY_FORMAT,  &ifd3 );
   if ( status != OK ) 
      print_error ("filename `%s' not found.", __FILE__, __LINE__, infilename3, 0,0,0,0);

   status = open_file( outfilename , WRITE_FILE, BINARY_FORMAT, &ofd );
   if ( status != OK ) 
      print_error ("cannot write to filename `%s'.", __FILE__, __LINE__, outfilename, 0,0,0,0);


  
   ALLOC1(status, data1,1,DATA); 
   ALLOC1(status, data2,1,DATA); 
   ALLOC1(status, data3,1,DATA); 
   if (status != OK) 
      PRINT("problems allocing 3 slices of data...");

   status = read_iffheader(ifd1,data1,&header_size1);
   if (status!=OK) {
      print_error("error reading header of <%s>.",infilename1,0,0,0,0);
      exit(-1);
   }

   status = read_iffheader(ifd2,data2,&header_size2);
   if (status!=OK) {
      print_error("error reading header of <%s>.",infilename2,0,0,0,0);
      exit(-1);
   }

   status = read_iffheader(ifd3,data3,&header_size3);
   if (status!=OK) {
      print_error("error reading header of <%s>.",infilename3,0,0,0,0);
      exit(-1);
   }
   
  if (!EQUAL3(data1->rows,           data2->rows,           data3->rows)   ||
      !EQUAL3(data1->cols,           data2->cols,           data3->cols)   ||
      !EQUAL3(data1->slices,         data2->slices,         data3->slices)  ||
      !EQUAL3(data1->pix_depth,      data2->pix_depth,      data3->pix_depth)    ||
      !EQUAL3(data1->bytes_per_voxel,data2->bytes_per_voxel,data3->bytes_per_voxel)  ) {
printf ("rows: %d %d %d\n",data1->rows,           data2->rows,           data3->rows);
printf ("cols: %d %d %d\n",data1->cols,           data2->cols,           data3->cols);
printf ("slices: %d %d %d\n",data1->slices,           data2->slices,           data3->slices);
printf ("pix_depth: %d %d %d\n",data1->pix_depth,           data2->pix_depth,           data3->pix_depth);
printf ("bytes_per_voxel: %d %d %d\n",data1->bytes_per_voxel,  data2->bytes_per_voxel,   data3->bytes_per_voxel);
      print_error("headers do not match.",__FILE__,__LINE__,0,0,0,0,0);
      exit(-1);
   }

   slice_size = data1->slice_size;
   ALLOC1(status, in_slice1, slice_size, char);
   ALLOC1(status, in_slice2, slice_size, char);
   ALLOC1(status, in_slice3, slice_size, char);
   if (status != OK) 
      PRINT("problems allocing 3 slices of data...");


   ALLOC1(status, data4,1,DATA);	/* output data struct  */
   if (status != OK) 
      PRINT("problems allocing DATA...");

   *data4 = *data1;

   ALLOC1(status, out_data, slice_size, char);
   if (status != OK) 
      PRINT("problems allocing out_data slice...");

   total_voxels = data1->slices * data1->rows * data1->cols;



   ALLOC1(status, fdata, total_voxels, float); 
      if (status != OK) 
      PRINT("problems allocing fDATA voxels...");

   initialize_progress_report( &progress, TRUE, data1->slices, "Gradient Magnitude:");
   
   f_ptr = fdata;
   *f_ptr = 0;
   data4->fp_max = -1000000.0;
   data4->fp_min =  1000000.0;
   if (data1->bytes_per_voxel==1) {
      for (slice = 0; slice < data1->slices; ++slice) {

	 status = io_binary_data( ifd1, READ_FILE, in_slice1, sizeof(char), slice_size);
	 if ( status != OK ) 
	    print_error ("problems reading slice %d of `%s'.", __FILE__, __LINE__, slice, 
			 infilename1, 0,0,0);
	 status = io_binary_data( ifd2, READ_FILE, in_slice2, sizeof(char), slice_size);
	 if ( status != OK ) 
	    print_error ("problems reading slice %d of `%s'.", __FILE__, __LINE__, slice, 
			 infilename2, 0,0,0);
	 status = io_binary_data( ifd3, READ_FILE, in_slice3, sizeof(char), slice_size);
	 if ( status != OK ) 
	    print_error ("problems reading slice %d of `%s'.", __FILE__, __LINE__, slice, 
			 infilename3, 0,0,0);
	 
	 c_ptr1 = in_slice1;
	 c_ptr2 = in_slice2;
	 c_ptr3 = in_slice3;

	 for (col = 0; col < data1->cols; ++col)
	    for (row = 0; row < data1->rows; ++ row) {

	       tx = get_float_DATA(*c_ptr1,data1);
	       ty = get_float_DATA(*c_ptr2,data2);
	       tz = get_float_DATA(*c_ptr3,data3);

	       if (slice==0  && col==0 && row==0) {
		  *f_ptr = 0;
	       }
	       else {
		  *f_ptr  = fsqrt(tx*tx + ty*ty + tz*tz);
	       }
	       if (data4->fp_max < *f_ptr) data4->fp_max = *f_ptr;
	       if (data4->fp_min > *f_ptr) data4->fp_min = *f_ptr;
	       
	       c_ptr1++;
	       c_ptr2++;
	       c_ptr3++;
	       f_ptr++;
	    }

	 update_progress_report( &progress, slice+1 );
      }

      status = write_iffheader(ofd, data4);  /* copy header of data4 to result file. */


      c_ptr1  = out_data;
      f_ptr = fdata;

      for (slice = 0; slice < data4->slices; ++slice) {

	 c_ptr1  = out_data;

	 for (col = 0; col < data4->cols; ++col)
	    for (row = 0; row < data4->rows; ++ row) {

	       *c_ptr1 = get_uchar_DATA(*f_ptr,data4);
	       c_ptr1++;
	       f_ptr++;
	    }

	 status = io_binary_data( ofd, WRITE_FILE, out_data, sizeof(char), slice_size);
	 if ( status != OK ) 
	    print_error ("problems writing slice %d of `%s'.", __FILE__, __LINE__, slice, 
			 outfilename, 0,0,0);
      }
      



   }
   else {
      for (slice = 0; slice < data1->slices; ++slice) {

	 status = io_binary_data( ifd1, READ_FILE, in_slice1, sizeof(char), slice_size);
	 if ( status != OK ) 
	    print_error ("problems reading slice %d of `%s'.", __FILE__, __LINE__, slice, 
			 infilename1, 0,0,0);
	 status = io_binary_data( ifd2, READ_FILE, in_slice2, sizeof(char), slice_size);
	 if ( status != OK ) 
	    print_error ("problems reading slice %d of `%s'.", __FILE__, __LINE__, slice, 
			 infilename2, 0,0,0);
	 status = io_binary_data( ifd3, READ_FILE, in_slice3, sizeof(char), slice_size);
	 if ( status != OK ) 
	    print_error ("problems reading slice %d of `%s'.", __FILE__, __LINE__, slice, 
			 infilename3, 0,0,0);
	 
	 s_ptr1 = (unsigned short *)in_slice1;
	 s_ptr2 = (unsigned short *)in_slice2;
	 s_ptr3 = (unsigned short *)in_slice3;
	 s_out = s_ptr3;

	 for (col = 0; col < data1->cols; ++col)
	    for (row = 0; row < data1->rows; ++ row) {
	       tx = (float)(*s_ptr1);
	       tx = tx - 128.0;
	       ty = (float)(*s_ptr2);
	       ty = ty - 128.0;
	       tz = (float)(*s_ptr3);
	       tz = tz - 128.0;
	       tv  = fsqrt(tx*tx + ty*ty + tz*tz);
	       value = (int)(tv+0.5);

	       if (value<0) 
		  *s_out = 0;
	       else if (value>(2<<16))
		  *s_out = (2<<16);
	       else
		  *s_out = (short)value;
	       s_ptr1++;
	       s_ptr2++;
	       s_ptr3++;
	    }

	 status = io_binary_data( ofd, WRITE_FILE, out_data, sizeof(char), slice_size);
	 if ( status != OK ) 
	    print_error ("problems writing slice %d of `%s'.", __FILE__, __LINE__, slice, 
			 outfilename, 0,0,0);
	 

	 update_progress_report( &progress, slice+1 );
      }
   }

   terminate_progress_report( &progress );


   FREE1(status, data1); 
   FREE1(status, data2); 
   FREE1(status, data3); 
   FREE(in_slice1);
   FREE(in_slice2);
   FREE(in_slice3);
   FREE(out_data);


   status = close_file(ifd1);
   status = close_file(ifd2);
   status = close_file(ifd3);
   status = close_file(ofd);

   FREE1(status, fdata);

   return(status);
   
}

