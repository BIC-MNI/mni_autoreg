/* ----------------------------- MNI Header -----------------------------------
   @NAME       : mincblur.c
   @INPUT      : name of file corresponding to msp volume data, voxel size aspect 
   ratio, size of filter kernel and output filename.
   
   mincblur inputfile.iff -x 1.0 -y 1.0 -z 2.0 -k 7.0 outputfile
   
   @OUTPUT     : five volumes, 
        1 - blurred volume - volume blurred by gaussian kernel (sigma = k/2.36)
	2 - d/dx           - derivative along x of blurred volume.
	3 - d/dy           - derivative along y of blurred volume.
	4 - d/dz           - derivative along z of blurred volume.
	5 - | d^3/dxdydz | = sqrt ( (d/dx)^2 + (d/dy)^2 + (d/dz)^2 ) 
	                   i.e. the derivative magnitude of the blurred volume.

   @RETURNS    : TRUE if ok, ERROR if error.

   @DESCRIPTION: This program will read in a volumetric dataset in .mnc format.
   This file is assumed to contain rectangular voxels.  

   @METHOD     : The blurred volume calculated in the FOURIER domain by multiplying the
   the FT(data) by FT(Gaussian).

           blurred vol = blur_z ( blur_y( blur_x(data) ) );

   the gradient volumes are calculated also in the fourier domain,
   for each direction,
                      -1
        d/dx data = FT  ( is *  FT(blurred_data) ) ; where i = imag, s = freq var.
                      -1
	d/dy data = FT  ( is *  FT(blurred_data) )
                      -1
	d/dz data = FT  ( is *  FT(blurred_data) )

   the gradient magnitude is calculated so that at each voxel of the gradient volume,
   the intensity is equal to:

       grad_mag = sqrt ( (d/dx)^2 + (d/dy)^2 + (d/dz)^2 ) 

   @GLOBALS    : char *prog_name - stores the name of the program.
                 int  debug      - prints out debugging info   
		 int  verbose    - prints out running info

   @CALLS      : 
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
   @CREATED    : January 25, 1992 louis louis (Original using .iff files)
   @MODIFIED   : $Log: mincblur.c,v $
   @MODIFIED   : Revision 1.12  1995-09-18 06:45:42  louis
   @MODIFIED   : this file is a working version of mincblur.  All references to numerical
   @MODIFIED   : recipes routines have been removed.  This version is included in the
   @MODIFIED   : package mni_reg-0.1i
   @MODIFIED   :
   Wed Jun 23 09:04:34 EST 1993 Louis Collins
        rewrite using mnc files and David Macdonald's libmni.a
   ---------------------------------------------------------------------------- */
#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/mincblur/mincblur.c,v 1.12 1995-09-18 06:45:42 louis Exp $";
#endif

#include <volume_io.h>
#include <ParseArgv.h>
#include <minc.h>
#include "kernel.h"
#include "mincblur.h"
#include <time_stamp.h>
#include <print_error.h>
#include <limits.h>

static char *default_dim_names[N_DIMENSIONS] = { MIzspace, MIyspace, MIxspace };
 


main (int argc, char *argv[] )
{   
  
  FILE 
    *ifd,*ofd;
 
  char 
    *infilename,
    *outfilename,
    blur_datafile[1024];
  
  Status 
    status;
  
  Volume
    data;
  Real
    min_value, max_value,
    step[3];
  int
    i,
    sizes[3];
  char *history;
  char *tname;

  /* set default values */
  apodize_data_flg = ms_volume_reals_flag = TRUE;
  prog_name = argv[0];
  gradonlyflg = curveonlyflg = curvatureflg = bluronlyflg = FALSE;
  fwhm = standard = 0.0;
  for_less(i,0,3) fwhm_3D[i] = -DBL_MAX;
  infilename =  outfilename = NULL;
  ifd = ofd = NULL;
  dimensions = 3;
  verbose = TRUE;
  debug   = FALSE;
  kernel_type = KERN_GAUSSIAN;

  history = time_stamp(argc, argv);

  /* Call ParseArgv to interpret all command line args */

  if (ParseArgv(&argc, argv, argTable, 0) || (argc!=3)) {
    (void) fprintf(stderr, 
		   "\nUsage: %s [<options>] <inputfile> <output_basename>\n", 
		   prog_name);
    (void) fprintf(stderr,"       %s [-help]\n\n", prog_name);
    exit(EXIT_FAILURE);
  }

  if (standard==0.0 && fwhm==0.0 && 
      fwhm_3D[0]==-DBL_MAX && fwhm_3D[1]==-DBL_MAX && fwhm_3D[2]==-DBL_MAX ) {
    print_error_and_line_num ("Must specify either -fwhm, -3D_fwhm or -standard on command line.", 
		 __FILE__, __LINE__);
  }

  if (fwhm==0.0) fwhm=standard;
  
  if (fwhm !=0.0 ) {
    for_less(i,0,3) fwhm_3D[i] = fwhm;
  };				/* else 3D_fwhm has the values from the command line */

  infilename  = argv[1];	/* set up necessary file names */
  outfilename = argv[2]; 

  ALLOC(tname, strlen(outfilename)+strlen("_blur.mnc")+2);
  tname = strcat(tname,outfilename);
  tname = strcat(tname,"_blur.mnc");

  if (!clobber_flag && file_exists(tname)) {
    print ("File %s exists.\n",tname);
    print ("Use -clobber to overwrite.\n");
    return ERROR;
  }

				/* check to see if the output file can be written */
  status = open_file( outfilename , WRITE_FILE, BINARY_FORMAT,  &ofd );
  if ( status != OK ) 
    print_error_and_line_num ("filename `%s' cannot be opened.", __FILE__, __LINE__, outfilename);
  status = close_file(ofd);
  remove(outfilename);   


  if (gradonlyflg  && !ms_volume_reals_flag) {
    print_error_and_line_num ("Must allow reals to be written to be able to calculate gradients.", 
		 __FILE__, __LINE__);
  }
  


  /******************************************************************************/
  /*             create blurred volume first                                    */
  /******************************************************************************/
  
  if (!curveonlyflg) {
    status = input_volume(infilename, 3, default_dim_names, NC_UNSPECIFIED,
			  FALSE, 0.0, 0.0, TRUE, &data, (minc_input_options *)NULL);
    if ( status != OK )
      print_error_and_line_num("problems reading `%s'.\n",__FILE__, __LINE__,infilename);
    
    get_volume_sizes(data, sizes);
    get_volume_separations(data, step);
    
    if (debug) {
      printf ( "===== Debugging information from %s =====\n", prog_name);
      printf ( "Data filename     = %s\n", infilename);
      printf ( "Output basename   = %s\n", outfilename);
      printf ( "Input volume      = %3d cols by %3d rows by %d slices\n",
	      sizes[INTERNAL_X], sizes[INTERNAL_Y], sizes[INTERNAL_Z]);
      printf ( "Input voxels are  = %8.3f %8.3f %8.3f\n", 
	      step[INTERNAL_X], step[INTERNAL_Y], step[INTERNAL_Z]);
      get_volume_real_range(data,&min_value, &max_value);
      printf ( "min/max value     = %8.3f %8.3f\n", min_value, max_value);
    }
    
    if (data->n_dimensions!=3) {
      print_error_and_line_num ("File %s has %d dimensions.  Only 3 dims supported.", 
		   __FILE__, __LINE__, infilename, data->n_dimensions);
    }
    
    
    if (apodize_data_flg) {
      if (debug) print ("Apodizing data at %f\n",fwhm);
      apodize_data(data, fwhm_3D[0], fwhm_3D[0], fwhm_3D[1], fwhm_3D[1], fwhm_3D[2], fwhm_3D[2] );
    }
    

    if (bluronlyflg || !gradonlyflg || curvatureflg) {
      
      status = blur3D_volume(data,
			     fwhm_3D[0],fwhm_3D[1],fwhm_3D[2],
			     infilename,
			     outfilename,
			     dimensions,kernel_type,history);
      
    }

  }    

  /******************************************************************************/
  /*             calculate d/dx,  d/dy and d/dz volumes                         */
  /******************************************************************************/


   if ( (!bluronlyflg || curvatureflg || gradonlyflg ) && (ms_volume_reals_flag)) {
      sprintf(blur_datafile,"%s_reals",outfilename);
      
      status = open_file( blur_datafile ,READ_FILE, BINARY_FORMAT,  &ifd );
      if (status!=OK)
	 print_error_and_line_num("Error opening <%s>.",__FILE__, __LINE__, blur_datafile);
      
      if (!curveonlyflg)
	gradient3D_volume(ifd,data,infilename,outfilename,dimensions,
			  history,FALSE);
      if (!curveonlyflg && curvatureflg) {
	gradient3D_volume(ifd,data,infilename,outfilename,dimensions,
			  history,TRUE);
      }

      status = close_file(ifd);
      if (status!=OK)
	 print_error_and_line_num("Error closing <%s>.",__FILE__, __LINE__, blur_datafile);
    }
      /*************************************************************************/
      /*           calculate magnitude of gradient                             */
      /*************************************************************************/
      
  if ((!bluronlyflg || gradonlyflg || curveonlyflg)) {
    calc_gradient_magnitude(outfilename, history, &min_value, &max_value);

    if (debug)
      print ("max, min grad mag = %f %f\n",min_value, max_value);
  }

  if (!bluronlyflg && !gradonlyflg && (curveonlyflg || curvatureflg) ) {
    print ("gaussian curvature:\n");
    calc_gaussian_curvature(outfilename, history, min_value, max_value);
  }
      
  
  return(status);
   
}


