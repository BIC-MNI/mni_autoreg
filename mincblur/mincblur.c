/* ----------------------------- MNI Header -----------------------------------
   @NAME       : mincblur.c
   @INPUT      : -name of file corresponding to minc volume data
                 -basename for output file
                 -size of blurring kernel
   
   @OUTPUT     : one or more of five volumes, depending on the command-line args:

        1 - blurred volume - volume blurred by gaussian kernel (sigma = k/2.36)
        2 - d/dx           - derivative along x of blurred volume.
        3 - d/dy           - derivative along y of blurred volume.
        4 - d/dz           - derivative along z of blurred volume.
        5 - | d^3/dxdydz | = sqrt ( (d/dx)^2 + (d/dy)^2 + (d/dz)^2 ) 
                           i.e. the derivative magnitude of the blurred volume.

   @RETURNS    : TRUE if ok, VIO_ERROR if error.

   @DESCRIPTION: This program will read in a volumetric dataset in
                 .mnc format.  This file is assumed to contain rectangular 
                 voxels.

   @METHOD     : The blurred volume calculated in the FOURIER domain by 
                 multiplying the FT(data) by the FT(Gaussian).

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
   @CREATED    : January 25, 1992 louis collins (Original using .iff files)
   @MODIFIED   : $Log: mincblur.c,v $
   @MODIFIED   : Revision 96.7  2009-07-23 22:34:00  claude
   @MODIFIED   : cleanup in mincblur for VIO_X, VIO_Y, VIO_Z
   @MODIFIED   :
   @MODIFIED   : Revision 96.6  2009/06/05 20:49:52  claude
   @MODIFIED   : Free memory after usage in mincblur
   @MODIFIED   :
   @MODIFIED   : Revision 96.5  2009/02/13 04:14:40  rotor
   @MODIFIED   :  * small updated to arguments of mincblur
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
   @MODIFIED   : Revision 96.2  2000/02/20 04:01:00  stever
   @MODIFIED   : * use new history_string() function to generate history strings
   @MODIFIED   :   when outputting MNI files (.mnc, .xfm)
   @MODIFIED   : * removed unused vax routines from Proglib
   @MODIFIED   : * tuned configure script; CPPFLAGS and LDFLAGS are now left alone,
   @MODIFIED   :   for the installer to use
   @MODIFIED   :
   @MODIFIED   : Revision 96.1  2000/01/27 18:03:52  louis
   @MODIFIED   : final checkin before switch to CVS
   @MODIFIED   :
 * Revision 96.0  1996/08/21  18:22:24  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:22:20  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:16:28  louis
 * Never released version 0.95
 *
 * Revision 1.14  1996/08/12  14:16:24  louis
 * Pre-release
 *
 * Revision 1.13  1995/09/18  09:02:42  collins
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
 * Revision 1.12  1995/09/18  06:45:42  collins
 * this file is a working version of mincblur.  All references to numerical
 * recipes routines have been removed.  This version is included in the
 * package mni_reg-0.1i
 *
   Wed Jun 23 09:04:34 EST 1993 Louis Collins
        rewrite using mnc files and David Macdonald's libmni.a
   ---------------------------------------------------------------------------- */
#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/mincblur/mincblur.c,v 96.7 2009-07-23 22:34:00 claude Exp $";
#endif

#include <config.h>
#include <float.h>
#include <volume_io.h>
#include <Proglib.h>
#include <ParseArgv.h>
#include <minc.h>
#include "kernel.h"
#include "mincblur.h"

/* confiscated from minctracc.c */
void get_volume_XYZV_indices(VIO_Volume data, int xyzv[]) {

  int     axis, i, vol_dims;
  char ** data_dim_names;

  vol_dims       = get_volume_n_dimensions(data);
  data_dim_names = get_volume_dimension_names(data);

  for(i=0; i<VIO_N_DIMENSIONS+1; i++) xyzv[i] = -1;
  for(i=0; i<vol_dims; i++) {
    if (convert_dim_name_to_spatial_axis(data_dim_names[i], &axis )) {
      xyzv[axis] = i;
    } else {     /* not a spatial axis */
      xyzv[VIO_Z+1] = i;
    }
  }

  delete_dimension_names(data, data_dim_names);
}


int main (int argc, char *argv[] )
{   
  
  FILE 
    *ofd, *reals_fp;
 
  char 
    *infilename,
    *partials_name,
    *output_basename,
    *temp_basename,
    tempname[1024],
    reals_filename[1024];
  VIO_Status 
    status;
  
  VIO_Volume
    data;
  VIO_Real
    min_value, max_value,
    step[3];
  int
    n_dimensions,
    i,
    sizes[3],
    xyzv[VIO_MAX_DIMENSIONS];
  char *history;
  char *tname;

  /* set default values */
  prog_name            = argv[0];
  apodize_data_flg     = TRUE;
  clobber_flag         = FALSE;
  verbose              = TRUE;
  debug                = FALSE;
  do_gradient_flag     = FALSE;
  do_partials_flag     = FALSE;
  infilename           = (char *)NULL;
  output_basename      = (char *)NULL;
  temp_basename        = (char *)NULL;
  partials_name        = (char *)NULL;
  ofd                  = (FILE *)NULL;
  reals_fp             = (FILE *)NULL;
  dimensions           = 3;
  kernel_type          = KERN_GAUSSIAN;
                                /* init kernel size */
  standard             =  0.0;
  fwhm                 =  0.0;
  for(i=0; i<3; i++){
    fwhm_3D[i] = -DBL_MAX;
  } 
  
  history = history_string(argc, argv);


   /******************************************************************************/
  /* Call ParseArgv to interpret all command line args                          */

  if (ParseArgv(&argc, argv, argTable, 0) || (argc!=3)) {
    (void) fprintf(stderr, 
                   "\nUsage: %s [<options>] <inputfile> <output_basename>\n", 
                   prog_name);
    (void) fprintf(stderr,"       %s [-help]\n\n", prog_name);
    exit(EXIT_FAILURE);
  }

  /******************************************************************************/
  /* find the size of the blurring kernel, from one of -std, -fwhm, -fwhm3d     */

  if (standard==0.0 && fwhm==0.0 && 
      fwhm_3D[0]==-DBL_MAX && fwhm_3D[1]==-DBL_MAX && fwhm_3D[2]==-DBL_MAX ) {
    print_error_and_line_num ("Must specify either -fwhm, -3D_fwhm or -standard on command line.\n", 
                 __FILE__, __LINE__);
  }

  if (fwhm==0.0) fwhm=standard*2.35;
  
  if (fwhm !=0.0 ) {
    for(i=0; i<3; i++) fwhm_3D[i] = fwhm;
  };                                

  /******************************************************************************/
  /*                   set up necessary file names                              */
  /******************************************************************************/
  infilename      = argv[1];        
  output_basename = argv[2]; 

                                /* check to see if the output file can be written */

  if (!clobber_flag ) {
    tname = malloc(strlen(output_basename)+strlen("_blur.mnc")+2 * sizeof(char*));
    strcpy(tname,output_basename);
    strcat(tname,"_blur.mnc");
    if( file_exists(tname)) {
      print ("File %s exists. -- %d\n", tname);
      print ("Use -clobber to overwrite.\n");
      free(tname);
      return VIO_ERROR;
    }
    free(tname);
  }

  status = open_file( output_basename , WRITE_FILE, BINARY_FORMAT,  &ofd );
  if ( status != VIO_OK ) 
    print_error_and_line_num ("filename `%s' cannot be opened.", 
                              __FILE__, __LINE__, output_basename);
  status = close_file(ofd);
  remove(output_basename);   

       /* if any gradient data is needed, then we have to save the blurred
          volume in float representation, otherwise quantization errors can
          mess up the derivatives.  So get temporary file to save 'real data' */

  if ((do_partials_flag || do_gradient_flag)) {

    temp_basename = tempnam(NULL,"mblur");
    
    if (temp_basename == (char *)NULL) {
      print_error_and_line_num ("Can't build temporary filename.", 
                                __FILE__, __LINE__);
    }

    sprintf(reals_filename,"%s_reals.raw",temp_basename);

    status = open_file( reals_filename, WRITE_FILE, BINARY_FORMAT,  &reals_fp );
    
    if (status != VIO_OK) {
      print_error_and_line_num ("Temporary file to save blurred volume cannot be opened.", 
                                __FILE__, __LINE__);
    }

  }


  /******************************************************************************/
  /*             create blurred volume first                                    */
  /******************************************************************************/

  status = input_volume(infilename, VIO_N_DIMENSIONS, 
                        get_default_dim_names( VIO_N_DIMENSIONS ),
                        NC_UNSPECIFIED, FALSE, 0.0, 0.0, TRUE, 
                        &data, (minc_input_options *)NULL);
  if ( status != VIO_OK )
    print_error_and_line_num("problems reading `%s'.\n",__FILE__, __LINE__,infilename);

  get_volume_XYZV_indices( data, xyzv );

  if (debug) {
    get_volume_sizes(data, sizes);
    get_volume_separations(data, step);
    printf ( "===== Debugging information from %s =====\n", prog_name);
    printf ( "Data filename     = %s\n", infilename);
    printf ( "Output basename   = %s\n", output_basename);
    printf ( "Input volume      = %3d cols by %3d rows by %d slices\n",
            sizes[xyzv[VIO_X]], sizes[xyzv[VIO_Y]], sizes[xyzv[VIO_Z]]);
    printf ( "Input voxels are  = %8.3f %8.3f %8.3f\n", 
            step[xyzv[VIO_X]], step[xyzv[VIO_Y]], step[xyzv[VIO_Z]]);
    get_volume_real_range(data,&min_value, &max_value);
    printf ( "min/max value     = %8.3f %8.3f\n", min_value, max_value);
  }
  
  n_dimensions = get_volume_n_dimensions (data);
  if (n_dimensions!=3) {
    print_error_and_line_num ("File %s has %d dimensions.  Only 3 dims supported.", 
                              __FILE__, __LINE__, infilename, n_dimensions);
  }
  
                                /* apodize data if needed */
  if (apodize_data_flg) {
    if (debug) print ("Apodizing data at (%f,%f) (%f,%f) (%f,%f)\n",
                      fwhm_3D[0], fwhm_3D[0], fwhm_3D[1], fwhm_3D[1], fwhm_3D[2], fwhm_3D[2] );
    apodize_data(data, xyzv, fwhm_3D[0], fwhm_3D[0], fwhm_3D[1], fwhm_3D[1], fwhm_3D[2], fwhm_3D[2] );
  }

  // Zero the kernel where we don't want to blur
  if ( dimensions == 2 )
      fwhm_3D[xyzv[2]] = 0;
  else if ( dimensions == 1 )
      fwhm_3D[xyzv[0]] = fwhm_3D[xyzv[1]] = 0;

                                /* now _BLUR_ the DATA! */
  status = blur3D_volume(data, xyzv,
                         fwhm_3D[0],fwhm_3D[1],fwhm_3D[2],
                         infilename,
                         output_basename,
                         reals_fp,
                         kernel_type,history);

  /******************************************************************************/
  /*             calculate d/dx,  d/dy and d/dz volumes                         */
  /******************************************************************************/

  if ((do_partials_flag || do_gradient_flag)) {

                                /* reopen REALS file for READ */
    status = close_file( reals_fp );
    if (status!=VIO_OK)
      print_error_and_line_num("Error closing <%s>.",__FILE__, __LINE__, reals_filename);
    
    status = open_file( reals_filename ,READ_FILE, BINARY_FORMAT,  &reals_fp );
    if (status!=VIO_OK)
      print_error_and_line_num("Error opening <%s>.",__FILE__, __LINE__, reals_filename);
      

    if (do_partials_flag) 
      partials_name = output_basename;
    else 
      partials_name = temp_basename;


    status = gradient3D_volume(reals_fp, data, xyzv, infilename, partials_name, dimensions,
                               history, FALSE);
    if (status!=VIO_OK)
      print_error_and_line_num("Can't calculate the gradient volumes.",__FILE__, __LINE__);


                                /* close and delete the REALS temp data file */
    status = close_file( reals_fp );
    if (status!=VIO_OK)
      print_error_and_line_num("Error closing <%s>.",__FILE__, __LINE__, reals_filename);

    remove_file( reals_filename );
    

    /*************************************************************************/
    /*           calculate and save the gradient magnitude volume            */
    /*************************************************************************/
    
    calc_gradient_magnitude(partials_name, output_basename, history, 
                            &min_value, &max_value);
      
                                /* remove the partial derivs, unless
                                   the user specifically wants to keep
                                   them */
    if (!do_partials_flag) {
      sprintf (tempname,"%s_dx.mnc",partials_name); remove_file(tempname);
      sprintf (tempname,"%s_dy.mnc",partials_name); remove_file(tempname);
      sprintf (tempname,"%s_dz.mnc",partials_name); remove_file(tempname);
    }

  }

  delete_volume( data );
  if( history ) free( history );
  if( temp_basename ) free( temp_basename );

  return(status);
   
}


