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
   @CREATED    : January 25, 1992 louis collins (Original using .iff files)
   @MODIFIED   : Wed Jun 23 09:04:34 EST 1993 Louis Collins
        rewrite using mnc files and David Macdonald's libmni.a
   ---------------------------------------------------------------------------- */

#include <def_mni.h>
#include <minc.h>
#include <ParseArgv.h>
#include <mincblur.h>

char *prog_name;
int  debug;
int  verbose;

extern int ms_volume_reals_flag;

main (int argc, char *argv[] )
{   
  int 
    dimensions,
    header_size,
    gradonlyflg,bluronlyflg,
    argn;
  
  char 
    *infilename,
    *outfilename,
    blur_datafile[250];
  FILE 
    *ifd,*ofd;
  double   
    fwhm,
    standard;
  float
    kernel1;
  
  Status 
    status;
  
  Volume
    data;
  real
    step[3];
  int
    sizes[3];

  static ArgvInfo argTable[] = {
    {"-fwhm", ARGV_FLOAT, (char *) 0, (char *) &fwhm,
       "Full-width-half-maximum of gaussian kernel"},
    {"-standarddev", ARGV_FLOAT, (char *) 0, (char *) &standard,
       "Standard deviation of gaussian kernel"},
    {"-dimensions", ARGV_INT, (char *) 0, (char *) &dimensions,
       "Number of dimensions to blur (either 2 or 3)."},
    {NULL, ARGV_HELP, NULL, NULL,
       "Program flags."},
    {"-blur_only", ARGV_CONSTANT, (char *) TRUE, (char *) &bluronlyflg,
       "Create only the blurred volume."},
    {"-grad_only", ARGV_CONSTANT, (char *) TRUE, (char *) &gradonlyflg, 
       "Create only the gradient volumes."},
    {"-no_reals", ARGV_CONSTANT, (char *) FALSE, (char *) &ms_volume_reals_flag, 
       "Do not write out the real (float) data."},
    {NULL, ARGV_HELP, NULL, NULL,
       "Options for logging progress. Default = -verbose."},
    {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
       "Write messages indicating progress"},
    {"-quiet", ARGV_CONSTANT, (char *) FALSE , (char *) &main_args.flags.verbose,
       "Do not write log messages"},
    {"-debug", ARGV_CONSTANT, (char *) TRUE, (char *) &main_args.flags.debug,
       "Print out debug info."},
    {NULL, ARGV_END, NULL, NULL, NULL}
  };
  

  
  /* set default values */
  
  prog_name = argv[0];
  gradonlyflg = bluronlyflg = FALSE;
  fwhm = std = 0.0;
  kernel1 = 0.0;
  infilename =  outfilename = NULL;
  ifd = ofd = NULL;
  dimensions = 3;
  verbose = TRUE;
  debug   = FALSE;


  /* Call ParseArgv to interpret all command line args */

  if (ParseArgv(&argc, argv, argTable, 0) || (argc!=3)) {
    (void) fprintf(stderr, 
		   "\nUsage: %s [<options>] <inputfile> <output_basename>\n", 
		   prog_name);
    (void) fprintf(stderr,"       %s [-help]\n\n", prog_name);
    exit(EXIT_FAILURE);
  }

  if (standard==0.0 && fwhm==0.0) {
    print_error ("Must specify either -fwhm or -standard on command line.", 
		 __FILE__, __LINE__, 0, 0,0,0,0);
  }

  infilename  = argv[1];	/* set up necessary file names */
  outfilename = argv[2]; 

				/* check to see if the output file can be written */
  status = open_file( outfilename , WRITE_FILE, BINARY_FORMAT,  &ofd );
  if ( status != OK ) 
    print_error ("filename `%s' cannot be opened.", __FILE__, __LINE__, outfilename, 0,0,0,0);
  status = close_file(ofd);
  remove(outfilename);   


  if (gradonlyflg  && !ms_volume_reals_flag) {
    print_error ("Must allow reals to be written to be able to calculate gradients.", 
		 __FILE__, __LINE__, 0, 0,0,0,0);
  }
  


  /******************************************************************************/
  /*             create blurred volume first                                    */
  /******************************************************************************/
  
  status = input_volume(infilename, &data);
  if ( status != OK )
    print_error("problems reading `%s'.\n",__FILE__, __LINE__,infilename, 0,0,0,0);
    
  DEBUG_PRINT1 ( "===== Debugging information from %s =====\n", prog_name);
  DEBUG_PRINT1 ( "Data filename     = %s\n", infilename);
  DEBUG_PRINT1 ( "Output basename   = %s\n", outfilename);
  get_volume_sizes(data, sizes);
  get_volume_separations(data, step);
  DEBUG_PRINT3 ( "Input volume      = %3d cols by %3d rows by %d slices\n",
		 sizes[X], sizes[Y], sizes[Z]);
  DEBUG_PRINT3 ( "Input voxels are  = %8.3f %8.3f %8.3f\n", 
		 step[X], step[Y], step[Z]);
  DEBUG_PRINT2 ( "min/max value     = %8.3f %8.3f\n", data->min_value, data->max_value);

  if (data->n_dimensions!=3) {
    print_error ("File %s has %d dimensions.  Only 3 dims supported.", 
		 __FILE__, __LINE__, infilename, data->n_dimensions,0,0,0);
  }

  if (bluronlyflg || !gradonlyflg) {
    status = blur3D_volume(data,
			   fwhm,
			   outfilename,
			   dimensions);
    
    if (status==OK) {
      status = close_file(ifd);
      if (status!=OK)
	print_error("Error closing <%s>.",infilename,0,0,0,0);
    }
    else
      print_error("Problems blurring <%s>.",infilename,0,0,0,0);
  }

  /******************************************************************************/
  /*             calculate d/dx,  d/dy and d/dz volumes                         */
  /******************************************************************************/


   if ((!bluronlyflg || gradonlyflg ) && (ms_volume_reals_flag)) {
      sprintf(blur_datafile,"%s_reals",outfilename);
      
      status = open_file( blur_datafile ,READ_FILE, BINARY_FORMAT,  &ifd );
      if (status!=OK)
	 print_error("Error opening <%s>.",blur_datafile,0,0,0,0);
      
      calc_float_gradient_3d(ifd,data,outfilename,dimensions);
      
      status = close_file(ifd);
      if (status!=OK)
	 print_error("Error closing <%s>.",blur_datafile,0,0,0,0);
      
      /*************************************************************************/
      /*           calculate magnitude of gradient                             */
      /*************************************************************************/
      
      calc_float_gradient_magnitude(outfilename);
   }

   status = close_file(ofd);

   return(status);
   
}


public void print_error(char *s, char * d1, int d2, int d3, int d4, int d5, int d6, int d7)
{
  (void) fprintf(stderr, "Error in %s in file %s, line %d\n",prog_name,d1,d2);
  (void) fprintf(stderr, "   %s\n", s, d3,d4,d5,d6,d7);
  exit(EXIT_FAILURE);
}
