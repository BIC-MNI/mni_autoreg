/* ----------------------------- MNI Header -----------------------------------
@NAME       : mincblur.h
@DESCRIPTION: header/prototype file for mincblur.c

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

@CREATED    : Wed Jun 23 09:04:34 EST 1993 Louis Collins
@MODIFIED   : $Log: mincblur.h,v $
@MODIFIED   : Revision 9.6  1996-08-21 18:22:20  louis
@MODIFIED   : Pre-release
@MODIFIED   :
 * Revision 9.5  1996/08/12  14:16:28  louis
 * Never released version 0.95
 *
 * Revision 1.11  1996/08/12  14:16:24  louis
 * Pre-release
 *
---------------------------------------------------------------------------- */

public Status blur3D_volume(Volume data,
			    double  kernel1, double  kernel2, double  kernel3, 
			    char *infile, 
			    char *outfile, 
			    FILE *reals_fp,
			    int ndim, int kernel_type, char *history);

public Status gradient3D_volume(FILE *ifd, 
				Volume data, 
				char *infile, 
				char *outfile, 
				int ndim,
				char *history,
				int curvature_flg);

public Status calc_Lvv_volume(FILE *ifd, 
			      Volume data, 
			      char *infile,
			      char *outfile, 
			      int ndim,
			      char *history);


public void apodize_data(Volume data, 
			 double xramp1,double xramp2,
			 double yramp1,double yramp2,
			 double zramp1,double zramp2);

public void calc_gradient_magnitude(char *infilename, 
				    char *output_basename,
				    char *history, 
				    Real *min_value, Real *max_value);

public void calc_gaussian_curvature(char *infilename, char *history,
				    Real min_value, Real max_value);


#define INTERNAL_X  2
#define INTERNAL_Y  1
#define INTERNAL_Z  0

char *prog_name;

double
  fwhm_3D[3],
  fwhm,
  standard;
int 
  verbose, 
  debug,
  clobber_flag, 
  apodize_data_flg,
  kernel_type,
  dimensions,
  do_Lvv_flag,
  do_gradient_flag,
  do_partials_flag;


static ArgvInfo argTable[] = {
  {"-fwhm", ARGV_FLOAT, (char *) 0, (char *) &fwhm, 
     "Full-width-half-maximum of gaussian kernel"},
  {"-standarddev", ARGV_FLOAT, (char *) 0, (char *) &standard,
     "Standard deviation of gaussian kernel"},
  {"-3dfwhm", ARGV_FLOAT, (char *) 3, (char *) fwhm_3D, 
     "Full-width-half-maximum of gaussian kernel"},
  {"-dimensions", ARGV_INT, (char *) 0, (char *) &dimensions,
     "Number of dimensions to blur (either 1,2 or 3)."},
  {NULL, ARGV_HELP, NULL, NULL,
     "Program flags."},
  {"-gaussian", ARGV_CONSTANT, (char *) KERN_GAUSSIAN, (char *) &kernel_type,
     "Use a gaussian smoothing kernel (default)."},
  {"-rect", ARGV_CONSTANT, (char *) KERN_RECT, (char *) &kernel_type,
     "Use a rect (box) smoothing kernel."},
  {"-gradient", ARGV_CONSTANT, (char *) TRUE, (char *) &do_gradient_flag, 
     "Create the gradient magnitude volume as well."},
  {"-Lvv", ARGV_CONSTANT, (char *) TRUE, (char *) &do_Lvv_flag, 
     "Create the Lvv volume as well."},
  {"-partial", ARGV_CONSTANT, (char *) TRUE, (char *) &do_partials_flag, 
     "Create the partial derivative and gradient magnitude volumes as well."},
  {"-no_apodize", ARGV_CONSTANT, (char *) FALSE, (char *) &apodize_data_flg, 
     "Do not apodize the data before blurring."},
  {"-no_clobber", ARGV_CONSTANT, (char *) FALSE, (char *) &clobber_flag,
     "Do not overwrite output file (default)."},
  {"-clobber", ARGV_CONSTANT, (char *) TRUE, (char *) &clobber_flag,
     "Overwrite output file."},
  {NULL, ARGV_HELP, NULL, NULL,
     "Options for logging progress. Default = -verbose."},
  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
     "Write messages indicating progress"},
  {"-quiet", ARGV_CONSTANT, (char *) FALSE , (char *) &verbose,
     "Do not write log messages"},
  {"-debug", ARGV_CONSTANT, (char *) TRUE, (char *) &debug,
     "Print out debug info."},
  {"-version", ARGV_FUNC, (char *) print_version_info, (char *)MNI_AUTOREG_LONG_VERSION,
     "Print out version info and exit."},
  {NULL, ARGV_END, NULL, NULL, NULL}
};

