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
@MODIFIED   : Revision 96.5  2009-07-23 22:34:00  claude
@MODIFIED   : cleanup in mincblur for VIO_X, VIO_Y, VIO_Z
@MODIFIED   :
@MODIFIED   : Revision 96.4  2009/02/13 04:14:40  rotor
@MODIFIED   :  * small updated to arguments of mincblur
@MODIFIED   :
@MODIFIED   : Revision 96.3  2006/11/28 09:12:21  rotor
@MODIFIED   :  * fixes to allow clean compile against minc 2.0
@MODIFIED   :
@MODIFIED   : Revision 96.2  2004/02/12 05:53:48  rotor
@MODIFIED   :  * removed public/private defs
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
 * Revision 1.11  1996/08/12  14:16:24  louis
 * Pre-release
 *
---------------------------------------------------------------------------- */

VIO_Status blur3D_volume(VIO_Volume data, int *xyzv,
                            double  kernel1, double  kernel2, double  kernel3, 
                            char *infile, 
                            char *outfile, 
                            FILE *reals_fp,
                            int kernel_type, char *history);

VIO_Status gradient3D_volume(FILE *ifd, 
                                VIO_Volume data, 
                                int *xyzv,
                                char *infile, 
                                char *outfile, 
                                int ndim,
                                char *history,
                                int curvature_flg);


void apodize_data(VIO_Volume data, int *xyzv,
                         double xramp1,double xramp2,
                         double yramp1,double yramp2,
                         double zramp1,double zramp2);

void calc_gradient_magnitude(char *infilename, 
                                    char *output_basename,
                                    char *history, 
                                    VIO_Real *min_value, VIO_Real *max_value);

void calc_gaussian_curvature(char *infilename, char *history,
                                    VIO_Real min_value, VIO_Real max_value);


char *prog_name;

double
  fwhm_3D[3],
  fwhm,
  standard;
int verbose;
int debug,
  clobber_flag, 
  apodize_data_flg,
  kernel_type,
  dimensions,
  do_gradient_flag,
  do_partials_flag;


ArgvInfo argTable[] = {
  {"-fwhm", ARGV_FLOAT, (char *) 1, (char *) &fwhm, 
     "Full-width-half-maximum of gaussian kernel"},
  {"-standarddev", ARGV_FLOAT, (char *) 1, (char *) &standard,
     "Standard deviation of gaussian kernel"},
  {"-3dfwhm", ARGV_FLOAT, (char *) 3, (char *) fwhm_3D, 
     "Full-width-half-maximum of gaussian kernel"},
  {"-dimensions", ARGV_INT, (char *) 0, (char *) &dimensions,
     "Number of dimensions to blur (either 1,2 or 3)."},
  
  {NULL, ARGV_HELP, NULL, NULL,
     "Program flags."},
  {"-no_clobber", ARGV_CONSTANT, (char *) FALSE, (char *) &clobber_flag,
     "Do not overwrite output file (default)."},
  {"-clobber", ARGV_CONSTANT, (char *) TRUE, (char *) &clobber_flag,
     "Overwrite output file."},
  {"-gaussian", ARGV_CONSTANT, (char *) KERN_GAUSSIAN, (char *) &kernel_type,
     "Use a gaussian smoothing kernel (default)."},
  {"-rect", ARGV_CONSTANT, (char *) KERN_RECT, (char *) &kernel_type,
     "Use a rect (box) smoothing kernel."},
  {"-gradient", ARGV_CONSTANT, (char *) TRUE, (char *) &do_gradient_flag, 
     "Create the gradient magnitude volume as well."},
  {"-partial", ARGV_CONSTANT, (char *) TRUE, (char *) &do_partials_flag, 
     "Create the partial derivative and gradient magnitude volumes as well."},
  {"-no_apodize", ARGV_CONSTANT, (char *) FALSE, (char *) &apodize_data_flg, 
     "Do not apodize the data before blurring."},
  
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














