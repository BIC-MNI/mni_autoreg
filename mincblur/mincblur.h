/* ----------------------------- MNI Header -----------------------------------
@NAME       : mincblur.h
@DESCRIPTION: header file for mincblur.c
@CREATED    : Wed Jun 23 09:04:34 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */

public Status blur3D_volume(Volume data,
			    float kernel1, 
			    char *infile, 
			    char *outfile, 
			    int ndim, int kernel_type, char *history);

public Status gradient3D_volume(FILE *ifd, 
				Volume data, 
				char *infile, 
				char *outfile, 
				int ndim,
				char *history,
				int curvature_flg);

public void apodize_data(Volume data, 
			 double xramp1,double xramp2,
			 double yramp1,double yramp2,
			 double zramp1,double zramp2);

public void calc_gradient_magnitude(char *infilename, char *history, 
				    Real *min_value, Real *max_value);

public void calc_gaussian_curvature(char *infilename, char *history,
				    Real min_value, Real max_value);


#define INTERNAL_X  2
#define INTERNAL_Y  1
#define INTERNAL_Z  0

char *prog_name;
int 
  verbose,
  debug,
  apodize_data_flg,
  kernel_type,
  dimensions,
  gradonlyflg,
  curvatureflg,
  curveonlyflg,
  gradonlyflg,
  bluronlyflg;
double   
  fwhm,
  standard;
int clobber_flag = FALSE;

extern int ms_volume_reals_flag;


static ArgvInfo argTable[] = {
  {"-fwhm", ARGV_FLOAT, (char *) 0, (char *) &fwhm, 
     "Full-width-half-maximum of gaussian kernel"},
  {"-standarddev", ARGV_FLOAT, (char *) 0, (char *) &standard,
     "Standard deviation of gaussian kernel"},
  {"-dimensions", ARGV_INT, (char *) 0, (char *) &dimensions,
     "Number of dimensions to blur (either 1,2 or 3)."},
  {NULL, ARGV_HELP, NULL, NULL,
     "Program flags."},
  {"-gaussian", ARGV_CONSTANT, (char *) KERN_GAUSSIAN, (char *) &kernel_type,
     "Use a gaussian smoothing kernel (default)."},
  {"-rect", ARGV_CONSTANT, (char *) KERN_RECT, (char *) &kernel_type,
     "Use a rect (box) smoothing kernel."},
  {"-blur_only", ARGV_CONSTANT, (char *) TRUE, (char *) &bluronlyflg,
     "Create only the blurred volume."},
  {"-grad_only", ARGV_CONSTANT, (char *) TRUE, (char *) &gradonlyflg, 
     "Create only the gradient volumes."},
  {"-gcur_only", ARGV_CONSTANT, (char *) TRUE, (char *) &curveonlyflg, 
     "Create only the curvature volume."},
  {"-curvature", ARGV_CONSTANT, (char *) TRUE, (char *) &curvatureflg, 
     "Create data to calculate curvature volumes."},
  {"-no_apodize", ARGV_CONSTANT, (char *) FALSE, (char *) &apodize_data_flg, 
     "Do not apodize the data before blurring."},
  {"-no_reals", ARGV_CONSTANT, (char *) FALSE, (char *) &ms_volume_reals_flag, 
     "Do not write out the real (float) data."},
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
  {NULL, ARGV_END, NULL, NULL, NULL}
};

