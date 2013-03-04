#define RECTANGLE  0
#define ELLIPSE    1

#define X         0
#define Y         1
#define Z         2

char *prog_name;
int  debug        = FALSE;
int  verbose      = TRUE;
int  clobber_flag = FALSE;
int  partial_flag  = TRUE;

VIO_Real   voxel_range[2] = { 1.0,  255.0 };
VIO_Real   real_range[2]  = { 0.0,    1.0 };
VIO_Real   step[3]        = { 2.0, 2.0, 2.0 };
VIO_Real   start[3]       = { 0.0, 0.0, 0.0 };
int    count[3]       = { 128, 128, 128 };
int    datatype       = NC_BYTE;
int    is_signed      = FALSE;

int    object         = RECTANGLE;
VIO_Real   center[3]      = { 128.0, 128.0, 128.0 };
VIO_Real   width[3]       = { 60.0, 80.0, 30.0 };
VIO_Real   fill_value     = 1.0;
VIO_Real   edge_value     = 1.0;
VIO_Real   background     = 0.0;

static ArgvInfo argTable[] = {
  {NULL, ARGV_HELP, NULL, NULL,
     "Object definition options."},
  {"-rectangle", ARGV_CONSTANT, (char *) RECTANGLE, (char *) &object,
     "Build rectangular object (default)."},
  {"-ellipse", ARGV_CONSTANT, (char *) ELLIPSE, (char *) &object,
     "Build ellipsoid object."},
  {"-center", ARGV_FLOAT, (char *) 3, 
     (char *) center,
     "Center point of object (in world coords x y z)."},
  {"-xcenter", ARGV_FLOAT, (char *) 0, 
     (char *) &center[X],
     "Center of object in X dimension"},
  {"-ycenter", ARGV_FLOAT, (char *) 0, 
     (char *) &center[Y],
     "Center of object in Y dimension"},
  {"-zcenter", ARGV_FLOAT, (char *) 0, 
     (char *) &center[Z],
     "Center of object in Z dimension"},
  {"-width", ARGV_FLOAT, (char *) 3, 
     (char *) width,
     "Width  of object (in world mm along X Y Z)."},
  {"-xwidth", ARGV_FLOAT, (char *) 0, 
     (char *) &width[X],
     "Width of object in X dimension"},
  {"-ywidth", ARGV_FLOAT, (char *) 0, 
     (char *) &width[Y],
     "Width of object in Y dimension"},
  {"-zwidth", ARGV_FLOAT, (char *) 0, 
     (char *) &width[Z],
     "Width of object in Z dimension"},
  {"-fill_value", ARGV_FLOAT, (char *) 0, 
     (char *) &fill_value,
     "Real value used to fill object"},
  {"-edge_value", ARGV_FLOAT, (char *) 0, 
     (char *) &edge_value,
     "Real value used to fill object"},
  {"-background", ARGV_FLOAT, (char *) 0, 
     (char *) &background,
     "Real value used to fill background"},
  {"-no_partial", ARGV_CONSTANT, (char *) FALSE, (char *) &partial_flag,
     "Do not account for partial volume effects."},

  {NULL, ARGV_HELP, NULL, NULL,
     "Volume definition options."},
  {"-nelements", ARGV_INT, (char *) 3, 
     (char *) count,
     "Number of elements along each dimension (X, Y, Z)"},
  {"-xnelements", ARGV_INT, (char *) 0, 
     (char *) &count[X],
     "Number of elements along the X dimension"},
  {"-ynelements", ARGV_INT, (char *) 0, 
     (char *) &count[Y],
     "Number of elements along the Y dimension"},
  {"-znelements", ARGV_INT, (char *) 0, 
     (char *) &count[Z],
     "Number of elements along the Z dimension"},
  {"-step", ARGV_FLOAT, (char *) 3, 
     (char *) step,
     "Step size along each dimension (X, Y, Z)"},
  {"-xstep", ARGV_FLOAT, (char *) 0, 
     (char *) &step[X],
     "Step size along the X dimension"},
  {"-ystep", ARGV_FLOAT, (char *) 0, 
     (char *) &step[Y],
     "Step size along the Y dimension"},
  {"-zstep", ARGV_FLOAT, (char *) 0, 
     (char *) &step[Z],
     "Step size along the Z dimension"},
  {"-start", ARGV_FLOAT, (char *) 3, 
     (char *) start,
     "Start point along each dimension (X, Y, Z)"},
  {"-xstart", ARGV_FLOAT, (char *) 0, 
     (char *) &start[X],
     "Start point along the X dimension"},
  {"-ystart", ARGV_FLOAT, (char *) 0, 
     (char *) &start[Y],
     "Start point along the Y dimension"},
  {"-zstart", ARGV_FLOAT, (char *) 0, 
     (char *) &start[Z],
     "Start point along the Z dimension"},
  {"-byte", ARGV_CONSTANT, (char *) NC_BYTE, (char *) &datatype,
     "Write out byte data (default)"},
  {"-short", ARGV_CONSTANT, (char *) NC_SHORT, (char *) &datatype,
     "Write out short integer data"},
  {"-long", ARGV_CONSTANT, (char *) NC_LONG, (char *) &datatype,
     "Write out long integer data"},
  {"-float", ARGV_CONSTANT, (char *) NC_FLOAT, (char *) &datatype,
     "Write out single-precision floating-point data"},
  {"-double", ARGV_CONSTANT, (char *) NC_DOUBLE, (char *) &datatype,
     "Write out double-precision floating-point data"},
  {"-signed", ARGV_CONSTANT, (char *) TRUE, (char *) &is_signed,
     "Write signed integer data"},
  {"-unsigned", ARGV_CONSTANT, (char *) FALSE, (char *) &is_signed,
     "Write unsigned integer data (default)"},
  {"-voxel_range", ARGV_FLOAT, (char *) 2, (char *) voxel_range,
     "Valid voxel range for output data"},
  {"-real_range", ARGV_FLOAT, (char *) 2, (char *) real_range,
     "Valid real range for output data"},
  
  
  
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
