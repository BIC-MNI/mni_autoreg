char *prog_name;
int  debug;
int  verbose;
VIO_Real threshold = 0.0;
int  mincresample = FALSE;
int  mincreshape  = FALSE;
int  minccrop     = FALSE;
int  two_lines    = FALSE; 

static ArgvInfo argTable[] = {
  {"-threshold", ARGV_FLOAT, (char *) FALSE, (char *) &threshold,
     "Real value threshold for bounding box."},
  {"-one_line", ARGV_CONSTANT, (char *) FALSE, (char *) &two_lines,
     "Output on one line (default): start_x y z width_x y z"},
  {"-two_lines", ARGV_CONSTANT, (char *) TRUE, (char *) &two_lines,
     "Output on two lines: start_x y z \\n width_x y z"},
  {"-mincresample", ARGV_CONSTANT, (char *) TRUE, (char *) &mincresample,
     "Output format for mincresample: (-step x y z -start x y z -nelements x y z"},
  {"-mincreshape", ARGV_CONSTANT, (char *) TRUE, (char *) &mincreshape,
     "Output format for mincreshape: (-start x,y,z -count dx,dy,dz"},
  {"-minccrop", ARGV_CONSTANT, (char *) TRUE, (char *) &minccrop,
     "Output format for minccrop: (-xlim x1 x2 -ylim y1 y2 -zlim z1 z2"},
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

