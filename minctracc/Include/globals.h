/*  ------------------------ global variables  ------------------------ */

char  *prog_name     = NULL;

volume_struct  *mask_data  = NULL;
volume_struct  *mask_model = NULL;
volume_struct  *model      = NULL;
volume_struct  *data       = NULL;


double  ftol         =  0.0001;
double  simplex_size = 20.0;

int  invert_mapping_flag = FALSE;

  
/*------------------------    Command line arguments --------------------*/

static ArgvInfo argTable[] = {
  {NULL, ARGV_HELP, NULL, NULL,
     "---Transformation maps one to volume two---"},
  {NULL, ARGV_HELP, NULL, NULL,
	"Initial transformation information."},
  {"-transformation", ARGV_FUNC, (char *) get_transformation, 
     (char *) &args.trans_info.transformation,
     "File giving initial world transformation. (Default = identity)."},
  
  {NULL, ARGV_HELP, NULL, NULL,
     "Output transformation type. Default = -procrustes."},
  {"-lsq6", ARGV_CONSTANT, (char *) TRANS_LSQ6, (char *) &args.trans_info.transform_type,
     "return 6 parameter transformation (3 translation, 3 rotation, scale=1.0)."},
  {"-lsq7", ARGV_CONSTANT, (char *) TRANS_PROCRUSTES, (char *) &args.trans_info.transform_type,
     "return 7 parameter transformation (lsq6 + 1 global scale, same as -procrustes)."},
  {"-lsq9", ARGV_CONSTANT, (char *) TRANS_LSQ9, (char *) &args.trans_info.transform_type,
     "return 9 parameter transformation (lsq6 + 3 scales)."},
  {"-lsq12", ARGV_CONSTANT, (char *) TRANS_LSQ, (char *) &args.trans_info.transform_type,
     "return full 12 parameter transformation (lsq9 + 3 shears)."},
  {"-lsq", ARGV_CONSTANT, (char *) TRANS_LSQ, (char *) &args.trans_info.transform_type,
     "return full 12 parameter transformation (same as -lsq12)."},
  {"-procrustes", ARGV_CONSTANT, (char *) TRANS_PROCRUSTES, (char *) &args.trans_info.transform_type,
     "return procrustes transformation (3 trans, 3 rots, 1 scale), same as -lsq7."},
  {"-forward", ARGV_CONSTANT, (char *) FALSE,
     (char *) &args.trans_info.invert_mapping_flag,
     "Recover forward transformation (source -> model, default).\n"},
  {"-invert", ARGV_CONSTANT, (char *) TRUE,
     (char *) &args.trans_info.invert_mapping_flag,
     "Recover inverted transformation (model -> source).\n"},
  
  {NULL, ARGV_HELP, NULL, NULL,
     "Options for mask volumes."},
  {"-model_mask", ARGV_FUNC, (char *) get_mask_file, 
     (char *) &args.filenames.mask_model,
     "Specifies a binary mask file for the target."},
  {"-source_mask", ARGV_FUNC, (char *) get_mask_file, 
     (char *) &args.filenames.mask_data,
     "Specifies a binary mask file for the source."},
  
  {NULL, ARGV_HELP, NULL, NULL,
     "Interpolation options."},
  {"-trilinear", ARGV_CONSTANT, (char *) trilinear_interpolant, 
     (char *) &args.interpolant,
     "Do trilinear interpolation"},
  {"-tricubic", ARGV_CONSTANT, (char *) tricubic_interpolant, 
     (char *) &args.interpolant,
     "Do tricubic interpolation"},
  {"-nearest_neighbour", ARGV_CONSTANT, 
     (char *) nearest_neighbour_interpolant, (char *) &args.interpolant,
     "Do nearest neighbour interpolation"},
  
  {NULL, ARGV_HELP, NULL, NULL,
     "Options for optimization."},
  {"-tol", ARGV_FLOAT, (char *) 0, 
     (char *) &ftol,
     "Stopping criteria tolerence"},
  {"-simplex", ARGV_FLOAT, (char *) 0, 
     (char *) &simplex_size,
     "Radius of simplex volume."},

  {NULL, ARGV_HELP, NULL, NULL,
     "Options for 3D lattice (default = target)."},
  {"-step", ARGV_FLOAT, (char *) 3, 
     (char *) args.step,
     "Step size along each dimension (X, Y, Z)"},
  {"-xstep", ARGV_FLOAT, (char *) 0, 
     (char *) &args.step[X],
     "Step size along the X dimension"},
  {"-ystep", ARGV_FLOAT, (char *) 0, 
     (char *) &args.step[Y],
     "Step size along the Y dimension"},
  {"-zstep", ARGV_FLOAT, (char *) 0, 
     (char *) &args.step[Z],
     "Step size along the Z dimension"},
  
  {NULL, ARGV_HELP, NULL, NULL,
     "Options for logging progress. Default = -verbose 1."},
  {"-verbose", ARGV_INT, (char *) 0, (char *) &args.flags.verbose,
     "Write messages indicating progress"},
  {"-quiet", ARGV_CONSTANT, (char *) 0 , (char *) &args.flags.verbose,
     "Do not write log messages"},
  {"-debug", ARGV_CONSTANT, (char *) TRUE, (char *) &args.flags.debug,
     "Print out debug info."},
  {NULL, ARGV_END, NULL, NULL, NULL}
};


static Linear_Transformation identity_matrix = {
  1.0, 0.0, 0.0, 0.0,
  0.0, 1.0, 0.0, 0.0,
  0.0, 0.0, 1.0, 0.0
};

static Transformation identity_transformation = {
  TRUE, do_linear_transformation, &identity_matrix
};

