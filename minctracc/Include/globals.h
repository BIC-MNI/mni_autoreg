/*  ------------------------ global variables  ------------------------ */

char  *prog_name     = NULL;

Volume  mask_data  = NULL;
Volume  mask_model = NULL;
Volume  model      = NULL;
Volume  data       = NULL;


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
     (char *) &main_args.trans_info.transformation,
     "Initial world transformation. (Default = identity)."},
  {"-est_center", ARGV_CONSTANT, (char *) TRUE, (char *) &main_args.trans_flags.estimate_center,
     "use center estimated from Principal axis trans."},
  {"-est_scales", ARGV_CONSTANT, (char *) TRUE, (char *) &main_args.trans_flags.estimate_scale,
     "use scales estimated from Principal axis trans."},
  {"-est_rotations", ARGV_CONSTANT, (char *) TRUE, (char *) &main_args.trans_flags.estimate_rots,
     "use rotations estimated from Principal axis trans."},
  {"-est_translations", ARGV_CONSTANT, (char *) TRUE, (char *) &main_args.trans_flags.estimate_trans,
     "use translations estimated from Principal axis trans."},
  {"-center", ARGV_FLOAT, (char *) 3, 
     (char *) main_args.trans_info.center,
     "Force center of rotation and scale."},
  {NULL, ARGV_HELP, NULL, NULL,
     "Output transformation type. Default = -procrustes."},
  {"-pat", ARGV_CONSTANT, (char *) TRANS_PAT, (char *) &main_args.trans_info.transform_type,
     "principal axis transformation matrix (input matrix ignored)."},
  {"-lsq6", ARGV_CONSTANT, (char *) TRANS_LSQ6, (char *) &main_args.trans_info.transform_type,
     "6 parameter transformation (3 translation, 3 rotation, scale=1.0)."},
  {"-lsq7", ARGV_CONSTANT, (char *) TRANS_PROCRUSTES, (char *) &main_args.trans_info.transform_type,
     "7 parameter transformation (lsq6 + 1 global scale, same as -procrustes)."},
  {"-lsq9", ARGV_CONSTANT, (char *) TRANS_LSQ9, (char *) &main_args.trans_info.transform_type,
     "9 parameter transformation (lsq6 + 3 scales)."},
  {"-lsq12", ARGV_CONSTANT, (char *) TRANS_LSQ, (char *) &main_args.trans_info.transform_type,
     "full 12 parameter transformation (lsq9 + 3 shears)."},
  {"-lsq", ARGV_CONSTANT, (char *) TRANS_LSQ, (char *) &main_args.trans_info.transform_type,
     "full 12 parameter transformation (same as -lsq12)."},
  {"-procrustes", ARGV_CONSTANT, (char *) TRANS_PROCRUSTES, (char *) &main_args.trans_info.transform_type,
     "procrustes transformation (3 trans, 3 rots, 1 scale), same as -lsq7."},
  {"-forward", ARGV_CONSTANT, (char *) FALSE,
     (char *) &main_args.trans_info.invert_mapping_flag,
     "Recover forward transformation (source -> model, default).\n"},
  {"-invert", ARGV_CONSTANT, (char *) TRUE,
     (char *) &main_args.trans_info.invert_mapping_flag,
     "Recover inverted transformation (model -> source).\n"},
  
  {NULL, ARGV_HELP, NULL, NULL,
     "Options for mask volumes."},
  {"-model_mask", ARGV_FUNC, (char *) get_mask_file, 
     (char *) &main_args.filenames.mask_model,
     "Specifies a binary mask file for the target."},
  {"-source_mask", ARGV_FUNC, (char *) get_mask_file, 
     (char *) &main_args.filenames.mask_data,
     "Specifies a binary mask file for the source."},
  
  {NULL, ARGV_HELP, NULL, NULL,
     "Interpolation options."},
  {"-trilinear", ARGV_CONSTANT, (char *) trilinear_interpolant, 
     (char *) &main_args.interpolant,
     "Do trilinear interpolation"},
  {"-tricubic", ARGV_CONSTANT, (char *) tricubic_interpolant, 
     (char *) &main_args.interpolant,
     "Do tricubic interpolation"},
  {"-nearest_neighbour", ARGV_CONSTANT, 
     (char *) nearest_neighbour_interpolant, (char *) &main_args.interpolant,
     "Do nearest neighbour interpolation"},
  
  {NULL, ARGV_HELP, NULL, NULL,
     "Optimization objective functions. (default = -xcorr)"},
  {"-xcorr", ARGV_CONSTANT, (char *) xcorr_objective, 
     (char *) &main_args.obj_function,
     "Use cross correlation."},
  {"-zscore", ARGV_CONSTANT, (char *) zscore_objective, 
     (char *) &main_args.obj_function,
     "Use normalized difference."},
  {"-ssc", ARGV_CONSTANT, (char *) ssc_objective, 
     (char *) &main_args.obj_function,
     "Use stocastic sign change."},
  {"-vr", ARGV_CONSTANT,      (char *) vr_objective, 
     (char *) &main_args.obj_function,
     "Use variance of ratio vol1/vol2"},
  {"-groups", ARGV_INT, (char *) 0, 
     (char *) &main_args.groups,
     "Number of groups for ratio calculations."},
  {"-threshold", ARGV_FLOAT, (char *) 0, 
     (char *) &main_args.threshold,
     "Lower limit for voxel threshold"},
  {"-speckle", ARGV_FLOAT, (char *) 0, 
     (char *) &main_args.speckle,
     "percent speckle noise to add to source"},
  
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
     (char *) main_args.step,
     "Step size along each dimension (X, Y, Z)"},
  {"-xstep", ARGV_FLOAT, (char *) 0, 
     (char *) &main_args.step[X],
     "Step size along the X dimension"},
  {"-ystep", ARGV_FLOAT, (char *) 0, 
     (char *) &main_args.step[Y],
     "Step size along the Y dimension"},
  {"-zstep", ARGV_FLOAT, (char *) 0, 
     (char *) &main_args.step[Z],
     "Step size along the Z dimension"},
  
  {NULL, ARGV_HELP, NULL, NULL,
     "Options for logging progress. Default = -verbose 1."},
  {"-verbose", ARGV_INT, (char *) 0, (char *) &main_args.flags.verbose,
     "Write messages indicating progress"},
  {"-quiet", ARGV_CONSTANT, (char *) 0 , (char *) &main_args.flags.verbose,
     "Do not write log messages"},
  {"-debug", ARGV_CONSTANT, (char *) TRUE, (char *) &main_args.flags.debug,
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

