/*  ------------------------ global variables  ------------------------ */

char  *prog_name     = NULL;

Volume  mask_data  = NULL;
Volume  mask_model = NULL;
Volume  model      = NULL;
Volume  model_dx   = NULL;
Volume  model_dy   = NULL;
Volume  model_dz   = NULL;
Volume  model_dxyz = NULL;
Volume  data       = NULL;
Volume  data_dx    = NULL;
Volume  data_dy    = NULL;
Volume  data_dz    = NULL;
Volume  data_dxyz  = NULL;


double  ftol         =  0.005;
double  simplex_size = 20.0;

int  invert_mapping_flag = FALSE;
int clobber_flag = FALSE;

  
/*------------------------    Command line arguments --------------------*/

static ArgvInfo argTable[] = {
  {NULL, ARGV_HELP, NULL, NULL,
     "---Transformation maps one to volume two---"},
  {NULL, ARGV_HELP, NULL, NULL,
	"Initial transformation information."},
  {"-no_clobber", ARGV_CONSTANT, (char *) FALSE, (char *) &clobber_flag,
     "Do not overwrite output file (default)."},
  {"-clobber", ARGV_CONSTANT, (char *) TRUE, (char *) &clobber_flag,
     "Overwrite output file."},
  {"-transformation", ARGV_FUNC, (char *) get_transformation, 
     (char *) &main_args.trans_info,
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
  {"-lsq3", ARGV_CONSTANT, (char *) TRANS_LSQ3, (char *) &main_args.trans_info.transform_type,
     "3 parameter transformation (3 translations only)."},
  {"-lsq6", ARGV_CONSTANT, (char *) TRANS_LSQ6, (char *) &main_args.trans_info.transform_type,
     "6 parameter transformation (3 translation, 3 rotation, scale=1.0)."},
  {"-lsq7", ARGV_CONSTANT, (char *) TRANS_LSQ7, (char *) &main_args.trans_info.transform_type,
     "7 parameter transformation (lsq6 + 1 global scale, same as -procrustes)."},
  {"-lsq9", ARGV_CONSTANT, (char *) TRANS_LSQ9, (char *) &main_args.trans_info.transform_type,
     "9 parameter transformation (lsq6 + 3 scales)."},
  {"-lsq10", ARGV_CONSTANT, (char *) TRANS_LSQ10, (char *) &main_args.trans_info.transform_type,
     "10 parameter transformation (lsq9 + 1 shear)."},
  {"-lsq12", ARGV_CONSTANT, (char *) TRANS_LSQ12, (char *) &main_args.trans_info.transform_type,
     "full 12 parameter transformation (lsq9 + 3 shears)."},
  {"-lsq", ARGV_CONSTANT, (char *) TRANS_LSQ, (char *) &main_args.trans_info.transform_type,
     "full 12 parameter transformation (same as -lsq12)."},
  {"-procrustes", ARGV_CONSTANT, (char *) TRANS_LSQ7, (char *) &main_args.trans_info.transform_type,
     "procrustes transformation (3 trans, 3 rots, 1 scale), same as -lsq7."},

  {"-nonlinear", ARGV_CONSTANT, (char *) TRANS_NONLIN, (char *) &main_args.trans_info.transform_type,
     "recover nonlinear deformation field."}, 

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
     "Use stochastic sign change."},
  {"-vr", ARGV_CONSTANT,      (char *) vr_objective, 
     (char *) &main_args.obj_function,
     "Use variance of ratio vol1/vol2"},
  {"-groups", ARGV_INT, (char *) 0, 
     (char *) &main_args.groups,
     "Number of groups for ratio calculations."},
  {"-threshold", ARGV_FLOAT, (char *) 2, 
     (char *) main_args.threshold,
     "Lower limit for voxel threshold"},
  {"-speckle", ARGV_FLOAT, (char *) 0, 
     (char *) &main_args.speckle,
     "percent speckle noise to add to source"},
  
  {NULL, ARGV_HELP, NULL, NULL,
     "Options for optimization."},
  {"-tol", ARGV_FLOAT, (char *) 0, 
     (char *) &ftol,
     "Stopping criteria tolerance"},
  {"-simplex", ARGV_FLOAT, (char *) 0, 
     (char *) &simplex_size,
     "Radius of simplex volume."},
  {"-w_translations", ARGV_FLOAT, (char *) 3, 
     (char *) &main_args.trans_info.weights[0],
     "Optimization weight of translation in x, y, z."},
  {"-w_rotations", ARGV_FLOAT, (char *) 3, 
     (char *) &main_args.trans_info.weights[3],
     "Optimization weight of rotations around x, y, z."},
  {"-w_scales", ARGV_FLOAT, (char *) 3, 
     (char *) &main_args.trans_info.weights[6],
     "Optimization weight of scaling along x, y, z."},
  {"-w_shear", ARGV_FLOAT, (char *) 3, 
     (char *) &main_args.trans_info.weights[9],
     "Optimization weight of shears a,b and c."},

  {NULL, ARGV_HELP, NULL, NULL,
     "Options for measurement comparison."},
  {"-matlab", ARGV_STRING, (char *) 0, 
     (char *) &main_args.filenames.matlab_file,
     "Output selected objective function value curves."},
  {"-measure", ARGV_STRING, (char *) 0, 
     (char *) &main_args.filenames.measure_file,
     "Output value of each obj. func. for given x-form."},

  {NULL, ARGV_HELP, NULL, NULL,
     "Options for 3D lattice (default = target)."},
  {"-step", ARGV_FLOAT, (char *) 3, 
     (char *) main_args.step,
     "Step size along each dimension (X, Y, Z)"},
  {"-xstep", ARGV_FLOAT, (char *) 0, 
     (char *) &main_args.step[0],
     "Step size along the column dimension"},
  {"-ystep", ARGV_FLOAT, (char *) 0, 
     (char *) &main_args.step[1],
     "Step size along the row dimension"},
  {"-zstep", ARGV_FLOAT, (char *) 0, 
     (char *) &main_args.step[2],
     "Step size along the slice dimension"},
  
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


