#include <config.h>
#include <float.h>
#include <ParseArgv.h>
#include <volume_io.h>
#include <minctracc.h>
#include <objectives.h>
#include "local_macros.h"


/* Globals !! :( */

char  *prog_name                 = NULL;

VIO_Volume  mask_data                = NULL;
VIO_Volume  mask_model               = NULL;
VIO_Volume  model                    = NULL;
VIO_Volume  model_dx                 = NULL;
VIO_Volume  model_dy                 = NULL;
VIO_Volume  model_dz                 = NULL;
VIO_Volume  model_dxyz               = NULL;
VIO_Volume  data                     = NULL;
VIO_Volume  data_dx                  = NULL;
VIO_Volume  data_dy                  = NULL;
VIO_Volume  data_dz                  = NULL;
VIO_Volume  data_dxyz                = NULL;

double  ftol                     = 0.005;
double  simplex_size             = 20.0;
int     iteration_limit          = 4;
double  iteration_weight         = 0.6;
double  smoothing_weight         = 0.5;
double  similarity_cost_ratio    = 0.5;
int     number_dimensions        = 3;
int     Matlab_num_steps         = 15;
int     Diameter_of_local_lattice= 5;

int     invert_mapping_flag      = FALSE;
int     clobber_flag             = FALSE;

VIO_Real initial_corr, final_corr;

Arg_Data main_argsX;

ArgvInfo argTable[] = {
  {NULL, ARGV_HELP, NULL, NULL,
     "---Transformation maps one to volume two---"},
  {NULL, ARGV_HELP, NULL, NULL,
        "Initial transformation information."},
  {"-no_clobber", ARGV_CONSTANT, (char *) FALSE, (char *) &clobber_flag,
     "Do not overwrite output file (default)."},
  {"-clobber", ARGV_CONSTANT, (char *) TRUE, (char *) &clobber_flag,
     "Overwrite output file."},
  {"-transformation", ARGV_FUNC, (char *) get_transformation, 
     (char *) &main_argsX.trans_info,
     "Initial world transformation. (Default = PAT)."},
  {"-identity", ARGV_CONSTANT, (char *) TRUE,
     (char *) &main_argsX.trans_info.use_identity,
     "Use identity transformation for starting point.\n"},
  {"-est_center", ARGV_CONSTANT, (char *) TRUE, (char *) &main_argsX.trans_flags.estimate_center,
     "use center estimated from Principal axis trans."},
  {"-est_scales", ARGV_CONSTANT, (char *) TRUE, (char *) &main_argsX.trans_flags.estimate_scale,
     "use scales estimated from Principal axis trans."},
/*  {"-est_rotations", ARGV_CONSTANT, (char *) TRUE, (char *) &main_argsX.trans_flags.estimate_rots,
     "use rotations estimated from Principal axis trans."},*/
  {"-est_translations", ARGV_CONSTANT, (char *) TRUE, (char *) &main_argsX.trans_flags.estimate_trans,
     "use translations estimated from Principal axis trans."},
  {"-center", ARGV_FLOAT, (char *) 3, 
     (char *) main_argsX.trans_info.center,
     "Force center of rotation and scale."},
  {NULL, ARGV_HELP, NULL, NULL,
     "\nOutput transformation type. Default = -procrustes."},
  {"-pat", ARGV_CONSTANT, (char *) TRANS_PAT, (char *) &main_argsX.trans_info.transform_type,
     "principal axis transformation matrix (input matrix ignored)."},
  {"-lsq3", ARGV_CONSTANT, (char *) TRANS_LSQ3, (char *) &main_argsX.trans_info.transform_type,
     "3 parameter transformation (3 translations only)."},
  {"-lsq6", ARGV_CONSTANT, (char *) TRANS_LSQ6, (char *) &main_argsX.trans_info.transform_type,
     "6 parameter transformation (3 translation, 3 rotation, scale=1.0)."},
  {"-lsq7", ARGV_CONSTANT, (char *) TRANS_LSQ7, (char *) &main_argsX.trans_info.transform_type,
     "7 parameter transformation (lsq6 + 1 global scale, same as -procrustes)."},
  {"-lsq9", ARGV_CONSTANT, (char *) TRANS_LSQ9, (char *) &main_argsX.trans_info.transform_type,
     "9 parameter transformation (lsq6 + 3 scales)."},
  {"-lsq10", ARGV_CONSTANT, (char *) TRANS_LSQ10, (char *) &main_argsX.trans_info.transform_type,
     "10 parameter transformation (lsq9 + 1 shear)."},
  {"-lsq12", ARGV_CONSTANT, (char *) TRANS_LSQ12, (char *) &main_argsX.trans_info.transform_type,
     "full 12 parameter transformation (lsq9 + 3 shears)."},
  {"-lsq", ARGV_CONSTANT, (char *) TRANS_LSQ, (char *) &main_argsX.trans_info.transform_type,
     "full 12 parameter transformation (same as -lsq12)."},
  {"-procrustes", ARGV_CONSTANT, (char *) TRANS_LSQ7, (char *) &main_argsX.trans_info.transform_type,
     "procrustes transformation (3 trans, 3 rots, 1 scale), same as -lsq7."},
  {"-forward", ARGV_CONSTANT, (char *) FALSE,
     (char *) &main_argsX.trans_info.invert_mapping_flag,
     "Recover forward transformation (source -> model, default).\n"},
  {"-invert", ARGV_CONSTANT, (char *) TRUE,
     (char *) &main_argsX.trans_info.invert_mapping_flag,
     "Recover inverted transformation (model -> source).\n"},
  {"-quaternions", ARGV_CONSTANT, (char *) TRANS_QUAT, (char *) &main_argsX.trans_info.rotation_type,
     "rotation with quaternion, in output give a vector and an angle.\n"},
  {"-rotations", ARGV_CONSTANT, (char *) TRANS_ROT, (char *) &main_argsX.trans_info.rotation_type,
     "rotation without quaternion. default type.\n"},

  
  {NULL, ARGV_HELP, NULL, NULL,
     "\nOptions for feature volumes."},
  {"-feature_vol", ARGV_GENFUNC, (char *) get_feature_volumes, 
     (char *) &main_argsX.filenames.mask_model,
     "Specify extra features <f1.mnc> <f2.mnc> [xcorr|diff|sqdiff|label|chamfer|corrcoeff|opticalflow] [weight]."},

  {NULL, ARGV_HELP, NULL, NULL,
     "\nOptions for mask volumes."},
  {"-model_mask", ARGV_FUNC, (char *) get_mask_file, 
     (char *) &main_argsX.filenames.mask_model,
     "Specifies a binary mask file for the target."},
  {"-source_mask", ARGV_FUNC, (char *) get_mask_file, 
     (char *) &main_argsX.filenames.mask_data,
     "Specifies a binary mask file for the source."},
  
  {NULL, ARGV_HELP, NULL, NULL,
     "\nInterpolation options. (Default = -trilinear)"},
  {"-trilinear", ARGV_CONSTANT, (char *) TRILINEAR,
     (char *) &main_argsX.interpolant_type,
     "Do trilinear interpolation"},
  {"-tricubic", ARGV_CONSTANT, (char *) TRICUBIC,
     (char *) &main_argsX.interpolant_type,
     "Do tricubic interpolation"},
  {"-nearest_neighbour", ARGV_CONSTANT, (char *) N_NEIGHBOUR,
     (char *) &main_argsX.interpolant_type,
     "Do nearest neighbour interpolation"},
  
  {NULL, ARGV_HELP, NULL, NULL,
     "\nLinear optimization objective functions. (default = -xcorr)"},
  {"-xcorr", ARGV_CONSTANT, (char *) XCORR,
     (char *) &main_argsX.obj_function_type,
     "Use cross correlation."},
  {"-zscore", ARGV_CONSTANT, (char *) ZSCORE,
     (char *) &main_argsX.obj_function_type,
     "Use normalized difference."},
  {"-ssc", ARGV_CONSTANT, (char *) SSC,
     (char *) &main_argsX.obj_function_type,
     "Use stochastic sign change (Minoshima)."},
  {"-vr", ARGV_CONSTANT,      (char *) VR,
     (char *) &main_argsX.obj_function_type,
     "Use variance of ratio vol1/vol2 (Woods)."},
  {"-mi", ARGV_CONSTANT,      (char *) MUTUAL_INFORMATION,
     (char *) &main_argsX.obj_function_type,
     "Use mutual information similarity measure (Collignon)"},
  {"-nmi", ARGV_CONSTANT,      (char *) NORMALIZED_MUTUAL_INFORMATION,
     (char *) &main_argsX.obj_function_type,
     "Use mutual information similarity measure (Studholme)"},
  {"-groups", ARGV_INT, (char *) 0, 
     (char *) &main_argsX.groups,
     "Number of groups for -vr and -mi."},
  {"-blur_pdf", ARGV_INT, (char *) 0, 
     (char *) &main_argsX.blur_pdf,
     "Size of pdf and jpdf blurring kernel for -mi."},
  {"-threshold", ARGV_FLOAT, (char *) 2, 
     (char *) main_argsX.threshold,
     "Lower limit for voxel real value threshold"},
  {"-speckle", ARGV_FLOAT, (char *) 0, 
     (char *) &main_argsX.speckle,
     "percent speckle noise to add to source"},
  
  {NULL, ARGV_HELP, NULL, NULL,
     "\nOptions for linear optimization."},
  {"-tol", ARGV_FLOAT, (char *) 0, 
     (char *) &ftol,
     "Stopping criteria tolerance"},
  {"-simplex", ARGV_FLOAT, (char *) 0, 
     (char *) &simplex_size,
     "Radius of simplex volume."},
  {"-w_translations", ARGV_FLOAT, (char *) 3, 
     (char *) &main_argsX.trans_info.weights[0],
     "Optimization weight of translation in x, y, z."},
  {"-w_rotations", ARGV_FLOAT, (char *) 3, 
     (char *) &main_argsX.trans_info.weights[3],
     "Optimization weight of rotations around x, y, z."},
  {"-w_scales", ARGV_FLOAT, (char *) 3, 
     (char *) &main_argsX.trans_info.weights[6],
     "Optimization weight of scaling along x, y, z."},
  {"-w_shear", ARGV_FLOAT, (char *) 3, 
     (char *) &main_argsX.trans_info.weights[9],
     "Optimization weight of shears a,b and c."},
  {"-use_bfgs", ARGV_CONSTANT, (char *) TRUE, (char *) &main_argsX.trans_info.use_bfgs,
     "use BFGS optimizer instead of amoeba "},

  {NULL, ARGV_HELP, NULL, NULL,
     "\nOptions for measurement comparison."},
  {"-matlab", ARGV_STRING, (char *) 0, 
     (char *) &main_argsX.filenames.matlab_file,
     "Output curves for selected objective function vs parameter."},
  {"-num_steps", ARGV_INT, (char *) 0, (char *) &Matlab_num_steps,
     "Number of steps at which to measure obj fn for matlab output."},
  {"-measure", ARGV_STRING, (char *) 0, 
     (char *) &main_argsX.filenames.measure_file,
     "Output value of each obj. func. for given x-form."},

  {NULL, ARGV_HELP, NULL, NULL,
     "\nOptions for 3D lattice."},
  {"-source_lattice", ARGV_CONSTANT, (char *) 1,
     (char *) &main_argsX.force_lattice,
     "Force lattice on source volume (def=auto selected).\n"},
  {"-model_lattice", ARGV_CONSTANT, (char *) 2,
     (char *) &main_argsX.force_lattice,
     "Force lattice on model volume (def=auto selected).\n"},
  {"-step", ARGV_FLOAT, (char *) 3, 
     (char *) main_argsX.step,
     "Step size along each dimension (X, Y, Z)"},
  {"-xstep", ARGV_FLOAT, (char *) 0, 
     (char *) &main_argsX.step[0],
     "Step size along the column dimension (mm)"},
  {"-ystep", ARGV_FLOAT, (char *) 0, 
     (char *) &main_argsX.step[1],
     "Step size along the row dimension"},
  {"-zstep", ARGV_FLOAT, (char *) 0, 
     (char *) &main_argsX.step[2],
     "Step size along the slice dimension"},

  {NULL, ARGV_HELP, NULL, NULL,
     "\nNon-linear transformation information:"},
  {"-nonlinear", ARGV_FUNC, (char*)get_nonlinear_objective, NULL,
      "recover nonlinear deformation field.  Optional arg {xcorr|diff|sqdiff|label|chamfer|corrcoeff|opticalflow} sets objective function."},
/*   {"-2D-non-lin", ARGV_CONSTANT, (char *) 2, (char *) &number_dimensions, */
/*      "Estimate the non-lin fit on a 2D slice only."}, */
/*   {"-3D-non-lin", ARGV_CONSTANT, (char *) 3, (char *) &number_dimensions, */
/*      "Estimate the non-lin fit on a 3D volume (default)."}, */
  {"-sub_lattice", ARGV_INT, (char *) 0, (char *) &Diameter_of_local_lattice,
     "number of nodes along diameter of local sub-lattice."},
  {"-lattice_diameter", ARGV_FLOAT, (char *) 3, 
     (char *) main_argsX.lattice_width,
     "widths of sub-lattice along each dimension (X, Y, Z)"},
  {"-max_def_magnitude", ARGV_FLOAT, (char *) 0, 
     (char *) &main_argsX.trans_info.max_def_magnitude,
     "maximum expected deformation magnitude (mm)"},

  {"-use_magnitude", ARGV_CONSTANT, (char *) TRUE, (char *) &main_argsX.trans_info.use_magnitude,
     "use magnitude data local deformation (default)."},

  {"-optical_flow", ARGV_CONSTANT, (char *) FALSE, (char *) &main_argsX.trans_info.use_magnitude,
     "use optical flow to compute deformation."},
  {"-use_simplex", ARGV_CONSTANT, (char *) TRUE, (char *) &main_argsX.trans_info.use_simplex,
     "use 3D simplex optimization for local deformation (default)."},
  {"-quadratic", ARGV_CONSTANT, (char *) FALSE, (char *) &main_argsX.trans_info.use_simplex,
     "use quadratic fit for local deformation."},
  {"-use_local", ARGV_CONSTANT, (char *) TRUE, (char *) &main_argsX.trans_info.use_local_smoothing,
     "Turn on local smoothing (default = global smoothing)."},
  {"-use_nonisotropic", ARGV_CONSTANT, (char *) FALSE, (char *) &main_argsX.trans_info.use_local_isotropic,
     "Turn on directionally sensitive smoothing (def=isotropic smoothing)."},
  {"-super", ARGV_INT, (char *) 0, (char *) &main_argsX.trans_info.use_super,
     "super sample deformation field during optimization (default)."},
  {"-no_super", ARGV_CONSTANT, (char *) 0, (char *) &main_argsX.trans_info.use_super,
     "do not super sample deformation field during optimization."},
  {"-iterations", ARGV_INT, (char *) 0, 
     (char *) &iteration_limit,
     "Number of iterations for non-linear optimization"},
  {"-weight", ARGV_FLOAT, (char *) 0, 
     (char *) &iteration_weight,
     "Weighting factor for each iteration in nl optimization"},
  {"-stiffness", ARGV_FLOAT, (char *) 0, 
     (char *) &smoothing_weight,
     "Weighting factor for smoothing between nl iterations"},
  {"-similarity_cost_ratio", ARGV_FLOAT, (char *) 0, 
     (char *) &similarity_cost_ratio,
     "Weighting factor for  r=similarity*w + cost(1*w)"},

  {NULL, ARGV_HELP, NULL, NULL,
     "\nOptions for logging progress. Default = -verbose 1."},
  {"-verbose", ARGV_INT, (char *) 0, (char *) &main_argsX.flags.verbose,
     "Write messages indicating progress"},
  {"-quiet", ARGV_CONSTANT, (char *) 0 , (char *) &main_argsX.flags.verbose,
     "Do not write log messages"},
  {"-debug", ARGV_CONSTANT, (char *) TRUE, (char *) &main_argsX.flags.debug,
     "Print out debug info."},
  {"-version", ARGV_FUNC, (char *) print_version_info, (char *)MNI_AUTOREG_LONG_VERSION,
     "Print out version info and exit."},
  {NULL, ARGV_END, NULL, NULL, NULL}
};




Arg_Data main_argsX = {
  {"","","","","","",""},        /* filenames           */
  {1,FALSE},                        /* verbose, debug      */
  {                                /* transformation info */
    FALSE,                        /*   use identity tranformation to start */
    TRUE,                        /*   do default tranformation (PAT) to start */
    TRUE,                        /*   use_mag=TRUE; do not use projections by default */
    50.0,
    TRUE,                        /*   use_simplex=TRUE ie use 3d simplex by default */
    FALSE,                       /*   use_bfgs=FALSE i.e don't use BFGS*/
    2,                                /*   use super sampling of deformation field  */
    FALSE,                        /* use local smoothing       */
    TRUE,                        /* use isotropic smoothing */
    "",                                /*   filename */
    NULL,                        /*   file_contents */
    0,                          /* buffer_length   */
    (VIO_General_transform *)NULL,        /*   General transform */
    (VIO_General_transform *)NULL,        /*   General transform copy of input */
    TRANS_LSQ7,                        /*   default type      */
    {-DBL_MAX, -DBL_MAX, -DBL_MAX},                /*   center            */
    {1.0, 1.0, 1.0},                /*   scale             */
    {0.0, 0.0, 0.0},                /*   shears            */
    {0.0, 0.0, 0.0},                /*   translations      */
    {0.0, 0.0, 0.0, 1.0},              /*   quaternions       */
    {0.0, 0.0, 0.0},            /*   rotations         */
    {1.0, 1.0, 1.0,  3.1415927/180.0, 3.1415927/180.0, 3.1415927/180.0, 0.02, 0.02, 0.02,  0.02, 0.02, 0.02}, /* optimization weights*/
    FALSE,                        /*   invert_mapping_flag                  */
    TRANS_ROT},                  /* default use normal rotation */
  {0,NULL, NULL, NULL, NULL, NULL, NULL},        /* FEATURE VOL */
  trilinear_interpolant,        /* use trilinear interpolation by default */
  TRILINEAR,                   /* use trilinear interpolation by default */
  xcorr_objective,              /* use cross-correlation by default       */
  XCORR,                        /* use cross-correlation by default       */
  OPT_SIMPLEX,                  /* use simplex optimization strategy      */
  0,                            /* do not force lattice on source or target */
  {4.0,4.0,4.0},                /* default step sizes for lattice         */
  {24.0,24.0,24.0},                /* default lattice diameter               */
  {0.0,0.0,0.0},                /* default start for lattice, reset in init_lattice */
  {0,0,0},                      /* default number of element in lattice, also reset */

  {{{1.0,0.0,0.0}},                /* default sampling lattice axes directions */
   {{0.0,1.0,0.0}},
   {{0.0,0.0,1.0}}},

  1,                            /* use first volume as default smallest volume      */
  {FALSE, FALSE, FALSE, FALSE},        /* VIO_Transform flags: est_cent, _scale, _rots, _trans */
  {0.0,0.0},                        /* lower limit of voxels considered                 */
  5.0,                                /* percent noise speckle                            */
  256,                                /* number of groups to use for ratio of variance    */
  3                                /* pdf blurring size for -mi                        */
};

Arg_Data *main_args = &main_argsX;

