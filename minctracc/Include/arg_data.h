#include "point_vector.h"

typedef struct Arg_Data_struct Arg_Data;

typedef float (*Objective_Function) (Volume d1,
				     Volume d2,
				     Volume m1,
				     Volume m2, 
				     Arg_Data *globals);

typedef void (*Transform_Function)(PointR *result,
                                   General_transform *trans_data, PointR *coordinate);

typedef int (*Interpolating_Function) 
     (Volume volume, PointR *coord, double *result);

typedef struct {
   int verbose;
   int debug;
} Program_Flags;

typedef struct {
   int estimate_center;
   int estimate_scale;
   int estimate_rots;
   int estimate_trans;
} Transform_Flags;

typedef struct {
  char *data;
  char *model;
  char *mask_data;
  char *mask_model;
  char *output_trans;
  char *measure_file;
  char *matlab_file;
} Program_Filenames;

typedef struct {
  int number_of_features;
  Volume *data;
  Volume *model;
  Volume *data_mask;
  Volume *model_mask;
  char **data_name;
  char **model_name;
  char **mask_data_name;
  char **mask_model_name;
  char *obj_func;
  Real *weight;
  Real *thresh_data;
  Real *thresh_model;
} Feature_volumes;

typedef struct {
  int use_identity;
  int use_default;
  int use_magnitude;
  int use_simplex;
  int use_super;
  int use_local_smoothing;
  int use_local_isotropic;
  char *file_name;
  char *file_contents;
  long buffer_length;
  General_transform *transformation;     /* optimized world to world transformation */
  General_transform *orig_transformation;/* input world to world transformation */
  int transform_type;		/* type of transformation to optimize */
  double center[3];		/* parameters corresponding to trans matrix */
  double scales[3];
  double shears[3];
  double rotations[3];
  double translations[3];
  double weights[12];		/* optimization weighting function */
  int invert_mapping_flag;	/* true if input transform maps model to source */
} Program_Transformation;

struct Arg_Data_struct {
  Program_Filenames      filenames;    /* names of all data filename to be used      */
  Program_Flags          flags;	       /* flags (debug, verbose etc...               */ 
  Program_Transformation trans_info;   /* world to world transformation information  */
  Feature_volumes        features;     /* struct contain extra feature info */

  Interpolating_Function interpolant;  /* point to interpolation funciton to be used */
  Objective_Function     obj_function; /* pointer to objective function to be used   */
  int                    optimize_type;/* Type of optimization strategy              */
  int                    force_lattice;/* =0, do not force; =1, force src; =2, frc target */
  double                 step[3];      /* step size for sampling lattice             */
  double                 start[3];     /* starting position for sampling lattice     */
  int                    count[3];     /* number of elements for sampling lattice    */
  VectorR                directions[3];/* directions for each axis of sampling lat   */
  int                    smallest_vol; /* either one or two, indicates the smaller vol */
  Transform_Flags        trans_flags;  /* flags defining which parameters to estimate*/
  double                 threshold[2]; /* lower limit of voxels considered           */
  double                 speckle;      /* percent noise speckle                      */
  int                    groups;       /* number of groups to use for ratio of variance */
  int                    blur_pdf;     /* number of voxels for blurring in -mi pdfs */
};

