#include "point_vector.h"

typedef struct Arg_Data_struct Arg_Data;

/* enums to define interpolants and objective functions */
enum Interpolating_Type { TRILINEAR, TRICUBIC, N_NEIGHBOUR };
enum Objective_Type { XCORR, ZSCORE, SSC, VR, MUTUAL_INFORMATION, NORMALIZED_MUTUAL_INFORMATION };


typedef float (*Objective_Function) (VIO_Volume d1,
                                     VIO_Volume d2,
                                     VIO_Volume m1,
                                     VIO_Volume m2, 
                                     Arg_Data *globals);

typedef void (*Transform_Function)(PointR *result,
                                   VIO_General_transform *trans_data, PointR *coordinate); /* transforme un point dans ces nouvelles coordonnees */

typedef int (*Interpolating_Function) 
     (VIO_Volume volume, PointR *coord, double *result);

typedef struct {
   int verbose;
   int debug;
} Program_Flags;

typedef struct {
   int estimate_center;
   int estimate_scale;
   int estimate_trans;
   int estimate_rots;
   int estimate_quats;
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
  VIO_Volume *data;
  VIO_Volume *model;
  VIO_Volume *data_mask;
  VIO_Volume *model_mask;
  char **data_name;
  char **model_name;
  char **mask_data_name;
  char **mask_model_name;
  char *obj_func;
  VIO_Real *weight;
  VIO_Real *thresh_data;
  VIO_Real *thresh_model;
} Feature_volumes;

typedef struct {
  int use_identity;
  int use_default;
  int use_magnitude;
  VIO_Real max_def_magnitude;        /* maximum size of deformation in def field */
  int use_simplex;
  int use_super;
  int use_local_smoothing;
  int use_local_isotropic;
  char *file_name;
  char *file_contents;
  long buffer_length;
  VIO_General_transform *transformation;     /* optimized world to world transformation */
  VIO_General_transform *orig_transformation;/* input world to world transformation */
  int transform_type;                /* type of transformation to optimize */
  double center[3];                /* parameters corresponding to trans matrix */
  double scales[3];
  double shears[3];
  double translations[3];
  double quaternions[4];  
  double rotations[3];
  double weights[12];        /* optimization weighting function with quaternions */
  int invert_mapping_flag;        /* true if input transform maps model to source */
  int rotation_type;            /* type of rotation quaternion used or not */
} Program_Transformation;

struct Arg_Data_struct {
  Program_Filenames      filenames;    /* names of all data filename to be used      */
  Program_Flags          flags;               /* flags (debug, verbose etc...               */ 
  Program_Transformation trans_info;   /* world to world transformation information  */
  Feature_volumes        features;     /* struct contain extra feature info */

  Interpolating_Function interpolant;  /* point to interpolation funciton to be used */
  enum Interpolating_Type interpolant_type;  /* enum defining interpolant to be used */
  Objective_Function     obj_function; /* pointer to objective function to be used   */
  enum Objective_Type    obj_function_type; /* enum defining objective function to be used */
  int                    optimize_type;/* Type of optimization strategy              */
  int                    force_lattice;/* =0, do not force; =1, force src; =2, frc target */
  double                 step[3];      /* step size for sampling lattice             */
  double                 lattice_width[3];   /* diameter (mm) of sub-lattice         */
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

