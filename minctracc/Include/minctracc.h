/* ----------------------------- MNI Header -----------------------------------
@NAME       : minctracc.h
@DESCRIPTION: Header file for minctracc.c
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Thu May 20 14:20:21 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */

#include <minc.h>
#include <ParseArgv.h>
#include "../make_xfm_file/def_tps.h"
#include "../make_xfm_file/def_main.h"


/* ------------------------ Constants used in program  ------------------------ */

#define VOL_NDIMS       3
#define WORLD_NDIMS     3   /* Number of world spatial dimensions */
#define SLICE_NDIMS     2   /* Number of slice dimensions */
#define MAT_NDIMS       WORLD_NDIMS+1   /* Number of dims for homogenous matrix */
#define NO_AXIS         -1
#define DEFAULT_MAX     1.0
#define DEFAULT_MIN     0.0
#define FILL_DEFAULT    DBL_MAX   /* Fillvalue indicating -nofill */
#define SMALL_VALUE     (100.0*FLT_MIN)   /* A small floating-point value */

#ifndef TRUE
#  define TRUE          1
#  define FALSE         0
#endif

#define NORMAL_STATUS   0
#define ERROR_STATUS    1

#define TRANS_PROCRUSTES  0
#define TRANS_LSQ         1
#define TRANS_LSQ6        2
#define TRANS_LSQ7        3
#define TRANS_LSQ9        4
#define TRANS_LSQ12       5
#define TRANS_PAT         6

#define OPT_SIMPLEX       0

/* ------------------------  Types used in program  ------------------------ */

typedef struct {
   double x, y, z;
} Coord_Vector;

typedef struct {
   double mat[WORLD_NDIMS][MAT_NDIMS];
} Linear_Transformation;

typedef void (*Transform_Function)(Coord_Vector *result,
                                   void *trans_data, Coord_Vector *coordinate);

typedef struct {
   int linear;			 /* ==TRUE if linear transformation */
   Transform_Function transform; /* name of the function used to apply the transformation */
   void *trans_data;		 /* pointer to data structure containing the transformation */
} Transformation;

typedef struct {
   char *name;
   int mincid;
   int icvid;
   int imgid;
   int maxid;
   int minid;
   int ndims;
   nc_type datatype;
   int is_signed;
   double vrange[2];                /* [0]=min, [1]=max */
   long nelements[MAX_VAR_DIMS];    /* Size of each dimension */
   int world_axes[MAX_VAR_DIMS];    /* Relates variable index to X, Y, Z 
                                       or NO_AXIS */
   int indices[VOL_NDIMS];        /* Indices of volume dimenions (subscripted
                                       from slowest to fastest) */
   int axes[WORLD_NDIMS];    /* Relates world X,Y,Z (index) to dimension 
                                order (value=0,1,2; 0=slowest varying) */
} File_Info;


typedef struct Volume_Data_Struct Volume_Data;

typedef int (*Interpolating_Function) 
     (Volume volume, Point *coord, double *result);

struct Volume_Data_Struct {
   nc_type datatype;         /* Type of data in volume */
   int is_signed;            /* Sign of data (TRUE if signed) */
   int use_fill;             /* TRUE if fill values should be used in
                                calculation of output image max/min */
   double fillvalue;         /* Value to return when out of bounds */
   int size[VOL_NDIMS];      /* Size of each dimension */
   void *data;               /* Pointer to volume data */
   double *scale;            /* Pointer to array of scales for slices */
   double *offset;           /* Pointer to array of offsets for slices */
   Interpolating_Function interpolant;
};

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
  int use_default;
  Transformation transformation;/* input world to world transformation */
  int transform_type;		/* type of transformation to optimize */
  double center[3];		/* parameters corresponding to trans matrix */
  double scales[3];
  double shears[3];
  double rotations[3];
  double translations[3];
  int invert_mapping_flag;		/* true if input transform maps model to source */
} Program_Transformation;


typedef struct Arg_Data Arg_Data_struct;

typedef float (*Objective_Function) (Volume d1,
				     Volume d2,
				     Volume m1,
				     Volume m2, 
				     Arg_Data_struct *globals);

typedef struct {
  Program_Filenames      filenames;    /* names of all data filename to be used      */
  Program_Flags          flags;	       /* flags (debug, verbose etc...               */ 
  Program_Transformation trans_info;   /* world to world transformation information  */
  Interpolating_Function interpolant;  /* point to interpolation funciton to be used */
  Objective_Function     obj_function; /* pointer to objective function to be used   */
  int                    optimize_type;/* Type of optimization strategy              */
  double                 step[3];      /* step size for sampling lattice             */
  double                 start[3];     /* starting position for sampling lattice     */
  int                    count[3];     /* number of elements for sampling lattice    */
  int                    smallest_vol; /* either one or two, indicates the smaller vol */
  Transform_Flags        trans_flags;  /* flags defining which parameters to estimate*/
  double                 threshold;    /* lower limit of voxels considered           */
  double                 speckle;      /* percent noise speckle                      */
  int                    groups;       /* number of groups to use for ratioo of variance */
} Arg_Data;


/*  ------------------------ Macros used in program  ------------------------ */

#define DO_TRANSFORM(result, transformation, coord) \
   (*transformation->transform) (result, transformation->trans_data, coord)

#define IS_LINEAR(transformation) (transformation->linear)

#define VECTOR_DIFF(result, first, second) { \
   result.x = first.x - second.x; \
   result.y = first.y - second.y; \
   result.z = first.z - second.z; \
}

#define VECTOR_ADD(result, first, second) { \
   result.x = first.x + second.x; \
   result.y = first.y + second.y; \
   result.z = first.z + second.z; \
}

#define VECTOR_SCALAR_MULT(result, vector, scalar) { \
   result.x = vector.x * scalar; \
   result.y = vector.y * scalar; \
   result.z = vector.z * scalar; \
}

#define VOLUME_VALUE(volume, ind0, ind1, ind2, value) \
{ \
   long _offset_; \
 \
   _offset_ = ((ind0)*volume->size[1] + (ind1))*volume->size[2] + (ind2); \
   switch (volume->datatype) { \
   case NC_BYTE: \
      if (volume->is_signed) \
         value = *((signed char *) volume->data + _offset_); \
      else \
         value = *((unsigned char *) volume->data + _offset_); \
      break; \
   case NC_SHORT: \
      if (volume->is_signed) \
         value = *((signed short *) volume->data + _offset_); \
      else \
         value = *((unsigned short *) volume->data + _offset_); \
      break; \
   case NC_LONG: \
      if (volume->is_signed) \
         value = *((signed long *) volume->data + _offset_); \
      else \
         value = *((unsigned long *) volume->data + _offset_); \
      break; \
   case NC_FLOAT: \
      value = *((float *) volume->data + _offset_); \
      break; \
   case NC_DOUBLE: \
      value = *((double *) volume->data + _offset_); \
      break; \
   } \
}

#ifndef DEBUG_PRINT
#   define DEBUG_PRINT(str) if (main_args.flags.debug) (void) fprintf (stderr,  str  );
#   define DEBUG_PRINT1(str,a1) if (main_args.flags.debug) (void) fprintf (stderr,  str ,a1 );
#   define DEBUG_PRINT2(str,a1,a2) if (main_args.flags.debug) (void) fprintf (stderr,  str ,a1,a2 );
#   define DEBUG_PRINT3(str,a1,a2,a3) if (main_args.flags.debug) (void) fprintf (stderr,  str ,a1,a2,a3 );
#   define DEBUG_PRINT4(str,a1,a2,a3,a4) if (main_args.flags.debug) (void) fprintf (stderr,  str ,a1,a2,a3,a4 );
#   define DEBUG_PRINT5(str,a1,a2,a3,a4,a5) if (main_args.flags.debug) (void) fprintf (stderr,  str ,a1,a2,a3,a4,a5 );
#   define DEBUG_PRINT6(str,a1,a2,a3,a4,a5,a6) if (main_args.flags.debug) (void) fprintf (stderr,  str ,a1,a2,a3,a4,a5,a6 );
#   define DEBUG_PRINT7(str,a1,a2,a3,a4,a5,a6,a7) if (main_args.flags.debug) (void) fprintf (stderr,  str ,a1,a2,a3,a4,a5,a6,a7 );
#endif


/*  ------------------------ Function prototypes  ------------------------ */

public void resample_volumes(Program_Flags *program_flags,
                             Volume in_vol, Volume out_vol, 
                             Transformation *transformation);

public void invert_transformation(Transformation *result, 
                                  Transformation *transformation);

public void do_linear_transformation(Coord_Vector *result, void *trans_data, 
                                     Coord_Vector *coordinate);

public void do_linear_transformation_point(Point *result, void *trans_data, 
                                     Point *coordinate);

public void do_non_linear_transformation(Coord_Vector *result, void *trans_data, 
                                     Coord_Vector *coordinate);

public int trilinear_interpolant(Volume volume, 
                                 Point *coord, double *result);

public int tricubic_interpolant(Volume volume, 
                                Point *coord, double *result);

public void do_Ncubic_interpolation(Volume volume, 
                                    long index[], int cur_dim, 
                                    double frac[], double *result);

public int nearest_neighbour_interpolant(Volume volume, 
                                         Point *coord, double *result);

public void mult_linear_transform(Transformation *result, 
                                  Transformation *transform1, 
                                  Transformation *transform2);

public int get_transformation(char *dst, char *key, char *nextArg);

public int get_mask_file(char *dst, char *key, char *nextArg);

public void procrustes(int npoints, int ndim, 
                       float **Apoints, float **Bpoints,
                       float *translation, float *centre_of_rotation,
                       float **rotation, float *scale);

public void transformations_to_homogeneous(int ndim, 
                  float *translation, float *centre_of_rotation,
                  float **rotation, float scale,
                  float **transformation);

public void translation_to_homogeneous(int ndim, float *translation,
                                       float **transformation);

public void rotation_to_homogeneous(int ndim, float **rotation,
                                       float **transformation);

public float fit_function(float *x);        /* apply cross correlation to the data sets    */

public float zscore_function(float *x);     /* calculate rms z-score difference.           */

public float check_function(float *x);      /* calculate the squared error between points2 */

public void invertmatrix(int n, float **mat, float **mat_invert);

public int   save_transform(
    global_data_struct   *main,
    char          filename[] );


public Boolean init_params(Volume d1,
			   Volume d2,
			   Volume m1,
			   Volume m2, 
			   Arg_Data *globals);

public void init_lattice(Volume d1,
			 Volume d2,
			 Volume m1,
			 Volume m2, 
			 Arg_Data *globals);

public Boolean optimize_linear_transformation(Volume d1,
					      Volume d2,
					      Volume m1,
					      Volume m2, 
					      Arg_Data *globals);

public float xcorr_objective(Volume d1,
			     Volume d2,
			     Volume m1,
			     Volume m2, 
			     Arg_Data *globals);

public float zscore_objective(Volume d1,
			      Volume d2,
			      Volume m1,
			      Volume m2, 
			      Arg_Data *globals);

public float vr_objective(Volume d1,
			  Volume d2,
			  Volume m1,
			  Volume m2, 
			  Arg_Data *globals);

public float ssc_objective(Volume d1,
			   Volume d2,
			   Volume m1,
			   Volume m2, 
			   Arg_Data *globals);

public void print_error(char *s, char * d1, int d2, int d3, int d4, int d5, int d6, int d7);

public float measure_fit(Volume d1,
			 Volume d2,
			 Volume m1,
			 Volume m2, 
			 Arg_Data *globals);

Arg_Data main_args = {
  {NULL,NULL,NULL,NULL,NULL,NULL,NULL},	/* filenames           */
  {1,FALSE},			/* verbose, debug      */
  {				/* transformation info */
    TRUE,			/*   do default Principal Axis Transformation start */
    {TRUE, NULL, NULL},		/*   transformation= linear, rest set in main, after ParseArg */
    TRANS_PROCRUSTES,		/*   default type      */
    {0.0, 0.0, 0.0},		/*   center            */
    {1.0, 1.0, 1.0},		/*   scale             */
    {0.0, 0.0, 0.0},		/*   shears            */
    {0.0, 0.0, 0.0},		/*   rotations         */
    {0.0, 0.0, 0.0},		/*   translations      */
    FALSE},			/*   invert_mapping_flag                  */
  trilinear_interpolant,	/* use trilinear interpolation by default */
  xcorr_objective,              /* use cross-correlation by default       */
  OPT_SIMPLEX,                  /* use simplex optimization strategy      */
  {4.0,4.0,4.0},		/* default step sizes for lattice         */
  {0.0,0.0,0.0},		/* default start for lattice, reset in init_lattice */
  {0,0,0},                      /* default number of element in lattice, also reset */
  1,                            /* use first volume as default smallest volume      */
  {FALSE, FALSE, FALSE, FALSE},
  0.0,				/* lower limit of voxels considered                 */
  5.0,				/* percent noise speckle                            */
  16				/* number of groups to use for ratio of variance    */
};

#define INTERPOLATE_VOXEL_VALUE(volume, coord, result) \
   (*(main_args.interpolant)) (volume, coord, result)

