/* ----------------------------- MNI Header -----------------------------------
@NAME       : gradmag_volume.h
@DESCRIPTION: Header/prototype file for gradmag_volume.c
@COPYRIGHT  :
              Copyright 1995 Louis Collins, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.

@CREATED    : Wed Jun 30 13:25:33 EST 1993 Louis Collins
                copied and modified from mincresample.h from Peter Neelin
@MODIFIED   : $Log: gradmag_volume.h,v $
@MODIFIED   : Revision 1.8  1996-08-12 14:16:22  louis
@MODIFIED   : Pre-release
@MODIFIED   :
---------------------------------------------------------------------------- */

/* Constants used in program */

#define VOL_NDIMS    3   /* Number of volume dimensions */
#define WORLD_NDIMS  3   /* Number of world spatial dimensions */
#define SLICE_NDIMS  2   /* Number of slice dimensions */
#define NO_AXIS -1
#define X 0
#define Y 1
#define Z 2
#define DEFAULT_MAX     1.0
#define DEFAULT_MIN     0.0
#define FILL_DEFAULT    DBL_MAX   /* Fillvalue indicating -nofill */
#define SMALL_VALUE     (100.0*FLT_MIN)   /* A small floating-point value */
#define PROCESSING_VAR "processing"
#define TEMP_IMAGE_VAR "gradmag_volume-temporary-image"
#ifndef TRUE
#  define TRUE 1
#  define FALSE 0
#endif

/* Types used in program */

typedef struct {
  char *name;
  int  mincid;
  int  icvid;
  int  imgid;
  int  maxid;
  int  minid;
  int  ndims;
  nc_type datatype;
  int  is_signed;
  double vrange[2];		/* [0]=min, [1]=max */
  long nelements[MAX_VAR_DIMS];	/* Size of each dimension */
  int  world_axes[MAX_VAR_DIMS];	/* Relates variable index to X, Y, Z 
				   or NO_AXIS */ 
  int  indices[VOL_NDIMS];	/* Indices of volume dimenions (subscripted
				   from slowest to fastest) */
  int  axes[WORLD_NDIMS];	/* Relates world X,Y,Z (index) to dimension 
				   order (value=0,1,2; 0=slowest varying) */
} File_Info;

typedef struct {
  nc_type datatype;		/* Type of data in volume */
  int is_signed;		/* Sign of data (TRUE if signed) */
  int use_fill;			/* TRUE if fill values should be used in
				   calculation of output image max/min */
  double fillvalue;		/* Value to return when out of bounds */
  int size[VOL_NDIMS];		/* Size of each dimension */
  void *data;			/* Pointer to volume data */
  double *scale;		/* Pointer to array of scales for slices */
  double *offset;		/* Pointer to array of offsets for slices */
} Volume_Data;

typedef struct {
   long size[SLICE_NDIMS];	/* Size of each dimension */
   double *data;		/* Pointer to slice data */
} Slice_Data;

typedef struct {
  File_Info   *file;		/* Information about associated file */
  Volume_Data *volume;		/* Volume data for (input volume) */
  Slice_Data  *slice;		/* Slice data for (output volume) */
} MincVolume;

typedef struct {
  int axes[WORLD_NDIMS];    /* Relates world X,Y,Z (index) to dimension 
			       order (value=0,1,2; 0=slowest varying) */
  long nelements[WORLD_NDIMS]; /* These are subscripted by X, Y and Z */
  double step[WORLD_NDIMS];
  double start[WORLD_NDIMS];
  double dircos[WORLD_NDIMS][WORLD_NDIMS];
  char units[WORLD_NDIMS][MI_MAX_ATTSTR_LEN];
  char spacetype[WORLD_NDIMS][MI_MAX_ATTSTR_LEN];
} Volume_Definition;

typedef struct {
  nc_type datatype;
  int is_signed;
  double vrange[2];
  double fillvalue;
  Volume_Definition volume_def;
} Default_Data;

/* Macros used in program */

#define MALLOC(size) malloc(size)

#define FREE(ptr) free(ptr)

/* Function prototypes */

public void build_vol_info(char *infilename, char *outfilename,
                        MincVolume *in_vol, MincVolume *out_vol, char *history);

public void load_vol_info(char *infilename, MincVolume *in_vol);

public void get_file_info(char *filename, 
                          Volume_Definition *volume_def,
                          File_Info *file_info);

public void create_output_file(char *filename, 
                               Volume_Definition *volume_def,
                               File_Info *in_file,
                               File_Info *out_file,
                               char *tm_stamp);

public double get_default_range(char *what, nc_type datatype, int is_signed);

public void finish_up(MincVolume *in_vol1, 
		      MincVolume *in_vol2, 
		      MincVolume *in_vol3, 
		      MincVolume *out_vol);

public void load_volume(File_Info *file, long start[], long count[], 
                        Volume_Data *volume);

public void get_mag_slice(Slice_Data *result,
			  Slice_Data *slice_dx,
			  Slice_Data *slice_dy,
			  Slice_Data *slice_dz,			  
			  double *minimum, double *maximum);


public void get_curvature_slice(Slice_Data *result,
				Slice_Data *slice_dx,
				Slice_Data *slice_dy,
				Slice_Data *slice_dz,
				Slice_Data *slice_dxx,
				Slice_Data *slice_dyy,
				Slice_Data *slice_dzz,
				double thresh,
				double *minimum, double *maximum);

public void make_gradmag_volumes(MincVolume *in_vol1, 
				 MincVolume *in_vol2, 
				 MincVolume *in_vol3, 
				 MincVolume *out_vol,
				 double *min_value, double *max_value);

public void make_curvature_volumes(MincVolume *in_vol1, 
				   MincVolume *in_vol2, 
				   MincVolume *in_vol3, 
				   MincVolume *in_volxx, 
				   MincVolume *in_volyy, 
				   MincVolume *in_volzz, 
				   MincVolume *out_vol,
				   double thresh);

public void make_vol_icv(MincVolume *in_vol);

public void calc_gradient_magnitude(char *infilename, char *history, char *output_basename,
				    double *min_value, double *max_value);

public void calc_gaussian_curvature(char *infilename, char *history,
				    double min_value, double max_value);
