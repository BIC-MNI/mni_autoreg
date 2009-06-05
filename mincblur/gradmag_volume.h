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
@MODIFIED   : Revision 96.3  2009-06-05 20:49:52  claude
@MODIFIED   : Free memory after usage in mincblur
@MODIFIED   :
@MODIFIED   : Revision 96.2  2006/11/28 09:12:21  rotor
@MODIFIED   :  * fixes to allow clean compile against minc 2.0
@MODIFIED   :
@MODIFIED   : Revision 96.1  2004/02/12 05:53:48  rotor
@MODIFIED   :  * removed public/private defs
@MODIFIED   :
@MODIFIED   : Revision 96.0  1996/08/21 18:22:24  louis
@MODIFIED   : Release of MNI_AutoReg version 0.96
@MODIFIED   :
 * Revision 9.6  1996/08/21  18:22:19  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:16:28  louis
 * Never released version 0.95
 *
 * Revision 1.8  1996/08/12  14:16:22  louis
 * Pre-release
 *
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
  double vrange[2];                /* [0]=min, [1]=max */
  long nelements[MAX_VAR_DIMS];        /* Size of each dimension */
  int  world_axes[MAX_VAR_DIMS];        /* Relates variable index to X, Y, Z 
                                   or NO_AXIS */ 
  int  indices[VOL_NDIMS];        /* Indices of volume dimenions (subscripted
                                   from slowest to fastest) */
  int  axes[WORLD_NDIMS];        /* Relates world X,Y,Z (index) to dimension 
                                   order (value=0,1,2; 0=slowest varying) */
} File_Info;

typedef struct {
  nc_type datatype;                /* Type of data in volume */
  int is_signed;                /* Sign of data (TRUE if signed) */
  int use_fill;                        /* TRUE if fill values should be used in
                                   calculation of output image max/min */
  double fillvalue;                /* Value to return when out of bounds */
  int size[VOL_NDIMS];                /* Size of each dimension */
  void *data;                        /* Pointer to volume data */
  double *scale;                /* Pointer to array of scales for slices */
  double *offset;                /* Pointer to array of offsets for slices */
} Volume_Data;

typedef struct {
   long size[SLICE_NDIMS];        /* Size of each dimension */
   double *data;                /* Pointer to slice data */
} Slice_Data;

typedef struct {
  File_Info   *file;                /* Information about associated file */
  Volume_Data *volume;                /* VIO_Volume data for (input volume) */
  Slice_Data  *slice;                /* Slice data for (output volume) */
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

void build_vol_info(char *infilename, char *outfilename,
                        MincVolume *in_vol, MincVolume *out_vol, char *history);

void load_vol_info(char *infilename, MincVolume *in_vol);
void free_vol_info( MincVolume * vol );

void get_file_info(char *filename, 
                          Volume_Definition *volume_def,
                          File_Info *file_info);

void create_output_file(char *filename, 
                               Volume_Definition *volume_def,
                               File_Info *in_file,
                               File_Info *out_file,
                               char *tm_stamp);

double get_default_range(char *what, nc_type datatype, int is_signed);

void finish_up(MincVolume *in_vol1, 
                      MincVolume *in_vol2, 
                      MincVolume *in_vol3, 
                      MincVolume *out_vol);

void load_volume(File_Info *file, long start[], long count[], 
                        Volume_Data *volume);

void get_mag_slice(Slice_Data *result,
                          Slice_Data *slice_dx,
                          Slice_Data *slice_dy,
                          Slice_Data *slice_dz,                          
                          double *minimum, double *maximum);


void get_curvature_slice(Slice_Data *result,
                                Slice_Data *slice_dx,
                                Slice_Data *slice_dy,
                                Slice_Data *slice_dz,
                                Slice_Data *slice_dxx,
                                Slice_Data *slice_dyy,
                                Slice_Data *slice_dzz,
                                double thresh,
                                double *minimum, double *maximum);

void make_gradmag_volumes(MincVolume *in_vol1, 
                                 MincVolume *in_vol2, 
                                 MincVolume *in_vol3, 
                                 MincVolume *out_vol,
                                 double *min_value, double *max_value);

void make_curvature_volumes(MincVolume *in_vol1, 
                                   MincVolume *in_vol2, 
                                   MincVolume *in_vol3, 
                                   MincVolume *in_volxx, 
                                   MincVolume *in_volyy, 
                                   MincVolume *in_volzz, 
                                   MincVolume *out_vol,
                                   double thresh);

void make_vol_icv(MincVolume *in_vol);

void calc_gradient_magnitude(char *infilename, char *history, char *output_basename,
                                    double *min_value, double *max_value);

void calc_gaussian_curvature(char *infilename, char *history,
                                    double min_value, double max_value);
