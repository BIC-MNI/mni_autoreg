/* ----------------------------- MNI Header -----------------------------------
@NAME       : gradmag_volume.c
@INPUT      : infilename - base file name of data
              history    - a comment referring to volume's history
@OUTPUT     : a gradient magnitude minc volume, saved to disk.
@RETURNS    : error status
@DESCRIPTION: calculates the gradient magnitude of three
              parital derivative volumes.
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

@CREATED    : some time in 1993
@MODIFIED   : $Log: gradmag_volume.c,v $
@MODIFIED   : Revision 96.5  2009-06-05 20:49:52  claude
@MODIFIED   : Free memory after usage in mincblur
@MODIFIED   :
@MODIFIED   : Revision 96.4  2006/11/28 09:12:21  rotor
@MODIFIED   :  * fixes to allow clean compile against minc 2.0
@MODIFIED   :
@MODIFIED   : Revision 96.3  2005/07/20 20:45:39  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 96.2  2005/01/29 16:34:05  rotor
@MODIFIED   :  * removed malloc.h include for smoother OSX build
@MODIFIED   :
@MODIFIED   : Revision 96.1  2004/02/12 05:53:48  rotor
@MODIFIED   :  * removed /static defs
@MODIFIED   :
@MODIFIED   : Revision 96.0  1996/08/21 18:22:24  louis
@MODIFIED   : Release of MNI_AutoReg version 0.96
@MODIFIED   :
 * Revision 9.6  1996/08/21  18:22:18  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:16:28  louis
 * Never released version 0.95
 *
 * Revision 1.12  1996/08/12  14:16:22  louis
 * Pre-release
 *
 * Revision 1.11  1995/09/18  09:02:42  collins
 * new functional version of mincblur with new and improved default behavior.
 * By default, only the blurred volume is created. If you want the gradient
 * magnitude volumes, you have to ask for it (-gradient).
 *
 * The temporary files (corresponding to the REAL data, and partial derivatives)
 * are removed by mincblur, so now it runs cleanly.  Unfortunately, these files
 * are _not_ deleted if mincblur fails, or is stopped.  I will have to use
 * unlink for this, but its a bit to much to do right now, since I would have
 * to change the way files are dealt with in gradient_volume.c, blur_volume.c
 * and gradmag_volume.c.
 *
 * Also, it is possible to keep the partial derivative volumes (-partial).
 *
 * this version is in mni_reg-0.1j
 *
 * Revision 1.10  1995/09/18  06:45:42  collins
 * this file is a working version of mincblur.  All references to numerical
 * recipes routines have been removed.  This version is included in the
 * package mni_reg-0.1i
 *
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/mincblur/gradmag_volume.c,v 96.5 2009-06-05 20:49:52 claude Exp $";
#endif

#include <config.h>             /* for EXIT_FAILURE (on some systems) */
#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <minc.h>
#include <ParseArgv.h>
#include <time_stamp.h>
#include "gradmag_volume.h"

/* Main program */

extern int verbose;
extern int debug;

void calc_gradient_magnitude(char *infilename, 
                                    char *output_basename,
                                    char *history, 
                                    double *min_value, double *max_value)
{   
  MincVolume 
    in_vol_struct1, 
    in_vol_struct2, 
    in_vol_struct3, 
    out_vol_struct;
  MincVolume 
    *in_vol1 = &in_vol_struct1, 
    *in_vol2 = &in_vol_struct2, 
    *in_vol3 = &in_vol_struct3, 
    *out_vol = &out_vol_struct;
  char
    infilename1[1024],
    infilename2[1024],
    infilename3[1024],
    fulloutfilename[1024];

  sprintf (infilename1,"%s_dx.mnc",infilename);
  sprintf (infilename2,"%s_dy.mnc",infilename);
  sprintf (infilename3,"%s_dz.mnc",infilename);
  sprintf (fulloutfilename,"%s_dxyz.mnc",output_basename);

  build_vol_info(infilename1, fulloutfilename,  in_vol1, out_vol, history);

  load_vol_info(infilename2,in_vol2);
  load_vol_info(infilename3,in_vol3);

  make_vol_icv(in_vol1);
  make_vol_icv(in_vol2);
  make_vol_icv(in_vol3);
  
  /* calculate the gradient magnitude */
  make_gradmag_volumes(in_vol1, in_vol2, in_vol3, out_vol, min_value, max_value);

  /* Finish up */
  finish_up(in_vol1, in_vol2, in_vol3, out_vol);

  /* Release memory */
  free_vol_info(out_vol);
  free_vol_info(in_vol1);
  free_vol_info(in_vol2);
  free_vol_info(in_vol3);
  
}

void calc_gaussian_curvature(char *infilename, char *history,
                                    double min_value, double max_value)
{
   
  MincVolume 
    in_vol_structx, 
    in_vol_structy, 
    in_vol_structz, 
    in_vol_structxx, 
    in_vol_structyy, 
    in_vol_structzz, 
    out_vol_struct;
  MincVolume 
    *in_volx  = &in_vol_structx, 
    *in_voly  = &in_vol_structy, 
    *in_volz  = &in_vol_structz, 
    *in_volxx = &in_vol_structxx, 
    *in_volyy = &in_vol_structyy, 
    *in_volzz = &in_vol_structzz, 
    *out_vol = &out_vol_struct;
  char
    infilenamex[1024],
    infilenamey[1024],
    infilenamez[1024],
    infilenamexx[1024],
    infilenameyy[1024],
    infilenamezz[1024],
    fulloutfilename[1024];
  double thresh;


  thresh = 0.10*(max_value - min_value) + min_value;

if (debug)
  (void)printf("min = %f, max = %f, thresh = %f\n",min_value, max_value, thresh);

  sprintf (infilenamex,"%s_dx.mnc",infilename);
  sprintf (infilenamey,"%s_dy.mnc",infilename);
  sprintf (infilenamez,"%s_dz.mnc",infilename);
  sprintf (infilenamexx,"%s_dxx.mnc",infilename);
  sprintf (infilenameyy,"%s_dyy.mnc",infilename);
  sprintf (infilenamezz,"%s_dzz.mnc",infilename);
  sprintf (fulloutfilename,"%s_gcur.mnc",infilename);

  build_vol_info(infilenamex, fulloutfilename,  in_volx, out_vol, history);

  load_vol_info(infilenamey, in_voly);
  load_vol_info(infilenamez, in_volz);
  load_vol_info(infilenamexx,in_volxx);
  load_vol_info(infilenameyy,in_volyy);
  load_vol_info(infilenamezz,in_volzz);

  make_vol_icv(in_volx);
  make_vol_icv(in_voly);
  make_vol_icv(in_volz);
  make_vol_icv(in_volxx);
  make_vol_icv(in_volyy);
  make_vol_icv(in_volzz);
  
  /* calculate the gradient magnitude */
  make_curvature_volumes(in_volx, in_voly, in_volz, 
                         in_volxx, in_volyy, in_volzz, 
                         out_vol,
                         thresh);
  
  /* Finish up */
  finish_up(in_volx, in_voly, in_volz, out_vol);
  finish_up(in_volxx, in_volyy, in_volzz, (MincVolume *)NULL);

  /* Release memory */
  free_vol_info(out_vol);
  free_vol_info(in_volx);
  free_vol_info(in_voly);
  free_vol_info(in_volz);
  free_vol_info(in_volxx);
  free_vol_info(in_volyy);
  free_vol_info(in_volzz);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_vol_icv
@INPUT      : in_vol - description of input volume.
@OUTPUT     : 
@RETURNS    : (nothing)
@DESCRIPTION: Build the image conversion variable and attach it
              to the image variable of the minc file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Wed Jun 30 09:01:51 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */
void make_vol_icv(MincVolume *in_vol)
{
  File_Info *fp;
  
  fp = in_vol->file;

  fp->icvid = miicv_create();

  (void) miicv_setint(fp->icvid, MI_ICV_TYPE,      NC_DOUBLE);
  (void) miicv_setstr(fp->icvid, MI_ICV_SIGN,      MI_SIGNED);
  (void) miicv_setdbl(fp->icvid, MI_ICV_VALID_MAX, DBL_MAX);
  (void) miicv_setdbl(fp->icvid, MI_ICV_VALID_MIN, DBL_MIN);
  (void) miicv_setint(fp->icvid, MI_ICV_USER_NORM, TRUE);
  (void) miicv_setint(fp->icvid, MI_ICV_DO_NORM,   TRUE);

  /* Attach image variable */

  (void) miicv_attach(fp->icvid, fp->mincid, fp->imgid);
  
}  



/* ----------------------------- MNI Header -----------------------------------
@NAME       : build_vol_info
@INPUT      : names of input and output filenames
@OUTPUT     : in_vol - description of input volume.
              out_vol - description of output volume.
@RETURNS    : (nothing)
@DESCRIPTION: Routine to get information input and 
              output files. Sets up all structures
              completely (including allocating space for data).
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Mon Jun 28 14:16:07 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */
void build_vol_info(char   *infile,
                           char   *outfile,
                           MincVolume *in_vol, 
                           MincVolume *out_vol, char *history)
{
   /* Argument parsing information */
   static Default_Data defaults={
      NC_SHORT,               /* Type will be modified anyway */
      INT_MIN,                /* Flag that is_signed has not been set */
      {-DBL_MAX, -DBL_MAX},   /* Flag that range not set */
      FILL_DEFAULT,           /* Flag indicating that fillvalue not set */
      {
          {2, 1, 0},          /* Axis order */
          {0, 0, 0},          /* nelements will be modified */
          {1.0, 1.0, 1.0},    /* Default step */
          {0.0, 0.0, 0.0},    /* Default start */
          {{1.0, 0.0, 0.0},   /* Default direction cosines */
           {0.0, 1.0, 0.0},
           {0.0, 0.0, 1.0}},
          {"", "", ""},       /* units */
          {"", "", ""}        /* spacetype */
       }
   };


   /* Other variables */
   int idim, index, ivar, varid;
   int ndims, dim[MAX_VAR_DIMS], imgdim[MAX_VAR_DIMS];
   int in_vindex, out_vindex;  /* VIO_Volume indices (0, 1 or 2) */
   int in_findex, out_findex;  /* File indices (0 to ndims-1) */
   long size, total_size, total_slice_size;
   File_Info *fp;
   

   /* Check input file for default argument information */
   in_vol->file = MALLOC(sizeof(File_Info));
   get_file_info(infile, &defaults.volume_def, in_vol->file);

   /* Get for input volume data information */
   in_vol->slice = NULL;
   in_vol->volume = MALLOC(sizeof(Volume_Data));
   in_vol->volume->datatype  = in_vol->file->datatype;
   in_vol->volume->is_signed = in_vol->file->is_signed;
   in_vol->volume->fillvalue = 0.0;
   in_vol->volume->use_fill = TRUE;
   
   /* Get space for volume data */
   total_slice_size = total_size = 1;
   for (idim=0; idim < WORLD_NDIMS; idim++) {
      index = defaults.volume_def.axes[idim];
      size = defaults.volume_def.nelements[idim];
      total_size *= size;
      in_vol->volume->size[index] = size;
      if (index != 0) 
        total_slice_size *= size;
   }


/* MALLOC OF DATA SPACE */

   in_vol->slice = MALLOC(sizeof(Slice_Data));
   in_vol->slice->data = MALLOC((size_t) total_slice_size * sizeof(double));

   /* Get space for slice scale and offset */
   in_vol->volume->scale  = MALLOC(sizeof(double) * in_vol->volume->size[0]);
   in_vol->volume->offset = MALLOC(sizeof(double) * in_vol->volume->size[0]);

   /* Check min/max variables */
   fp = in_vol->file;
   if ((fp->maxid != MI_ERROR) && (fp->minid != MI_ERROR)) {
     
     /* Get MIimage dimensions */
     (void) ncvarinq(fp->mincid, fp->imgid, NULL, NULL, &ndims, imgdim, NULL);
     
     /* Check MIimagemax/min dimensions */
     for (ivar=0; ivar<2; ivar++) {
       varid = (ivar==0 ? fp->maxid : fp->minid);
       (void) ncvarinq(fp->mincid, varid, NULL, NULL, &ndims, dim, NULL);
       for (idim=0; idim < ndims; idim++) {
         if ((dim[idim] == imgdim[fp->indices[1]]) ||
             (dim[idim] == imgdim[fp->indices[2]])) {
           (void) fprintf(stderr, 
                          "MIimagemax/min vary over slice dimensions.\n");
           exit(EXIT_FAILURE);
         }
       }        /* End loop over MIimagemax/min dimensions */
     }        /* End loop over variables MIimagemax/min */
     
   }        /* If both MIimagemax/min exist */
   
   /* Set the default output file datatype */
   defaults.datatype = in_vol->file->datatype;




   /* Set up the file description for the output file */

   out_vol->file = MALLOC(sizeof(File_Info));
   out_vol->file->ndims     = in_vol->file->ndims;
   out_vol->file->datatype  = in_vol->file->datatype;
   out_vol->file->is_signed = in_vol->file->is_signed;
   out_vol->file->vrange[0] = in_vol->file->vrange[0];
   out_vol->file->vrange[1] = in_vol->file->vrange[1];
   for (idim=0; idim < out_vol->file->ndims; idim++) {
      out_vol->file->nelements[idim] = in_vol->file->nelements[idim];
      out_vol->file->world_axes[idim] = in_vol->file->world_axes[idim];
   }

   /* Get space for output slice */
   out_vol->volume = NULL;
   out_vol->slice = MALLOC(sizeof(Slice_Data));

   /* Loop through list of axes, getting size of volume and slice */
   total_size = 1;
   for (idim=0; idim < WORLD_NDIMS; idim++) {
      
      /* Get the index for input and output volumes */
      in_vindex = in_vol->file->axes[idim];       /* 0, 1 or 2 */
      out_vindex = defaults.volume_def.axes[idim];    /* 0, 1 or 2 */
      in_findex = in_vol->file->indices[in_vindex];     /* 0 to ndims-1 */
      out_findex = in_vol->file->indices[out_vindex];   /* 0 to ndims-1 */
      size = defaults.volume_def.nelements[idim];

      /* Update output axes and indices and nelements */
      out_vol->file->nelements[out_findex] = size;
      out_vol->file->world_axes[out_findex] = idim;
      out_vol->file->axes[idim] = out_vindex;
      out_vol->file->indices[out_vindex] = out_findex;
      
      /* Update slice size */
      if (out_vindex != 0) {
         out_vol->slice->size[out_vindex-1] = size;
         total_size *= size;
      }
   }
   out_vol->slice->data = MALLOC((size_t) total_size * sizeof(double));

   /* Create the output file */
   create_output_file(outfile,  &defaults.volume_def, 
                      in_vol->file, out_vol->file,
                      history);
   
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : load_vol_info
@INPUT      : names of input filenames
@OUTPUT     : in_vol - description of input volume.
@RETURNS    : (nothing)
@DESCRIPTION: Routine to get information input and 
              output files. Sets up all structures
              completely (including allocating space for data).
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Mon Jun 28 14:16:07 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */
void load_vol_info(char   *infile,
                          MincVolume *in_vol)
{
  int idim, index, ivar, varid;
  int ndims, dim[MAX_VAR_DIMS], imgdim[MAX_VAR_DIMS];
  long size, total_size, total_slice_size;
  File_Info *fp;
  Volume_Definition volume_def;

  /* Get for input volume data information */
  
  in_vol->file = MALLOC(sizeof(File_Info));

  get_file_info(infile, &volume_def, in_vol->file);
  
  in_vol->slice = NULL;
  in_vol->volume = MALLOC(sizeof(Volume_Data));
  in_vol->volume->datatype  = in_vol->file->datatype;
  in_vol->volume->is_signed = in_vol->file->is_signed;
  in_vol->volume->fillvalue = 0.0;
  in_vol->volume->use_fill = TRUE;
  
  /* Get space for volume data */
  total_slice_size = total_size = 1;
  for (idim=0; idim < WORLD_NDIMS; idim++) {
    index = volume_def.axes[idim];
    size = volume_def.nelements[idim];
    total_size *= size;
    in_vol->volume->size[index] = size;
    if (index != 0) 
      total_slice_size *= size;
  }
  
  /* MALLOC OF DATA SPACE */
  
  in_vol->slice = MALLOC(sizeof(Slice_Data));
  in_vol->slice->data = MALLOC((size_t) total_slice_size * sizeof(double));


  /* Get space for slice scale and offset */
  in_vol->volume->scale  = MALLOC(sizeof(double) * in_vol->volume->size[0]);
  in_vol->volume->offset = MALLOC(sizeof(double) * in_vol->volume->size[0]);
  
  /* Check min/max variables */
  fp = in_vol->file;
  if ((fp->maxid != MI_ERROR) && (fp->minid != MI_ERROR)) {
    
    /* Get MIimage dimensions */
    (void) ncvarinq(fp->mincid, fp->imgid, NULL, NULL, &ndims, imgdim, NULL);
    
    /* Check MIimagemax/min dimensions */
    for (ivar=0; ivar<2; ivar++) {
      varid = (ivar==0 ? fp->maxid : fp->minid);
      (void) ncvarinq(fp->mincid, varid, NULL, NULL, &ndims, dim, NULL);
      for (idim=0; idim < ndims; idim++) {
        if ((dim[idim] == imgdim[fp->indices[1]]) ||
            (dim[idim] == imgdim[fp->indices[2]])) {
          (void) fprintf(stderr, 
                         "MIimagemax/min vary over slice dimensions.\n");
          exit(EXIT_FAILURE);
        }
      }        /* End loop over MIimagemax/min dimensions */
    }        /* End loop over variables MIimagemax/min */
    
  }        /* If both MIimagemax/min exist */
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : free_vol_info
@INPUT      : volume
@OUTPUT     : 
@RETURNS    : (nothing)
@DESCRIPTION: Routine to free memory allocated in input_vol_info.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Fri Jun  5 14:20:54 EDT 2009  Claude Lepage
@MODIFIED   : 
---------------------------------------------------------------------------- */
void free_vol_info( MincVolume * vol ) {

  if( vol ) {
    if( vol->slice ) {
      if( vol->slice->data ) FREE( vol->slice->data );
      FREE( vol->slice );
      vol->slice = NULL;
    }
  
    if( vol->volume ) {
      if( vol->volume->scale ) FREE( vol->volume->scale );
      if( vol->volume->offset ) FREE( vol->volume->offset );
      FREE( vol->volume );
      vol->volume = NULL;
    }

    if( vol->file ) {
      miicv_free( vol->file->icvid );
      FREE( vol->file );
      vol->file = NULL;
    }
  }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_file_info
@INPUT      : filename - name of file to read
@OUTPUT     : volume_def - description of volume
              file_info - description of file
@RETURNS    : (nothing)
@DESCRIPTION: Routine to get information about the volume definition of
              a minc file. 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : February 9, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void get_file_info(char              *filename, 
                          Volume_Definition *volume_def,
                          File_Info         *file_info)
{
   int dim[MAX_VAR_DIMS], dimid, status, length;
   int axis_counter, idim, jdim, cur_axis;
   char attstr[MI_MAX_ATTSTR_LEN];
   char dimname[MAX_NC_NAME];
   double vrange[2];

   /* Open the minc file */
   file_info->mincid = ncopen(filename, NC_NOWRITE);
   file_info->name = filename;

   /* Get variable identifiers */
   file_info->imgid = ncvarid(file_info->mincid, MIimage);
   ncopts = 0;
   file_info->maxid = ncvarid(file_info->mincid, MIimagemax);
   file_info->minid = ncvarid(file_info->mincid, MIimagemin);
   ncopts = NC_VERBOSE | NC_FATAL;

   /* Get information about datatype dimensions of variable */
   (void) ncvarinq(file_info->mincid, file_info->imgid, NULL, 
                   &file_info->datatype, &file_info->ndims, dim, NULL);

   /* Get sign attriubte */
   if (file_info->datatype == NC_BYTE)
      file_info->is_signed = FALSE;
   else
      file_info->is_signed = TRUE;
   ncopts = 0;
   if ((miattgetstr(file_info->mincid, file_info->imgid, MIsigntype, 
                    MI_MAX_ATTSTR_LEN, attstr) != NULL)) {
      if (strcmp(attstr, MI_SIGNED) == 0)
         file_info->is_signed = TRUE;
      else if (strcmp(attstr, MI_UNSIGNED) == 0)
         file_info->is_signed = FALSE;
   }
   ncopts = NC_VERBOSE | NC_FATAL;

   /* Get valid max and min */
   ncopts = 0;
   status=miattget(file_info->mincid, file_info->imgid, MIvalid_range, 
                   NC_DOUBLE, 2, vrange, &length);
   if ((status!=MI_ERROR) && (length==2)) {
      if (vrange[1] > vrange[0]) {
         file_info->vrange[0] = vrange[0];
         file_info->vrange[1] = vrange[1];
      }
      else {
         file_info->vrange[0] = vrange[1];
         file_info->vrange[1] = vrange[0];
      }
   }
   else {
      status=miattget1(file_info->mincid, file_info->imgid, MIvalid_max, 
                       NC_DOUBLE, &(file_info->vrange[1]));
      if (status==MI_ERROR)
         file_info->vrange[1] =
            get_default_range(MIvalid_max, file_info->datatype, 
                              file_info->is_signed);
  
      status=miattget1(file_info->mincid, file_info->imgid, MIvalid_min, 
                       NC_DOUBLE, &(file_info->vrange[0]));
      if (status==MI_ERROR)
         file_info->vrange[0] =
            get_default_range(MIvalid_min, file_info->datatype, 
                              file_info->is_signed);
   }
   ncopts = NC_VERBOSE | NC_FATAL;

   /* Set variables for keeping track of spatial dimensions */
   axis_counter = 0;                   /* Keeps track of values for axes */

   /* Initialize volume definition variables */
   for (idim=0; idim < WORLD_NDIMS; idim++) {
      volume_def->axes[idim] = NO_AXIS;
      volume_def->step[idim] = 1.0;
      volume_def->start[idim] = 0.0;
      for (jdim=0; jdim < WORLD_NDIMS; jdim++) {
         if (jdim==idim)
            volume_def->dircos[idim][jdim] = 1.0;
         else
            volume_def->dircos[idim][jdim] = 0.0;
      }
      (void) strcpy(volume_def->units[idim], "mm");
      (void) strcpy(volume_def->spacetype[idim], MI_NATIVE);
   }

   /* Loop through dimensions, getting dimension information */

   for (idim=0; idim < file_info->ndims; idim++) {

      /* Get size of dimension */
      (void) ncdiminq(file_info->mincid, dim[idim], dimname, 
                      &file_info->nelements[idim]);

      /* Check variable name */
      cur_axis = NO_AXIS;
      if (strcmp(dimname, MIxspace)==0)
         cur_axis = X;
      else if (strcmp(dimname, MIyspace)==0)
         cur_axis = Y;
      else if (strcmp(dimname, MIzspace)==0)
         cur_axis = Z;

      /* Save world axis info */
      file_info->world_axes[idim] = cur_axis;

      /* Check for spatial dimension */
      if (cur_axis == NO_AXIS) continue;

      /* Set axis */
      if (volume_def->axes[cur_axis] != NO_AXIS) {
         (void) fprintf(stderr, "Repeated spatial dimension %s in file %s.\n",
                 dimname, filename);
         exit(EXIT_FAILURE);
      }
      volume_def->axes[cur_axis] = axis_counter++;

      /* Save spatial axis specific info */
      file_info->axes[cur_axis] = volume_def->axes[cur_axis];
      file_info->indices[volume_def->axes[cur_axis]] = idim;
      volume_def->nelements[cur_axis] = file_info->nelements[idim];

      /* Check for existence of variable */
      ncopts = 0;
      dimid = ncvarid(file_info->mincid, dimname);
      ncopts = NC_VERBOSE | NC_FATAL;
      if (dimid == MI_ERROR) continue;
             
      /* Get attributes from variable */
      ncopts = 0;
      (void) miattget1(file_info->mincid, dimid, MIstep, 
                       NC_DOUBLE, &volume_def->step[cur_axis]);
      (void) miattget1(file_info->mincid, dimid, MIstart, 
                       NC_DOUBLE, &volume_def->start[cur_axis]);
      (void) miattget(file_info->mincid, dimid, MIdirection_cosines, 
                      NC_DOUBLE, WORLD_NDIMS, 
                      volume_def->dircos[cur_axis], NULL);
      (void) miattgetstr(file_info->mincid, dimid, MIunits, 
                         MI_MAX_ATTSTR_LEN, volume_def->units[cur_axis]);
      (void) miattgetstr(file_info->mincid, dimid, MIspacetype, 
                         MI_MAX_ATTSTR_LEN, volume_def->spacetype[cur_axis]);
      ncopts = NC_VERBOSE | NC_FATAL;

   }   /* End of loop over dimensions */

   /* Check that we have the correct number of spatial dimensions */
   if (axis_counter != WORLD_NDIMS) {
      (void) fprintf(stderr, 
                     "Incorrect number of spatial dimensions in file %s.\n",
                     filename);
         exit(EXIT_FAILURE);
   }

   return;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_output_file
@INPUT      : filename - name of file to create
              volume_def - description of volume
              in_file - description of input file
              out_file - description of output file
              tm_stamp - string describing program invocation
@OUTPUT     : (nothing) 
@RETURNS    : (nothing)
@DESCRIPTION: Routine to create an minc output file and set up its parameters
              properly.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : February 9, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void create_output_file(char *filename,
                               Volume_Definition *volume_def,
                               File_Info *in_file,
                               File_Info *out_file,
                               char *tm_stamp)
{
   int ndims, in_dims[MAX_VAR_DIMS], out_dims[MAX_VAR_DIMS];
   char dimname[MAX_NC_NAME];
   int cur_dim, axis, idim, dimid, varid;
   int att_length;
   nc_type datatype;
   int nexcluded, excluded_vars[10];
   char *string;
   int dim_exists, is_volume_dimension;
   int in_index, out_index;


   axis = 0;
   /* Save the file name */
   out_file->name = filename;

   /* Create the list of excluded variables */
   nexcluded = 0;
   ncopts = 0;
   if ((varid=ncvarid(in_file->mincid, MIxspace)) != MI_ERROR)
      excluded_vars[nexcluded++] = varid;
   if ((varid=ncvarid(in_file->mincid, MIyspace)) != MI_ERROR)
      excluded_vars[nexcluded++] = varid;
   if ((varid=ncvarid(in_file->mincid, MIzspace)) != MI_ERROR)
      excluded_vars[nexcluded++] = varid;
   if ((varid=ncvarid(in_file->mincid, MIimage)) != MI_ERROR)
      excluded_vars[nexcluded++] = varid;
   if ((varid=ncvarid(in_file->mincid, MIimagemax)) != MI_ERROR)
      excluded_vars[nexcluded++] = varid;
   if ((varid=ncvarid(in_file->mincid, MIimagemin)) != MI_ERROR)
      excluded_vars[nexcluded++] = varid;
   ncopts = NC_VERBOSE | NC_FATAL;

   /* Create the file */
   out_file->mincid = nccreate(filename, NC_CLOBBER );
 
   /* Copy all other variable definitions */
   (void) micopy_all_var_defs(in_file->mincid, out_file->mincid, 
                              nexcluded, excluded_vars);

   /* Add the time stamp */
   ncopts=0;
   if ((ncattinq(out_file->mincid, NC_GLOBAL, MIhistory, &datatype,
                 &att_length) == MI_ERROR) ||
       (datatype != NC_CHAR))
      att_length = 0;
   att_length += strlen(tm_stamp) + 1;
   string = MALLOC(att_length);
   string[0] = '\0';
   (void) miattgetstr(out_file->mincid, NC_GLOBAL, MIhistory, att_length, 
                      string);
   ncopts=NC_VERBOSE | NC_FATAL;
   (void) strcat(string, tm_stamp);
   (void) miattputstr(out_file->mincid, NC_GLOBAL, MIhistory, string);
   FREE(string);

   /* Get the dimension ids from the input file */
   (void) ncvarinq(in_file->mincid, in_file->imgid, NULL, NULL, 
                   &ndims, in_dims, NULL);

   /* Check for screw-up on number of dimensions */
   if (ndims != out_file->ndims) {
      (void) fprintf(stderr, 
                     "Error in number of dimensions for output file.\n");
      exit(EXIT_FAILURE);
   }

   /* Create the dimensions for the image variable */
   for (out_index=0; out_index<ndims; out_index++) {

      /* Check to see if this is a volume dimension */
      is_volume_dimension = (out_file->world_axes[out_index] != NO_AXIS);

      /* Get the input index */
      if (!is_volume_dimension) 
         in_index = out_index;
      else {
         axis = out_file->world_axes[out_index];
         if ((axis<0) || (axis>=WORLD_NDIMS)) {
            (void) fprintf(stderr,
                           "Error creating dimensions for output file.\n");
            exit(EXIT_FAILURE);
         }
         in_index = in_file->indices[in_file->axes[axis]];
      }

      /* Get the dimension name from the input file */
      (void) ncdiminq(in_file->mincid, in_dims[in_index], dimname, NULL);

      /* Check to see if the dimension already exists */
      ncopts = 0;
      out_dims[out_index] = ncdimid(out_file->mincid, dimname);
      ncopts = NC_VERBOSE | NC_FATAL;
      dim_exists = (out_dims[out_index] != MI_ERROR);

      /* If we have a volume dimension and it exists already with the wrong
         size, then we must rename it */
      if (is_volume_dimension && dim_exists && 
          (out_file->nelements[out_index] != in_file->nelements[in_index])) {
         string = MALLOC(MAX_NC_NAME);
         ncopts = 0;
         idim = 0;
         do {
            (void) sprintf(string, "%s%d", dimname, idim);
            idim++;
         } while (ncdimid(out_file->mincid, string) != MI_ERROR);
         ncopts = NC_VERBOSE | NC_FATAL;
         (void) ncdimrename(out_file->mincid, out_dims[out_index], string);
         FREE(string);
         out_dims[out_index] = ncdimdef(out_file->mincid, dimname, 
                                        out_file->nelements[out_index]);
      }
      else if (!dim_exists)
         out_dims[out_index] = ncdimdef(out_file->mincid, dimname, 
                                        out_file->nelements[out_index]);

      /* If this is a volume dimension, create a variable */
      if (is_volume_dimension) {

         /* Create the variable */
         dimid = micreate_std_variable(out_file->mincid, dimname, NC_DOUBLE,
                                       0, NULL);
         (void) miattputdbl(out_file->mincid, dimid, MIstep, 
                            volume_def->step[axis]);
         (void) miattputdbl(out_file->mincid, dimid, MIstart, 
                            volume_def->start[axis]);
         (void) ncattput(out_file->mincid, dimid, MIdirection_cosines, 
                         NC_DOUBLE, WORLD_NDIMS, volume_def->dircos[axis]);
         (void) miattputstr(out_file->mincid, dimid, MIunits, 
                            volume_def->units[axis]);
         (void) miattputstr(out_file->mincid, dimid, MIspacetype, 
                            volume_def->spacetype[axis]);

      }       /* If volume dimension */
   }       /* Loop over dimensions */

   /* Create the image variable */
   out_file->imgid = micreate_std_variable(out_file->mincid, MIimage, 
                                           out_file->datatype,
                                           ndims, out_dims);
   (void) micopy_all_atts(in_file->mincid, in_file->imgid,
                          out_file->mincid, out_file->imgid);
   (void) miattputstr(out_file->mincid, out_file->imgid, MIcomplete,
                      MI_FALSE);
   (void) ncattput(out_file->mincid, out_file->imgid, MIvalid_range, 
                   NC_DOUBLE, 2, out_file->vrange);
   if (out_file->is_signed)
      (void) miattputstr(out_file->mincid, out_file->imgid,
                         MIsigntype, MI_SIGNED);
   else
      (void) miattputstr(out_file->mincid, out_file->imgid,
                         MIsigntype, MI_UNSIGNED);

   /* Create the image max and min variables. We have to make sure that
      these vary over spatial slices (this will violate the MINC standard if
      there are non-spatial dimensions that vary faster than the fastest two
      spatial dimensions, but the way the resample code is set up, there's 
      no easy way around it, except to fix it when finishing up) */
   cur_dim = out_file->indices[1];
   for (idim=out_file->indices[1]+1; idim<ndims; idim++) {
      if (idim != out_file->indices[2]) {
         out_dims[cur_dim] = out_dims[idim];
         cur_dim++;
      }
   }
   /* To commit our subterfuge, we check for an error creating the variable
      (it is illegal). If there is an error, we rename MIimage, create the
      variable, restore MIimage and add the pointer */
   ncopts = 0;
   out_file->maxid = micreate_std_variable(out_file->mincid, MIimagemax,
                                           NC_DOUBLE, ndims-2, out_dims);
   ncopts = NC_VERBOSE | NC_FATAL;
   if (out_file->maxid == MI_ERROR) {
      (void) ncvarrename(out_file->mincid, out_file->imgid, TEMP_IMAGE_VAR);
      out_file->maxid = micreate_std_variable(out_file->mincid, MIimagemax,
                                              NC_DOUBLE, ndims-2, out_dims);
      (void) ncvarrename(out_file->mincid, out_file->imgid, MIimage);
      (void) miattput_pointer(out_file->mincid, out_file->imgid, 
                              MIimagemax, out_file->maxid);
   }
   if (in_file->maxid != MI_ERROR)
      (void) micopy_all_atts(in_file->mincid, in_file->maxid,
                             out_file->mincid, out_file->maxid);
   /* Repeat for min variable */
   ncopts = 0;
   out_file->minid = micreate_std_variable(out_file->mincid, MIimagemin,
                                           NC_DOUBLE, ndims-2, out_dims);
   ncopts = NC_VERBOSE | NC_FATAL;
   if (out_file->minid == MI_ERROR) {
      (void) ncvarrename(out_file->mincid, out_file->imgid, TEMP_IMAGE_VAR);
      out_file->minid = micreate_std_variable(out_file->mincid, MIimagemin,
                                              NC_DOUBLE, ndims-2, out_dims);
      (void) ncvarrename(out_file->mincid, out_file->imgid, MIimage);
      (void) miattput_pointer(out_file->mincid, out_file->imgid, 
                              MIimagemin, out_file->minid);
   }
   if (in_file->minid != MI_ERROR)
      (void) micopy_all_atts(in_file->mincid, in_file->minid,
                             out_file->mincid, out_file->minid);

   ncopts = 0;

   /* Get id of processing variable (create it if needed) */
   varid = ncvarid(out_file->mincid, PROCESSING_VAR);
   if (varid == MI_ERROR) {
      varid = ncvardef(out_file->mincid, PROCESSING_VAR, NC_LONG, 0, NULL);
      (void) miadd_child(out_file->mincid, 
                         ncvarid(out_file->mincid, MIrootvariable), varid);
   }

   /* Get into data mode */
   (void) ncendef(out_file->mincid);

   /* Copy all the other data */
   (void) micopy_all_var_values(in_file->mincid, out_file->mincid,
                                nexcluded, excluded_vars);

   /* Create and attach an icv */
   out_file->icvid = miicv_create();
   (void) miicv_setint(out_file->icvid, MI_ICV_TYPE, NC_DOUBLE);
   (void) miicv_setint(out_file->icvid, MI_ICV_DO_NORM, TRUE);
   (void) miicv_setint(out_file->icvid, MI_ICV_USER_NORM, TRUE);
   (void) miicv_attach(out_file->icvid, out_file->mincid, out_file->imgid);

   return;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_default_range
@INPUT      : what     - MIvalid_min means get default min, MIvalid_min means 
                 get default min
              datatype - type of variable
              is_signed   - TRUE if variable is signed
@OUTPUT     : (none)
@RETURNS    : default maximum or minimum for the datatype
@DESCRIPTION: Return the defaults maximum or minimum for a given datatype
              and sign.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 10, 1992 (Peter Neelin)
@MODIFIED   : February 10, 1993 (Peter Neelin)
                 - ripped off from MINC code
---------------------------------------------------------------------------- */
double get_default_range(char *what, nc_type datatype, int is_signed)
{
   double limit;

   limit = 0;

   if (strcmp(what, MIvalid_max)==0) {
      switch (datatype) {
      case NC_LONG:  limit = (is_signed) ? LONG_MAX : ULONG_MAX; break;
      case NC_SHORT: limit = (is_signed) ? SHRT_MAX : USHRT_MAX; break;
      case NC_BYTE:  limit = (is_signed) ? SCHAR_MAX : UCHAR_MAX; break;
      default: limit = DEFAULT_MAX; break;
      }
   }
   else if (strcmp(what, MIvalid_min)==0) {
      switch (datatype) {
      case NC_LONG:  limit = (is_signed) ? LONG_MIN : 0; break;
      case NC_SHORT: limit = (is_signed) ? SHRT_MIN : 0; break;
      case NC_BYTE:  limit = (is_signed) ? SCHAR_MIN : 0; break;
      default: limit = DEFAULT_MIN; break;
      }
   }
   else {
      limit = 0.0;
   }

   return limit;
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : finish_up
@INPUT      : in_vol - input volume
              out_vol - output volume
@OUTPUT     : (nothing) 
@RETURNS    : (nothing)
@DESCRIPTION: Routine to finish up at end of program, closing files, etc.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : February 15, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void finish_up(MincVolume *in_vol1, 
                      MincVolume *in_vol2, 
                      MincVolume *in_vol3, 
                      MincVolume *out_vol)
{
   File_Info *in_file, *out_file;


   /* Close the output file */
   if (out_vol != (MincVolume *)NULL) {
     out_file = out_vol->file;
     (void) miattputstr(out_file->mincid, out_file->imgid, MIcomplete, MI_TRUE);
     (void) ncclose(out_file->mincid);
   }

   /* Close the input files */
   in_file = in_vol1->file;
   (void) ncclose(in_file->mincid);
   in_file = in_vol2->file;
   (void) ncclose(in_file->mincid);
   in_file = in_vol3->file;
   (void) ncclose(in_file->mincid);

   return;
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_gradmag_volumes
@INPUT      : in_vol1 - description of input dx volume
              in_vol2 - description of input dy volume
              in_vol3 - description of input dz volume
              out_vol - description of output dxyz volume
@OUTPUT     : (none)
@RETURNS    : (none)
@DESCRIPTION: dxyz = sqrt(dx*dx + dy*dy + dz*dz), on a voxel by voxel basis.
                this is calculated one slice at a time.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Wed Jun 30 09:01:51 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */
void make_gradmag_volumes(MincVolume *in_vol1, 
                                 MincVolume *in_vol2, 
                                 MincVolume *in_vol3, 
                                 MincVolume *out_vol,
                                 double *min_val, double *max_val)
{
   long in_start[MAX_VAR_DIMS], in_count[MAX_VAR_DIMS], in_end[MAX_VAR_DIMS];
   long out_start[MAX_VAR_DIMS], out_count[MAX_VAR_DIMS];
   long mm_start[MAX_VAR_DIMS];   /* VIO_Vector for min/max variables */
   long nslice, islice;
   int idim, slice_dim, index, slice_index;
   double maximum, minimum, valid_range[2];
   File_Info *ifp,*ofp;

   nslice = slice_index = 0;

   /* Set pointers to file information */
   ifp = in_vol1->file;
   ofp = out_vol->file;

   /* Set input file start, count and end vectors for reading a volume
      at a time */
   (void) miset_coords(ifp->ndims, (long) 0, in_start);
   (void) miset_coords(ifp->ndims, (long) 1, in_count);
   for (idim=0; idim < ifp->ndims; idim++) {
      in_end[idim] = ifp->nelements[idim];
   }
   for (idim=0; idim < VOL_NDIMS; idim++) {
      index = ifp->indices[idim];
      in_count[index] = ifp->nelements[index];
   }

   /* Set output file count for writing a slice and get the number of 
      output slices */
   (void) miset_coords(ifp->ndims, (long) 1, out_count);
   for (idim=0; idim < VOL_NDIMS; idim++) {
      index = ofp->indices[idim];
      if (idim==0) {
         slice_index = index;
         nslice = ofp->nelements[index];
      }
      else {
         out_count[index] = ofp->nelements[index];
      }
   }

   slice_dim = ifp->ndims - 3;

   /* Initialize global max and min */
   valid_range[0] =  DBL_MAX;
   valid_range[1] = -DBL_MAX;

   /* Print log message */
   if (verbose) {
      (void) fprintf(stderr, "Transforming slices:");
      (void) fflush(stderr);
   }


   /* Copy the start vector */
   for (idim=0; idim < ifp->ndims; idim++)
     out_start[idim] = in_start[idim];
   
   /* Loop over slices */
   for (islice=0; islice < nslice; islice++) {
     
     /* Print log message */
     if (verbose) {
       (void) fprintf(stderr, ".");
       (void) fflush(stderr);
     }
     
     /* Read in one slice from each volume */
     in_start[slice_dim] = islice;
     in_count[slice_dim] = 1;
     
     (void) miicv_get(in_vol1->file->icvid, in_start, in_count, in_vol1->slice->data); 
     (void) miicv_get(in_vol2->file->icvid, in_start, in_count, in_vol2->slice->data); 
     (void) miicv_get(in_vol3->file->icvid, in_start, in_count, in_vol3->slice->data); 
     
     /* Set slice number in out_start */
     out_start[slice_index] = islice;
     
     /* Get the slice */
     
     get_mag_slice(out_vol->slice, in_vol1->slice, in_vol2->slice, in_vol3->slice, 
                   &minimum, &maximum);

     /* Update global max and min */
     if (maximum > valid_range[1]) valid_range[1] = maximum;
     if (minimum < valid_range[0]) valid_range[0] = minimum;
     
     /* Write the max, min and slice */
     (void) mivarput1(ofp->mincid, ofp->maxid, 
                      mitranslate_coords(ofp->mincid, 
                                         ofp->imgid, out_start,
                                         ofp->maxid, mm_start),
                      NC_DOUBLE, NULL, &maximum);
     (void) mivarput1(ofp->mincid, ofp->minid, 
                      mitranslate_coords(ofp->mincid, 
                                         ofp->imgid, out_start,
                                         ofp->minid, mm_start),
                      NC_DOUBLE, NULL, &minimum);
     (void) miicv_put(ofp->icvid, out_start, out_count,
                      out_vol->slice->data);
     
   }    /* End loop over slices */
   
   
   /* Print end of log message */
   if (verbose) {
      (void) fprintf(stderr, "Done\n");
      (void) fflush(stderr);
   }

   /* If output volume is floating point, write out global max and min */
   if ((ofp->datatype == NC_FLOAT) || (ofp->datatype == NC_DOUBLE)) {
      (void) ncattput(ofp->mincid, ofp->imgid, MIvalid_range, 
                      NC_DOUBLE, 2, valid_range);
   }

   *min_val = valid_range[0]; 
   *max_val = valid_range[1];

}
/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_curvature_volumes
@INPUT      : in_vol1 - description of input dx volume
              in_vol2 - description of input dy volume
              in_vol3 - description of input dz volume
              out_vol - description of output dxyz volume
@OUTPUT     : (none)
@RETURNS    : (none)
@DESCRIPTION: dxyz = sqrt(dx*dx + dy*dy + dz*dz), on a voxel by voxel basis.
                this is calculated one slice at a time.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Wed Jun 30 09:01:51 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */
void make_curvature_volumes(MincVolume *in_vol1, 
                                   MincVolume *in_vol2, 
                                   MincVolume *in_vol3, 
                                   MincVolume *in_volxx, 
                                   MincVolume *in_volyy, 
                                   MincVolume *in_volzz, 
                                   MincVolume *out_vol,
                                   double thresh)
{
   long in_start[MAX_VAR_DIMS], in_count[MAX_VAR_DIMS], in_end[MAX_VAR_DIMS];
   long out_start[MAX_VAR_DIMS], out_count[MAX_VAR_DIMS];
   long mm_start[MAX_VAR_DIMS];   /* VIO_Vector for min/max variables */
   long nslice, islice;
   int idim, slice_dim, index, slice_index;
   double maximum, minimum, valid_range[2];
   File_Info *ifp,*ofp;

   nslice = slice_index = 0;

   /* Set pointers to file information */
   ifp = in_vol1->file;
   ofp = out_vol->file;

   /* Set input file start, count and end vectors for reading a volume
      at a time */
   (void) miset_coords(ifp->ndims, (long) 0, in_start);
   (void) miset_coords(ifp->ndims, (long) 1, in_count);
   for (idim=0; idim < ifp->ndims; idim++) {
      in_end[idim] = ifp->nelements[idim];
   }
   for (idim=0; idim < VOL_NDIMS; idim++) {
      index = ifp->indices[idim];
      in_count[index] = ifp->nelements[index];
   }

   /* Set output file count for writing a slice and get the number of 
      output slices */
   (void) miset_coords(ifp->ndims, (long) 1, out_count);
   for (idim=0; idim < VOL_NDIMS; idim++) {
      index = ofp->indices[idim];
      if (idim==0) {
         slice_index = index;
         nslice = ofp->nelements[index];
      }
      else {
         out_count[index] = ofp->nelements[index];
      }
   }

   slice_dim = ifp->ndims - 3;

   /* Initialize global max and min */
   valid_range[0] =  DBL_MAX;
   valid_range[1] = -DBL_MAX;

   /* Print log message */
   if (verbose) {
      (void) fprintf(stderr, "Transforming slices:");
      (void) fflush(stderr);
   }


   /* Copy the start vector */
   for (idim=0; idim < ifp->ndims; idim++)
     out_start[idim] = in_start[idim];
   
   /* Loop over slices */
   for (islice=0; islice < nslice; islice++) {
     
     /* Print log message */
     if (verbose) {
       (void) fprintf(stderr, ".");
       (void) fflush(stderr);
     }
     
     /* Read in one slice from each volume */
     in_start[slice_dim] = islice;
     in_count[slice_dim] = 1;
     
     (void) miicv_get(in_vol1->file->icvid, in_start, in_count, in_vol1->slice->data); 
     (void) miicv_get(in_vol2->file->icvid, in_start, in_count, in_vol2->slice->data); 
     (void) miicv_get(in_vol3->file->icvid, in_start, in_count, in_vol3->slice->data); 
     (void) miicv_get(in_volxx->file->icvid, in_start, in_count, in_volxx->slice->data); 
     (void) miicv_get(in_volyy->file->icvid, in_start, in_count, in_volyy->slice->data); 
     (void) miicv_get(in_volzz->file->icvid, in_start, in_count, in_volzz->slice->data); 
     
     /* Set slice number in out_start */
     out_start[slice_index] = islice;
     
     /* Get the slice */
     
     get_curvature_slice(out_vol->slice, in_vol1->slice, in_vol2->slice, in_vol3->slice, 
                         in_volxx->slice, in_volyy->slice, in_volzz->slice, 
                         thresh,&minimum, &maximum);
     
     /* Update global max and min */
     if (maximum > valid_range[1]) valid_range[1] = maximum;
     if (minimum > valid_range[0]) valid_range[0] = minimum;
     
     /* Write the max, min and slice */
     (void) mivarput1(ofp->mincid, ofp->maxid, 
                      mitranslate_coords(ofp->mincid, 
                                         ofp->imgid, out_start,
                                         ofp->maxid, mm_start),
                      NC_DOUBLE, NULL, &maximum);
     (void) mivarput1(ofp->mincid, ofp->minid, 
                      mitranslate_coords(ofp->mincid, 
                                         ofp->imgid, out_start,
                                         ofp->minid, mm_start),
                      NC_DOUBLE, NULL, &minimum);
     (void) miicv_put(ofp->icvid, out_start, out_count,
                      out_vol->slice->data);
     
   }    /* End loop over slices */
   
   
   /* Print end of log message */
   if (verbose) {
      (void) fprintf(stderr, "Done\n");
      (void) fflush(stderr);
   }

   /* If output volume is floating point, write out global max and min */
   if ((ofp->datatype == NC_FLOAT) || (ofp->datatype == NC_DOUBLE)) {
      (void) ncattput(ofp->mincid, ofp->imgid, MIvalid_range, 
                      NC_DOUBLE, 2, valid_range);
   }

}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_mag_slice
@INPUT      : slice_dx - slice data struct contains partial dirivitive in x dir
              slice_dy - slice data struct contains partial dirivitive in y dir
              slice_dz - slice data struct contains partial dirivitive in z dir
@OUTPUT     : result - slice data struct contains new slice
              minimum - slice minimum (excluding data from outside volume)
              maximum - slice maximum (excluding data from outside volume)
@RETURNS    : (none)
@DESCRIPTION: result = sqrt( dx^2 + dy^2 + dz^2) for each pixel.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Tue Jun 29 10:46:01 EST 1993 (louis collins)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void get_mag_slice(Slice_Data *result,
                          Slice_Data *slice_dx,
                          Slice_Data *slice_dy,
                          Slice_Data *slice_dz,
                          double *minimum, double *maximum)
{
   double *dptr;
   double *dptr1;
   double *dptr2;
   double *dptr3;
   long irow, icol;
   
   /* Get slice and volume pointers */
   
   dptr  = result->data;
   dptr1 = slice_dx->data;
   dptr2 = slice_dy->data;
   dptr3 = slice_dz->data;

   *maximum = -DBL_MAX;
   *minimum =  DBL_MAX;

   for (irow=0; irow < result->size[0]; irow++) {   /* Loop over rows of slice */
     
     for (icol=0; icol < result->size[1]; icol++) {      /* Loop over columns */
       
       *dptr = sqrt( (*dptr1 * *dptr1) + (*dptr2 * *dptr2) + (*dptr3 * *dptr3) );

       if (*dptr > *maximum) *maximum = *dptr;
       if (*dptr < *minimum) *minimum = *dptr;

       dptr++;
       dptr1++;
       dptr2++;
       dptr3++;
       
     }     /* Loop over columns */
   }        /* Loop over rows */



   if ((*maximum == -DBL_MAX) && (*minimum ==  DBL_MAX)) {
     *minimum = 0.0;
     *maximum = SMALL_VALUE;
   }
   else if (*maximum <= *minimum) {
     if (*minimum == 0.0) 
       *maximum = SMALL_VALUE;
     else if (*minimum < 0.0)
       *maximum = 0.0;
     else
       *maximum = 2.0 * (*minimum);
   }

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_curvature_slice
@INPUT      : slice_dx - slice data struct contains 1st partial dirivitive in x dir
              slice_dy - slice data struct contains 1st partial dirivitive in y dir
              slice_dz - slice data struct contains 1st partial dirivitive in z dir
              slice_dxx - slice data struct contains 2nd partial dirivitive in x dir
              slice_dyy - slice data struct contains 2nd partial dirivitive in y dir
              slice_dzz - slice data struct contains 2nd partial dirivitive in z dir
@OUTPUT     : result - slice data struct contains new slice
              minimum - slice minimum (excluding data from outside volume)
              maximum - slice maximum (excluding data from outside volume)
@RETURNS    : (none)
@DESCRIPTION: result = (dx^2)*dyy*dzz + (dy^2)*dxx*dzz + (dz^2)*dxx*dyy
                       ------------------------------------------------
                                     (dx^2+ dy^2 + dz^2)^2
                   for each pixel.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Apr 6, 1994 lc
@MODIFIED   : 
---------------------------------------------------------------------------- */
void get_curvature_slice(Slice_Data *result,
                                Slice_Data *slice_dx,
                                Slice_Data *slice_dy,
                                Slice_Data *slice_dz,
                                Slice_Data *slice_dxx,
                                Slice_Data *slice_dyy,
                                Slice_Data *slice_dzz,
                                double thresh,
                                double *minimum, double *maximum)
{
   double *dptr;
   double *dptrx;
   double *dptry;
   double *dptrz;
   double *dptrxx;
   double *dptryy;
   double *dptrzz;
   double t,numerator;
   long irow, icol;
   
   /* Get slice and volume pointers */
   
   dptr  = result->data;
   dptrx = slice_dx->data;  dptry =slice_dy->data;  dptrz =slice_dz->data;
   dptrxx= slice_dxx->data; dptryy=slice_dyy->data; dptrzz=slice_dzz->data;

   *maximum = -DBL_MAX;
   *minimum =  DBL_MAX;

   for (irow=0; irow < result->size[0]; irow++) {   /* Loop over rows of slice */
     
     for (icol=0; icol < result->size[1]; icol++) {      /* Loop over columns */
       
       t = (*dptrx * *dptrx) + (*dptry * *dptry) + (*dptrz * *dptrz);

       if (t < thresh*thresh) {
         *dptr = 0.0;
       }
       else {

         numerator = (*dptrx * *dptrx) * *dptryy * *dptrzz + 
                     (*dptry * *dptry) * *dptrxx * *dptrzz +
                     (*dptrz * *dptrz) * *dptrxx * *dptryy;

         *dptr = numerator / (t*t);

       }

       if (*dptr > *maximum) *maximum = *dptr;
       if (*dptr < *minimum) *minimum = *dptr;

       dptr++;
       dptrx++;  dptry++;  dptrz++;
       dptrxx++; dptryy++; dptrzz++;
       
     }     /* Loop over columns */
   }        /* Loop over rows */
   

   if ((*maximum == -DBL_MAX) && (*minimum ==  DBL_MAX)) {
     *minimum = 0.0;
     *maximum = SMALL_VALUE;
   }
   else if (*maximum <= *minimum) {
     if (*minimum == 0.0) 
       *maximum = SMALL_VALUE;
     else if (*minimum < 0.0)
       *maximum = 0.0;
     else
       *maximum = 2.0 * (*minimum);
   }

}

