/* ----------------------------- MNI Header -----------------------------------
@NAME       : read_data_files.c

@DESCRIPTION: routines to read in volumetric data files, using a basename
              for the filename as input and appending extensions to it
              to get the necessary names.
@COPYRIGHT  :
              Copyright 1993 Louis Collins, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.

@CREATED    : Tue Nov 16 13:56:34 EST 1993 LC
@MODIFIED   : $Log: read_data_files.c,v $
@MODIFIED   : Revision 1.6  2006-11-29 09:09:32  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 1.5  2005/07/20 20:45:48  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 1.4  2004/02/12 05:54:11  rotor
@MODIFIED   :  * removed /static defs
@MODIFIED   :
@MODIFIED   : Revision 1.3  2004/02/04 20:42:02  lenezet
@MODIFIED   : *** empty log message ***
@MODIFIED   :
@MODIFIED   : Revision 1.2  2002/03/26 14:15:36  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 1.1  1999/10/25 19:52:13  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :
 * Revision 1.3  1995/02/22  08:56:06  collins
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.2  94/04/06  11:48:48  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D 
 * 
 * Revision 1.1  94/02/21  16:36:21  louis
 * Initial revision
 * 
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/minctracc/Files/read_data_files.c,v 1.6 2006-11-29 09:09:32 rotor Exp $";
#endif

#include <config.h>
#include <volume_io.h>
#include <Proglib.h>

static char *default_dim_names[VIO_N_DIMENSIONS] =
   { MIzspace, MIyspace, MIxspace };



VIO_Status read_deform_data(VIO_Volume *dx,
                        VIO_Volume *dy,
                        VIO_Volume *dz,
                        char *name)
{
  char fullname[500];
  VIO_Status status;
  VIO_Volume tx, ty, tz;
  
  status = VIO_OK;
  
  
  /* ------------------------------------------------------------ dx vol --*/
  (void)sprintf(fullname,"%s_dx.mnc",name);
  
  if (!file_exists(fullname))
    print_error_and_line_num("Cannot find %s\n",__FILE__, __LINE__, fullname);
  
  status = input_volume( fullname, 3, default_dim_names, 
                        NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                        TRUE, &tx, (minc_input_options *)NULL );
  if (status != VIO_OK)
    print_error_and_line_num("problems reading in dx volume, probably not enough memory!\n",__FILE__, __LINE__);
  
  *dx = tx;
  
  
  
  /* ------------------------------------------------------------ dy vol --*/
  (void)sprintf(fullname,"%s_dy.mnc",name);
  
  if (!file_exists(fullname))
    print_error_and_line_num("Cannot find %s\n",__FILE__, __LINE__, fullname);
  
  status = input_volume( fullname, 3, default_dim_names, 
                        NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                        TRUE, &ty, (minc_input_options *)NULL );
  if (status != VIO_OK)
    print_error_and_line_num("problems reading in dy volume, probably not enough memory!\n",__FILE__, __LINE__);
  
  *dy = ty;
  
  
  
  /* ------------------------------------------------------------ dz vol --*/
  (void)sprintf(fullname,"%s_dz.mnc",name);
  
  if (!file_exists(fullname))
    print_error_and_line_num("Cannot find %s\n",__FILE__, __LINE__, fullname);
  
  status = input_volume( fullname, 3, default_dim_names, 
                        NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                        TRUE, &tz, (minc_input_options *)NULL );
  if (status != VIO_OK)
    print_error_and_line_num("problems reading in dz volume, probably not enough memory!\n",__FILE__, __LINE__);
  
  *dz = tz;
    
  return(status);
}


VIO_Status read_all_data(VIO_Volume *dblur,
                             VIO_Volume *dx,
                             VIO_Volume *dy,
                             VIO_Volume *dz,
                             VIO_Volume *dxyz, 
                             char *name)
{
  char fullname[500];
  VIO_Status status;
  VIO_Volume tx, ty, tz, txyz, tblur;
  
  status = VIO_OK;
  
  status = read_deform_data(&tx, &ty, &tz, name);
  if (status != VIO_OK) {
    print_error_and_line_num("problems reading in dx,dy or dz volume.\n",__FILE__, __LINE__);
  }
  else {
    *dx = tx;
    *dy = ty;
    *dz = tz;
    
    
    (void)sprintf(fullname,"%s_dxyz.mnc",name);
    
    if (!file_exists(fullname))
      print_error_and_line_num("Cannot find %s\n",__FILE__, __LINE__, fullname);
    
    status = input_volume( fullname, 3, default_dim_names, 
                          NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                          TRUE, &txyz, (minc_input_options *)NULL );
    if (status != VIO_OK)
      print_error_and_line_num("problems reading in dxyz volume, maybe not enough memory!\n",__FILE__, __LINE__);
    

    *dxyz = txyz;


    (void)sprintf(fullname,"%s_blur.mnc",name);
    
    if (!file_exists(fullname))
      print_error_and_line_num("Cannot find %s\n",__FILE__, __LINE__, fullname);
    
    status = input_volume( fullname, 3, default_dim_names, 
                          NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                          TRUE, &tblur, (minc_input_options *)NULL );
    if (status != VIO_OK)
      print_error_and_line_num("problems reading in dxyz volume, maybe not enough memory!\n",__FILE__, __LINE__);
    
    *dblur = tblur;
  }
  
  return(status);
  
}


VIO_Status save_deform_data(VIO_Volume dx,
                               VIO_Volume dy,
                               VIO_Volume dz,
                               char *name,
                               char *history)
{
  char fullname[500];
  VIO_Status status;
  
  status = VIO_OK;
  
  
  /* ------------------------------------------------------------ dx vol --*/
  (void)sprintf(fullname,"%s_dx.mnc",name);
  
/*  if (!file_exists(fullname))
    print_error_and_line_num("Cannot find %s\n",__FILE__, __LINE__, fullname);
  */

  status = output_volume(fullname, 
                         NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                         dx, history, (minc_output_options *)NULL );
  if (status != VIO_OK)
    print_error_and_line_num("problems saving in dx volume.\n",__FILE__, __LINE__);
    
  
  /* ------------------------------------------------------------ dy vol --*/
  (void)sprintf(fullname,"%s_dy.mnc",name);
  
/*  if (!file_exists(fullname))
    print_error_and_line_num("Cannot find %s\n",__FILE__, __LINE__, fullname);
  */

  status = output_volume( fullname, 
                        NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                        dy, history, (minc_output_options *)NULL );
  if (status != VIO_OK)
    print_error_and_line_num("problems saving in dy volume.\n",__FILE__, __LINE__);
  
  
  /* ------------------------------------------------------------ dz vol --*/
  (void)sprintf(fullname,"%s_dz.mnc",name);
  
/*
  if (!file_exists(fullname))
    print_error_and_line_num("Cannot find %s\n",__FILE__, __LINE__, fullname);
*/
  
  status = output_volume( fullname, 
                        NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                        dz, history, (minc_output_options *)NULL );
  if (status != VIO_OK)
    print_error_and_line_num("problems saving in dz volume.\n",__FILE__, __LINE__);
  
  return(status);
}

