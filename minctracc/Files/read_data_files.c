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
@MODIFIED   : Revision 1.4  2004-02-12 05:54:11  rotor
@MODIFIED   :  * removed public/private defs
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
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Files/read_data_files.c,v 1.4 2004-02-12 05:54:11 rotor Exp $";
#endif

#include <volume_io/internal_volume_io.h>
#include <print_error.h>

static char *default_dim_names[N_DIMENSIONS] =
   { MIzspace, MIyspace, MIxspace };



Status read_deform_data(Volume *dx,
			Volume *dy,
			Volume *dz,
			char *name)
{
  char fullname[500];
  Status status;
  Volume tx, ty, tz;
  
  status = OK;
  
  
  /* ------------------------------------------------------------ dx vol --*/
  (void)sprintf(fullname,"%s_dx.mnc",name);
  
  if (!file_exists(fullname))
    print_error_and_line_num("Cannot find %s\n",__FILE__, __LINE__, fullname);
  
  status = input_volume( fullname, 3, default_dim_names, 
			NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			TRUE, &tx, (minc_input_options *)NULL );
  if (status != OK)
    print_error_and_line_num("problems reading in dx volume, probably not enough memory!\n",__FILE__, __LINE__);
  
  *dx = tx;
  
  
  
  /* ------------------------------------------------------------ dy vol --*/
  (void)sprintf(fullname,"%s_dy.mnc",name);
  
  if (!file_exists(fullname))
    print_error_and_line_num("Cannot find %s\n",__FILE__, __LINE__, fullname);
  
  status = input_volume( fullname, 3, default_dim_names, 
			NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			TRUE, &ty, (minc_input_options *)NULL );
  if (status != OK)
    print_error_and_line_num("problems reading in dy volume, probably not enough memory!\n",__FILE__, __LINE__);
  
  *dy = ty;
  
  
  
  /* ------------------------------------------------------------ dz vol --*/
  (void)sprintf(fullname,"%s_dz.mnc",name);
  
  if (!file_exists(fullname))
    print_error_and_line_num("Cannot find %s\n",__FILE__, __LINE__, fullname);
  
  status = input_volume( fullname, 3, default_dim_names, 
			NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			TRUE, &tz, (minc_input_options *)NULL );
  if (status != OK)
    print_error_and_line_num("problems reading in dz volume, probably not enough memory!\n",__FILE__, __LINE__);
  
  *dz = tz;
    
  return(status);
}


Status read_all_data(Volume *dblur,
			     Volume *dx,
			     Volume *dy,
			     Volume *dz,
			     Volume *dxyz, 
			     char *name)
{
  char fullname[500];
  Status status;
  Volume tx, ty, tz, txyz, tblur;
  
  status = OK;
  
  status = read_deform_data(&tx, &ty, &tz, name);
  if (status != OK) {
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
    if (status != OK)
      print_error_and_line_num("problems reading in dxyz volume, maybe not enough memory!\n",__FILE__, __LINE__);
    

    *dxyz = txyz;


    (void)sprintf(fullname,"%s_blur.mnc",name);
    
    if (!file_exists(fullname))
      print_error_and_line_num("Cannot find %s\n",__FILE__, __LINE__, fullname);
    
    status = input_volume( fullname, 3, default_dim_names, 
			  NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			  TRUE, &tblur, (minc_input_options *)NULL );
    if (status != OK)
      print_error_and_line_num("problems reading in dxyz volume, maybe not enough memory!\n",__FILE__, __LINE__);
    
    *dblur = tblur;
  }
  
  return(status);
  
}


Status save_deform_data(Volume dx,
			       Volume dy,
			       Volume dz,
			       char *name,
			       char *history)
{
  char fullname[500];
  Status status;
  
  status = OK;
  
  
  /* ------------------------------------------------------------ dx vol --*/
  (void)sprintf(fullname,"%s_dx.mnc",name);
  
/*  if (!file_exists(fullname))
    print_error_and_line_num("Cannot find %s\n",__FILE__, __LINE__, fullname);
  */

  status = output_volume(fullname, 
			 NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			 dx, history, (minc_output_options *)NULL );
  if (status != OK)
    print_error_and_line_num("problems saving in dx volume.\n",__FILE__, __LINE__);
    
  
  /* ------------------------------------------------------------ dy vol --*/
  (void)sprintf(fullname,"%s_dy.mnc",name);
  
/*  if (!file_exists(fullname))
    print_error_and_line_num("Cannot find %s\n",__FILE__, __LINE__, fullname);
  */

  status = output_volume( fullname, 
			NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			dy, history, (minc_output_options *)NULL );
  if (status != OK)
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
  if (status != OK)
    print_error_and_line_num("problems saving in dz volume.\n",__FILE__, __LINE__);
  
  return(status);
}

