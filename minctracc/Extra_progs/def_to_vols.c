/* ----------------------------- MNI Header -----------------------------------
@NAME       : def_to_vols
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : status
@DESCRIPTION: Program to extract the 3 componants of a deformation vol.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Mon May  8 11:29:40 MET DST 1995  LC
@MODIFIED   : $Log: def_to_vols.c,v $
@MODIFIED   : Revision 1.5  2006-11-29 09:09:31  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 1.4  2005/07/20 20:45:46  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 1.3  2004/02/12 05:54:05  rotor
@MODIFIED   :  * removed public/private defs
@MODIFIED   :
@MODIFIED   : Revision 1.2  2002/03/26 14:15:29  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 1.1  1999/10/25 19:52:06  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :

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
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Extra_progs/def_to_vols.c,v 1.5 2006-11-29 09:09:31 rotor Exp $";
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <volume_io.h>

/* Constants */
#ifndef TRUE
#  define TRUE 1
#  define FALSE 0
#endif

static char *my_ZYX_dim_names[] = { MIzspace, MIyspace, MIxspace };

void get_volume_XYZV_indices(VIO_Volume data, int xyzv[])
{
  
  int 
    axis, i, vol_dims;
  char 
    **data_dim_names;

  vol_dims       = get_volume_n_dimensions(data);
  data_dim_names = get_volume_dimension_names(data);

  for(i=0; i<N_DIMENSIONS+1; i++) xyzv[i] = -1;
  for(i=0; i<vol_dims; i++) {
    if (convert_dim_name_to_spatial_axis(data_dim_names[i], &axis )) {
      xyzv[axis] = i; 
    } 
    else {     /* not a spatial axis */
      xyzv[Z+1] = i;
    }
  }
  delete_dimension_names(data_dim_names);

}


/* Main program */

int main(int argc, char *argv[])
{
   VIO_Volume
     component_vol,
     current_vol;
   FILE 
     *fp;
   VIO_Real
     min_range, max_range,
     value,
     voxel[VIO_MAX_DIMENSIONS],
     steps[VIO_MAX_DIMENSIONS],
     new_steps[VIO_MAX_DIMENSIONS],
     start[VIO_MAX_DIMENSIONS],
     new_start[VIO_MAX_DIMENSIONS];
   int
     i,
     ind[VIO_MAX_DIMENSIONS],
     count[VIO_MAX_DIMENSIONS],
     new_count[VIO_MAX_DIMENSIONS],
     new_xyzv[VIO_MAX_DIMENSIONS],
     xyzv[VIO_MAX_DIMENSIONS];
   
   char name[500];

   VIO_Status status;

   /* Check arguments */
   if (argc != 3 && argc != 5) {
      (void) fprintf(stderr, "Usage: %s <input.mnc> <output_basename> [min_real_range max_real_range]\n",
                     argv[0]);
      exit(EXIT_FAILURE);
   }

   if( input_volume( argv[1], 4, NULL, NC_UNSPECIFIED, FALSE,
                    0.0, 0.0, TRUE, &current_vol, (minc_input_options *)NULL ) != OK ) {
     (void)fprintf(stderr, "Error opening input volume file %s.\n",
                   argv[1]);
     exit(EXIT_FAILURE);
   }

   if (argc==5) {
     min_range = atof(argv[3]);
     max_range = atof(argv[4]);
   }
   else {
     min_range = -50.0;
     max_range = 50;
   }
   get_volume_sizes(       current_vol, count);
   get_volume_separations( current_vol, steps);
   get_volume_XYZV_indices(current_vol, xyzv);
   
   component_vol = create_volume(3, my_ZYX_dim_names, NC_SHORT, TRUE, 0.0, 0.0);
   get_volume_XYZV_indices(component_vol, new_xyzv);

   for(i=0; i<3; i++) {
     new_count[ new_xyzv[i] ] = count[ xyzv[i] ];
     new_steps[ new_xyzv[i] ] = steps[ xyzv[i] ];
   }
   for(i=0; i<MAX_DIMENSIONS; i++) {
     start[i] = 0.0;
     voxel[i] = 0.0;
     new_start[i] = 0.0;
   }

   convert_voxel_to_world(current_vol,
                          voxel,
                          &start[0], &start[1], &start[2]);
   
   set_volume_sizes(component_vol, new_count);
   set_volume_separations(component_vol, new_steps);
   set_volume_starts(component_vol,  start);
   alloc_volume_data(component_vol);
   set_volume_real_range(component_vol, min_range, max_range);

   
   for(ind[xyzv[Z+1]]=0; ind[xyzv[Z+1]]<count[xyzv[Z+1]]; ind[xyzv[Z+1]]++) {

     for(ind[xyzv[VIO_X]]=0; ind[xyzv[VIO_X]]<count[xyzv[VIO_X]]; ind[xyzv[VIO_X]]++) 
       for(ind[xyzv[VIO_Y]]=0; ind[xyzv[VIO_Y]]<count[xyzv[VIO_Y]]; ind[xyzv[VIO_Y]]++) 
         for(ind[xyzv[VIO_Z]]=0; ind[xyzv[VIO_Z]]<count[xyzv[VIO_Z]]; ind[xyzv[VIO_Z]]++) {
           
           value = get_volume_real_value(current_vol,
                                         ind[0],ind[1],ind[2],ind[3],ind[4]);
           set_volume_real_value(component_vol,
                               ind[ xyzv[VIO_Z] ], ind[ xyzv[VIO_Y] ], ind[ xyzv[VIO_X] ], 0, 0,
                                 value);
           
         }
     
     switch (ind[ xyzv[Z+1] ]) {
     case X: (void)sprintf(name,"%s_dx.mnc",argv[2]); break;
     case Y: (void)sprintf(name,"%s_dy.mnc",argv[2]); break;
     case Z: (void)sprintf(name,"%s_dz.mnc",argv[2]); break;
     }
     if (output_volume(name, NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                       component_vol,
                       NULL, (minc_output_options *)NULL) != OK) {
       
       (void) fprintf(stderr,"Cannot write %s\n",name);
       exit( EXIT_FAILURE );
       
       }

   }
    
   delete_volume(component_vol);

   exit(EXIT_SUCCESS);
}
