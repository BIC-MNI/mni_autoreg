/* ----------------------------- MNI Header -----------------------------------
@NAME       : mincbbox.c
@INPUT      : minc volume, 
@OUTPUT     : bounding box coordinates
@RETURNS    : TRUE if ok, ERROR if error.
@DESCRIPTION: This program will calculate the bounding box of the data in
              the minc volume and return the coordinates of the box in terms
	      of startx, starty, startz, widthx, widthy, widthz 
	      all in mm
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

@CREATED    : Thu Jun  2 10:21:00 EST 1994   Louis Collins
@MODIFIED   : $Log: mincbbox.c,v $
@MODIFIED   : Revision 1.1  2000-01-27 16:40:02  louis
@MODIFIED   : Initial revision
@MODIFIED   :
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/mincbbox/mincbbox.c,v 1.1 2000-01-27 16:40:02 louis Exp $";
#endif



#include <config.h>
#include <volume_io.h>
#include <Proglib.h>
#include "mincbbox.h"


static char *default_dim_names[N_DIMENSIONS] = { MIxspace, MIyspace, MIzspace };
static char *My_File_order_dimension_names[MAX_DIMENSIONS] = { "", "", "", "", "" };

main (int argc, char *argv[] )
{   
  FILE 
    *ofd;
  char 
    *infilename;
  Status 
    status;
  Volume
    data;
  Real
    voxel, value, 
    x,y,z,
    minx, miny, minz,
    maxx, maxy, maxz;
  int
    v_minx, v_miny, v_minz,
    v_maxx, v_maxy, v_maxz;

  int
    i,j,k, 
    sizes[3];

  /* set default values */
  
  prog_name = argv[0];
  infilename =  NULL;

  verbose = TRUE;
  debug   = FALSE;

  /* Call ParseArgv to interpret all command line args */

  if (ParseArgv(&argc, argv, argTable, 0) || (argc!=2)) {
    (void) fprintf(stderr, 
		   "\nUsage: %s [<options>] <inputfile> \n", 
		   prog_name);
    (void) fprintf(stderr,"       %s [-help]\n\n", prog_name);
    exit(EXIT_FAILURE);
  }

  infilename  = argv[1];	/* set up necessary file names */

  if(debug) {
    print ("input -> %s\n",infilename);
    print ("thres -> %f\n",threshold);
  }

  /******************************************************************************/
  /*             read in input volume data                                      */
  /******************************************************************************/
  
  if (mincreshape)
      status = input_volume(infilename, 3, My_File_order_dimension_names , 
			    NC_UNSPECIFIED, FALSE, 
			    0.0, 0.0, TRUE, &data, 
			    (minc_input_options *)NULL);
  else
      status = input_volume(infilename, 3, default_dim_names , 
			    NC_UNSPECIFIED, FALSE, 
			    0.0, 0.0, TRUE, &data, 
			    (minc_input_options *)NULL);


  if ( status != OK )
    print_error_and_line_num("problems reading `%s'.\n",__FILE__, __LINE__,infilename, 0,0,0,0);

  if (get_volume_n_dimensions(data)!=3) {
    print_error_and_line_num ("File %s has %d dimensions.  Only 3 dims supported.", 
		 __FILE__, __LINE__, infilename, get_volume_n_dimensions(data),0,0,0);
  }
    
  get_volume_sizes(data, sizes);

  minx = miny = minz = DBL_MAX;
  maxx = maxy = maxz = -DBL_MAX;

  v_minx = v_miny = v_minz = INT_MAX;
  v_maxx = v_maxy = v_maxz = -INT_MAX;

  for_less(i, 0, sizes[0]) {
    for_less(j, 0, sizes[1]) {
      for_less(k, 0, sizes[2]) {
	GET_VOXEL_3D(voxel, data, i,j,k);
	value = CONVERT_VOXEL_TO_VALUE(data, voxel);

	if (value>threshold) {

	  convert_3D_voxel_to_world(data,(Real)i,(Real)j,(Real)k, &x,&y,&z);
	  
	  if (x<minx) minx = x;
	  if (y<miny) miny = y;
	  if (z<minz) minz = z;

	  if (x>maxx) maxx = x;
	  if (y>maxy) maxy = y;
	  if (z>maxz) maxz = z;

	  if (i<v_minx) v_minx = i;
	  if (j<v_miny) v_miny = j;
	  if (k<v_minz) v_minz = k;

	  if (i>v_maxx) v_maxx = i;
	  if (j>v_maxy) v_maxy = j;
	  if (k>v_maxz) v_maxz = k;

	}
      }
    }
  }

  
  if (minccrop) {
    print ("-xlim %d %d -ylim %d %d -zlim %d %d\n",
           v_minx, v_maxx, v_miny,v_maxy, v_minz,v_maxz);
  }
  else
  if (mincreshape) {
    print ("-start %d,%d,%d -count %d,%d,%d\n",
	   v_minx,v_miny,v_minz,
	   v_maxx-v_minx+1,v_maxy-v_miny+1,v_maxz-v_minz+1);
  }
  else
  if (mincresample)
    print ("-step 1.0 1.0 1.0 -start %f %f %f -nelements %d %d %d\n",
	   minx, miny, minz, ROUND(maxx-minx)+1, ROUND(maxy-miny)+1, ROUND(maxz-minz)+1);
  else
    if (two_lines)
      print ("%f %f %f\n%f %f %f\n",minx, miny, minz, maxx-minx+1, maxy-miny+1, maxz-minz+1);
    else
      print ("%f %f %f    %f %f %f\n",minx, miny, minz, maxx-minx+1, maxy-miny+1, maxz-minz+1);

  return(OK);
}


