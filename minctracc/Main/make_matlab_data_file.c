/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_matlab_data_file.c
@INPUT      : d1,d2,m1,m2 - data and mask volumes.
@OUTPUT     : creates an output file readable by matlab.
@RETURNS    : nothing
@DESCRIPTION: 
              this routines calculates the objective function value for 
	      the current transformation and variants thereof.

	      each parameter is varied in turn, one at a time, from 
	      -simplex to +simplex around the parameter.
	      
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

@CREATED    : Mon Oct  4 13:06:17 EST 1993 Louis
@MODIFIED   : $Log: make_matlab_data_file.c,v $
@MODIFIED   : Revision 1.2  1994-02-21 16:35:40  louis
@MODIFIED   : version before feb 22 changes
@MODIFIED   :
 * Revision 1.1  93/11/15  16:26:47  louis
 * Initial revision
 * 
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Main/make_matlab_data_file.c,v 1.2 1994-02-21 16:35:40 louis Exp $";
#endif


#include <volume_io.h>
#include <recipes.h>
#include <limits.h>

#include "constants.h"
#include "arg_data.h"
#include "objectives.h"
#include "segment_table.h"
#include "print_error.h"

extern Arg_Data main_args;

extern   Volume   Gdata1, Gdata2, Gmask1, Gmask2;
extern   int      Ginverse_mapping_flag, Gndim;
extern   double   simplex_size ;
extern   Segment_Table  *segment_table;

public float fit_function(float *params);

public void make_zscore_volume(Volume d1, Volume m1, 
			       float threshold); 

public void add_speckle_to_volume(Volume d1, 
				  float speckle,
				  double  *start, int *count, VectorR directions[]);

public void make_matlab_data_file(Volume d1,
				  Volume d2,
				  Volume m1,
				  Volume m2, 
				  char *comments,
				  Arg_Data *globals)
{

#define NUM_STEPS 15

  Status
    status;
  float 
    *p;
  int 
    i,j, 
    ndim;
  FILE
    *ofd;
  Real
    start,step;

  if (globals->obj_function == zscore_objective) { /* replace volume d1 and d2 by zscore volume  */
    make_zscore_volume(d1,m1,(float)globals->threshold[0]);
    make_zscore_volume(d2,m2,(float)globals->threshold[1]);
  } 
  else  if (globals->obj_function == ssc_objective) {	/* add speckle to the data set */
    
    make_zscore_volume(d1,m1,(float)globals->threshold[0]); /* need to make data sets comparable */
    make_zscore_volume(d2,m2,(float)globals->threshold[1]); /* in mean and sd...                 */
    
    if (globals->smallest_vol == 1)
      add_speckle_to_volume(d1, 
			    globals->speckle,
			    globals->start, globals->count, globals->directions);
    else
      add_speckle_to_volume(d2, 
			    globals->speckle,
			    globals->start, globals->count, globals->directions);    
  } else if (globals->obj_function == vr_objective) {
    if (globals->smallest_vol == 1) {
      if (!build_segment_table(&segment_table, d1, globals->groups))
	print_error("%s",__FILE__, __LINE__,"Could not build segment table for source volume\n");
    }
    else {
      if (!build_segment_table(&segment_table, d2, globals->groups))
	print_error("%s",__FILE__, __LINE__,"Could not build segment table for target volume\n");
    }
    if (globals->flags.debug && globals->flags.verbose>1) {
      print ("groups = %d\n",segment_table->groups);
      for_less(i, segment_table->min, segment_table->max+1) {
	print ("%5d: table = %5d, function = %5d\n",i,segment_table->table[i],
	       (segment_table->segment)(i,segment_table) );
      }
    }
    
  }

  switch (globals->trans_info.transform_type) {
  case TRANS_PROCRUSTES: 
    ndim = 7;
    break;
  case TRANS_LSQ: 
    ndim = 12;
    break;
  case TRANS_LSQ6: 
    ndim = 6;
    break;
  case TRANS_LSQ7: 
    ndim = 7;
    break;
  case TRANS_LSQ9: 
    ndim = 9;
    break;
  case TRANS_LSQ12: 
    ndim = 12;
    break;
  default:
    (void)fprintf(stderr, "Unknown type of transformation requested (%d)\n",
		   globals->trans_info.transform_type);
    (void)fprintf(stderr, "Error in line %d, file %s\n",__LINE__, __FILE__);
    print_error   ("%s",__FILE__, __LINE__,"Exiting to system.\n");
  }

				/* set GLOBALS to communicate with the
				   function to be fitted!              */
  Gndim = ndim;
  if (globals->smallest_vol == 1) {
    Gdata1 = d1;
    Gdata2 = d2;
    Gmask1 = m1;
    Gmask2 = m2;
    Ginverse_mapping_flag = FALSE;
  }
  else {
    Gdata1 = d2;
    Gdata2 = d1;
    Gmask1 = m2;
    Gmask2 = m1;
    Ginverse_mapping_flag = TRUE;
  }



  p = vector(1,ndim); /* parameter values */
  
  p[1]=globals->trans_info.translations[0];
  p[2]=globals->trans_info.translations[1];
  p[3]=globals->trans_info.translations[2];
  
  p[4]=globals->trans_info.rotations[0]; 
  p[5]=globals->trans_info.rotations[1]; 
  p[6]=globals->trans_info.rotations[2];
  
  if (ndim >= 7) p[7]=globals->trans_info.scales[0];
  if (ndim >7) {
    p[8]=globals->trans_info.scales[1];
    p[9]=globals->trans_info.scales[2];
  }
  
  if (ndim==12) {
    for_less( i, 0, 3 )		/* set shears */
      p[10+i] = globals->trans_info.shears[i];
  }
  
				/*  translation +/- simplex_size
				    rotation    +/- simplex_size*DEG_TO_RAD
				    scale       +/- simplex_size/50
				*/
  
  

  status = open_file(  globals->filenames.matlab_file, WRITE_FILE, BINARY_FORMAT,  &ofd );
  if ( status != OK ) 
    print_error ("filename `%s' cannot be opened.", 
		 __FILE__, __LINE__, globals->filenames.matlab_file);


  (void)fprintf (ofd,"%% %s\n",comments);

				/* do translations */
  for_inclusive(j,1,ndim) {
    switch (j) {
    case  1: (void)fprintf (ofd,"tx = [\n"); step=simplex_size/NUM_STEPS; break;
    case  2: (void)fprintf (ofd,"ty = [\n"); step=simplex_size/NUM_STEPS; break;
    case  3: (void)fprintf (ofd,"tz = [\n"); step=simplex_size/NUM_STEPS; break;
    case  4: (void)fprintf (ofd,"rx = [\n"); step=simplex_size*DEG_TO_RAD/NUM_STEPS; break;
    case  5: (void)fprintf (ofd,"ry = [\n"); step=simplex_size*DEG_TO_RAD/NUM_STEPS; break;
    case  6: (void)fprintf (ofd,"rz = [\n"); step=simplex_size*DEG_TO_RAD/NUM_STEPS; break;
    case  7: (void)fprintf (ofd,"sx = [\n"); step=simplex_size/(50*NUM_STEPS); break;
    case  8: (void)fprintf (ofd,"sy = [\n"); step=simplex_size/(50*NUM_STEPS); break;
    case  9: (void)fprintf (ofd,"sz = [\n"); step=simplex_size/(50*NUM_STEPS); break; 
    case 10: (void)fprintf (ofd,"hx = [\n"); step=simplex_size*DEG_TO_RAD/NUM_STEPS; break;
    case 11: (void)fprintf (ofd,"hy = [\n"); step=simplex_size*DEG_TO_RAD/NUM_STEPS; break;
    case 12: (void)fprintf (ofd,"hz = [\n"); step=simplex_size*DEG_TO_RAD/NUM_STEPS; break;
    }
    
    start = p[j];
    for_inclusive(i,-NUM_STEPS,NUM_STEPS) {
      p[j] = start + i*step;
      (void)fprintf (ofd, "%f %f\n",p[j], fit_function(p));
    }
    p[j] = start;
    
    (void)fprintf (ofd,"];\n"); 
    
  }
  
  status = close_file(ofd);
  if ( status != OK ) 
    print_error ("filename `%s' cannot be closed.", 
		 __FILE__, __LINE__, globals->filenames.matlab_file);
  
  
  free_vector(p,1,ndim);
  
}
