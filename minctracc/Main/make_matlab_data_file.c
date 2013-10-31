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
@MODIFIED   : Revision 96.10  2010-04-01 04:49:16  rotor
@MODIFIED   :  * fixed time.h include
@MODIFIED   :
@MODIFIED   : Revision 96.9  2008/10/08 15:17:49  louis
@MODIFIED   : added -nmi option for linear normalized mutual information
@MODIFIED   :
@MODIFIED   : Revision 96.8  2006/11/30 09:07:31  rotor
@MODIFIED   :  * many more changes for clean minc 2.0 build
@MODIFIED   :
@MODIFIED   : Revision 96.6  2005/07/20 20:45:48  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 96.5  2004/02/12 05:54:21  rotor
@MODIFIED   :  * removed /static defs
@MODIFIED   :
@MODIFIED   : Revision 96.4  2002/08/14 19:54:42  lenezet
@MODIFIED   :  quaternion option added for the rotation
@MODIFIED   :
@MODIFIED   : Revision 96.3  2002/03/26 14:15:38  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 96.2  2000/02/15 19:02:07  stever
@MODIFIED   : Add tests for param2xfm, minctracc -linear.
@MODIFIED   :
@MODIFIED   : Revision 96.1  1999/10/25 19:52:18  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :
 * Revision 96.0  1996/08/21  18:21:51  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:21:49  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:15:42  louis
 * Never released version 0.95
 *
 * Revision 1.6  1996/08/12  14:15:40  louis
 * Pre-release
 *
 * Revision 1.5  1995/02/22  08:56:06  collins
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.4  94/04/26  12:54:23  louis
 * updated with new versions of make_rots, extract2_parameters_from_matrix 
 * that include proper interpretation of skew.
 * 
 * Revision 1.3  94/04/06  11:48:39  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.2  94/02/21  16:35:40  louis
 * version before feb 22 changes
 * 
 * Revision 1.1  93/11/15  16:26:47  louis
 * Initial revision
 * 
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/minctracc/Main/make_matlab_data_file.c,v 96.10 2010-04-01 04:49:16 rotor Exp $";
#endif


#include <config.h>
#include <volume_io.h>

#include "constants.h"
#include "arg_data.h"
#include "objectives.h"
#include "segment_table.h"
#include "Proglib.h"

#include "local_macros.h"

extern Arg_Data main_args;

extern   VIO_Volume   Gdata1, Gdata2, Gmask1, Gmask2;
extern   int      Ginverse_mapping_flag, Gndim;
extern   double   simplex_size ;
extern   Segment_Table  *segment_table;

extern VIO_Real            **prob_hash_table; 
extern VIO_Real            *prob_fn1;         
extern VIO_Real            *prob_fn2;         

extern int Matlab_num_steps;

float fit_function(float *params);
float fit_function_quater(float *params);

void make_zscore_volume(VIO_Volume d1, VIO_Volume m1, 
                               VIO_Real *threshold); 

void add_speckle_to_volume(VIO_Volume d1, 
                                  float speckle,
                                  double  *start, int *count, VectorR directions[]);

void parameters_to_vector(double *trans, 
                                  double *rots, 
                                  double *scales,
                                  double *shears,
                                  float  *op_vector,
                                  double *weights);

void parameters_to_vector_quater(double *trans, 
                                        double *quats, 
                                        double *scales,
                                        double *shears,
                                        float  *op_vector,
                                        double *weights);

VIO_BOOL replace_volume_data_with_ubyte(VIO_Volume data);

void make_matlab_data_file(VIO_Volume d1,
                                  VIO_Volume d2,
                                  VIO_Volume m1,
                                  VIO_Volume m2, 
                                  char *comments,
                                  Arg_Data *globals)
{


  VIO_Status
    status;
  float 
    *p;
  int 
    i,j,stat, 
    ndim;
  FILE
    *ofd;
  VIO_Real
    start,step;
  double trans[3], quats[4], shears[3], scales[3],rots[3];
  VIO_Data_types
    data_type;

  start = 0.0;
  if (globals->obj_function == zscore_objective) { /* replace volume d1 and d2 by zscore volume  */
    make_zscore_volume(d1,m1,&globals->threshold[0]);
    make_zscore_volume(d2,m2,&globals->threshold[1]);
  } 
  else  if (globals->obj_function == ssc_objective) {        /* add speckle to the data set */
    
    make_zscore_volume(d1,m1,&globals->threshold[0]); /* need to make data sets comparable */
    make_zscore_volume(d2,m2,&globals->threshold[1]); /* in mean and sd...                 */
    
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
        print_error_and_line_num("%s",__FILE__, __LINE__,"Could not build segment table for source volume\n");
    }
    else {
      if (!build_segment_table(&segment_table, d2, globals->groups))
        print_error_and_line_num("%s",__FILE__, __LINE__,"Could not build segment table for target volume\n");
    }
    if (globals->flags.debug && globals->flags.verbose>1) {
      print ("groups = %d\n",segment_table->groups);
      for(i=segment_table->min; i<segment_table->max+1; i++) {
        print ("%5d: table = %5d, function = %5d\n",i,segment_table->table[i],
               (segment_table->segment)(i,segment_table) );
      }
    }
    
  } else if (globals->obj_function == mutual_information_objective || globals->obj_function == normalized_mutual_information_objective)
                                /* Collignon's mutual information */
    {

      if ( globals->groups != 256 ) {
        print ("WARNING: -groups was %d, but will be forced to 256 in this run\n",globals->groups);
        globals->groups = 256;
      }

      data_type = get_volume_data_type (d1);
      if (data_type != VIO_UNSIGNED_BYTE) {
        print ("WARNING: source volume not UNSIGNED BYTE, will do conversion now.\n");
        if (!replace_volume_data_with_ubyte(d1)) {
          print_error_and_line_num("Can't replace volume data with unsigned bytes\n",
                             __FILE__, __LINE__);
        }
      }

      data_type = get_volume_data_type (d2);
      if (data_type != VIO_UNSIGNED_BYTE) {
        print ("WARNING: target volume not UNSIGNED BYTE, will do conversion now.\n");
        if (!replace_volume_data_with_ubyte(d2)) {
          print_error_and_line_num("Can't replace volume data with unsigned bytes\n",
                             __FILE__, __LINE__);
        }
      }

      ALLOC(   prob_fn1,   globals->groups);
      ALLOC(   prob_fn2,   globals->groups);
      VIO_ALLOC2D( prob_hash_table, globals->groups, globals->groups);

    } 


if(globals->trans_info.rotation_type == TRANS_ROT)
  {
  /* ---------------- prepare the weighting array for for the objective function  ---------*/

  stat = TRUE;
  switch (globals->trans_info.transform_type) {
  case TRANS_LSQ3: 
    for(i=3; i<12; i++) globals->trans_info.weights[i] = 0.0;
    for(i=0; i<3; i++) {
      globals->trans_info.scales[i] = 1.0;
      globals->trans_info.rotations[i] = 0.0;
      globals->trans_info.shears[i] = 0.0;
    }
    break;
  case TRANS_LSQ6: 
    for(i=6; i<12; i++) globals->trans_info.weights[i] = 0.0;
    for(i=0; i<3; i++) {
      globals->trans_info.scales[i] = 1.0;
      globals->trans_info.shears[i] = 0.0;
    }
    break;
  case TRANS_LSQ7: 
    for(i=7; i<12; i++) globals->trans_info.weights[i] = 0.0;
    for(i=0; i<3; i++) {
      globals->trans_info.shears[i] = 0.0;
    }
    break;
  case TRANS_LSQ9: 
    for(i=9; i<12; i++) globals->trans_info.weights[i] = 0.0;
    for(i=0; i<3; i++) {
      globals->trans_info.shears[i] = 0.0;
    }
    break;
  case TRANS_LSQ10: 
    for(i=10; i<12; i++) globals->trans_info.weights[i] = 0.0;
    for(i=1; i<3; i++) {
      globals->trans_info.shears[i] = 0.0;
    }
    break;
  case TRANS_LSQ: 
                                /* nothing to be zeroed */
    break;
  case TRANS_LSQ12: 
                                /* nothing to be zeroed */
    break;
  default:
    (void)fprintf(stderr, "Unknown type of transformation requested (%d)\n",
                   globals->trans_info.transform_type);
    (void)fprintf(stderr, "Error in line %d, file %s\n",__LINE__, __FILE__);
    stat = FALSE;
  }

  if ( !stat ) 
    print_error_and_line_num ("Can't calculate measure (stat is false).", 
                              __FILE__, __LINE__);


                                /* find number of dimensions for obj function */
  ndim = 0;
  for(i=0; i<12; i++)
    if (globals->trans_info.weights[i] != 0.0) ndim++;


                                /* set GLOBALS to communicate with the
                                   function to be fitted!              */
  Gndim = ndim;
  Gdata1 = d1;
  Gdata2 = d2;
  Gmask1 = m1;
  Gmask2 = m2;
  Ginverse_mapping_flag = FALSE;

print ("trans: %10.5f %10.5f %10.5f \n",
       globals->trans_info.translations[0],globals->trans_info.translations[1],globals->trans_info.translations[2]);
print ("rots : %10.5f %10.5f %10.5f \n",
       globals->trans_info.rotations[0],globals->trans_info.rotations[1],globals->trans_info.rotations[2]);
print ("scale: %10.5f %10.5f %10.5f \n",
       globals->trans_info.scales[0],globals->trans_info.scales[1],globals->trans_info.scales[2]);


  if (ndim>0) {

    ALLOC(p,ndim+1); /* parameter values */
    
    /*  translation +/- simplex_size
        rotation    +/- simplex_size*DEG_TO_RAD
        scale       +/- simplex_size/50
        */
    
    
    
    status = open_file(  globals->filenames.matlab_file, WRITE_FILE, BINARY_FORMAT,  &ofd );
    if ( status != VIO_OK ) 
      print_error_and_line_num ("filename `%s' cannot be opened.", 
                   __FILE__, __LINE__, globals->filenames.matlab_file);
    
    
    (void)fprintf (ofd,"%% %s\n",comments);
    
    /* do translations */
    for(j=1; j<=12; j++) {
      if (globals->trans_info.weights[j-1] != 0.0) {
        switch (j) {
        case  1: (void)fprintf (ofd,"tx = [\n"); start = globals->trans_info.translations[0]; break;
        case  2: (void)fprintf (ofd,"ty = [\n"); start = globals->trans_info.translations[1]; break;
        case  3: (void)fprintf (ofd,"tz = [\n"); start = globals->trans_info.translations[2]; break;
        case  4: (void)fprintf (ofd,"rx = [\n"); start = globals->trans_info.rotations[0]; break;
        case  5: (void)fprintf (ofd,"ry = [\n"); start = globals->trans_info.rotations[1]; break;
        case  6: (void)fprintf (ofd,"rz = [\n"); start = globals->trans_info.rotations[2]; break;
        case  7: (void)fprintf (ofd,"sx = [\n"); start = globals->trans_info.scales[0]; break;
        case  8: (void)fprintf (ofd,"sy = [\n"); start = globals->trans_info.scales[1]; break;
        case  9: (void)fprintf (ofd,"sz = [\n"); start = globals->trans_info.scales[2]; break;
        case 10: (void)fprintf (ofd,"shx = [\n");start = globals->trans_info.shears[0]; break;
        case 11: (void)fprintf (ofd,"shy = [\n");start = globals->trans_info.shears[1]; break;
        case 12: (void)fprintf (ofd,"shz = [\n");start = globals->trans_info.shears[2]; break; 
        }
        
        for(i=0; i<3; i++) {
          trans[i]  = globals->trans_info.translations[i];
          scales[i] = globals->trans_info.scales[i];
          shears[i] = globals->trans_info.shears[i];
          rots[i]   = globals->trans_info.rotations[i];
        }

        step =  globals->trans_info.weights[j-1] * simplex_size/ Matlab_num_steps;
        
        for(i=-Matlab_num_steps; i<=Matlab_num_steps; i++) {
          
          switch (j) {
          case  1: trans[0] =start + i*step; break;
          case  2: trans[1] =start + i*step; break;
          case  3: trans[2] =start + i*step; break;
          case  4: rots[0]  =start + i*step; break;
          case  5: rots[1]  =start + i*step; break;
          case  6: rots[2]  =start + i*step; break;
          case  7: scales[0]=start + i*step; break;
          case  8: scales[1]=start + i*step; break;
          case  9: scales[2]=start + i*step; break;
          case 10: shears[0]=start + i*step; break;
          case 11: shears[1]=start + i*step; break;
          case 12: shears[2]=start + i*step; break; 
          }
          
          parameters_to_vector(trans,
                               rots,
                               scales,
                               shears,
                               p,
                               globals->trans_info.weights);
    
          (void)fprintf (ofd, "%f %f %f\n",i*step, start+i*step, fit_function(p));
        }

        (void)fprintf (ofd,"];\n"); 
      }
    }
    
    status = close_file(ofd);
    if ( status != VIO_OK ) 
      print_error_and_line_num ("filename `%s' cannot be closed.", 
                   __FILE__, __LINE__, globals->filenames.matlab_file);
    
    
    FREE(p);
  }
  }


    
if(globals->trans_info.rotation_type == TRANS_QUAT)
  {
  /* ---------------- prepare the weighting array for for the objective function  ---------*/

  stat = TRUE;
  switch (globals->trans_info.transform_type) {
  case TRANS_LSQ3: 
    for(i=3; i<13; i++) globals->trans_info.weights[i] = 0.0;
    for(i=0; i<3; i++) {
      globals->trans_info.scales[i] = 1.0;
      globals->trans_info.quaternions[i] = 0.0;
      globals->trans_info.shears[i] = 0.0;
    }
    globals->trans_info.quaternions[3] = 0.1;
    break;
  case TRANS_LSQ6: 
    for(i=7; i<13; i++) globals->trans_info.weights[i] = 0.0;
    for(i=0; i<3; i++) {
      globals->trans_info.scales[i] = 1.0;
      globals->trans_info.shears[i] = 0.0;
    }
    break;
  case TRANS_LSQ7: 
    for(i=8; i<13; i++) globals->trans_info.weights[i] = 0.0;
    for(i=0; i<3; i++) {
      globals->trans_info.shears[i] = 0.0;
    }
    break;
  case TRANS_LSQ9: 
    for(i=10; i<13; i++) globals->trans_info.weights[i] = 0.0;
    for(i=0; i<3; i++) {
      globals->trans_info.shears[i] = 0.0;
    }
    break;
  case TRANS_LSQ10: 
    for(i=11; i<13; i++) globals->trans_info.weights[i] = 0.0;
    for(i=1; i<3; i++) {
      globals->trans_info.shears[i] = 0.0;
    }
    break;
  case TRANS_LSQ: 
                                /* nothing to be zeroed */
    break;
  case TRANS_LSQ12: 
                                /* nothing to be zeroed */
    break;
  default:
    (void)fprintf(stderr, "Unknown type of transformation requested (%d)\n",
                   globals->trans_info.transform_type);
    (void)fprintf(stderr, "Error in line %d, file %s\n",__LINE__, __FILE__);
    stat = FALSE;
  }

  if ( !stat ) 
    print_error_and_line_num ("Can't calculate measure (stat is false).", 
                              __FILE__, __LINE__);


                                /* find number of dimensions for obj function */
  ndim = 0;
  for(i=0; i<13; i++)
    if (globals->trans_info.weights[i] != 0.0) ndim++;


                                /* set GLOBALS to communicate with the
                                   function to be fitted!              */
  Gndim = ndim;
  Gdata1 = d1;
  Gdata2 = d2;
  Gmask1 = m1;
  Gmask2 = m2;
  Ginverse_mapping_flag = FALSE;

print ("trans: %10.5f %10.5f %10.5f \n",
       globals->trans_info.translations[0],globals->trans_info.translations[1],globals->trans_info.translations[2]);
print ("quats : %10.5f %10.5f %10.5f  %10.5f\n",
       globals->trans_info.quaternions[0],globals->trans_info.quaternions[1],globals->trans_info.quaternions[2],globals->trans_info.quaternions[3]);
print ("scale: %10.5f %10.5f %10.5f \n",
       globals->trans_info.scales[0],globals->trans_info.scales[1],globals->trans_info.scales[2]);


  if (ndim>0) {

    ALLOC(p,ndim+1); /* parameter values */
    
    /*  translation +/- simplex_size
        rotation    +/- simplex_size*DEG_TO_RAD
        scale       +/- simplex_size/50
        */
    
    
    
    status = open_file(  globals->filenames.matlab_file, WRITE_FILE, BINARY_FORMAT,  &ofd );
    if ( status != VIO_OK ) 
      print_error_and_line_num ("filename `%s' cannot be opened.", 
                   __FILE__, __LINE__, globals->filenames.matlab_file);
    
    
    (void)fprintf (ofd,"%% %s\n",comments);
    
    /* do translations */
    for(j=1; j<=13; j++) {
      if (globals->trans_info.weights[j-1] != 0.0) {
        switch (j) {
        case  1: (void)fprintf (ofd,"tx = [\n"); start = globals->trans_info.translations[0]; break;
        case  2: (void)fprintf (ofd,"ty = [\n"); start = globals->trans_info.translations[1]; break;
        case  3: (void)fprintf (ofd,"tz = [\n"); start = globals->trans_info.translations[2]; break;
        case  4: (void)fprintf (ofd,"rx = [\n"); start = globals->trans_info.quaternions[0]; break;
        case  5: (void)fprintf (ofd,"ry = [\n"); start = globals->trans_info.quaternions[1]; break;
        case  6: (void)fprintf (ofd,"rz = [\n"); start = globals->trans_info.quaternions[2]; break;
        case  7: (void)fprintf (ofd,"rz = [\n"); start = globals->trans_info.quaternions[3]; break;
        case  8: (void)fprintf (ofd,"sx = [\n"); start = globals->trans_info.scales[0]; break;
        case  9: (void)fprintf (ofd,"sy = [\n"); start = globals->trans_info.scales[1]; break;
        case 10: (void)fprintf (ofd,"sz = [\n"); start = globals->trans_info.scales[2]; break;
        case 11: (void)fprintf (ofd,"shx = [\n");start = globals->trans_info.shears[0]; break;
        case 12: (void)fprintf (ofd,"shy = [\n");start = globals->trans_info.shears[1]; break;
        case 13: (void)fprintf (ofd,"shz = [\n");start = globals->trans_info.shears[2]; break; 
        }
        
        for(i=0; i<3; i++) {
          trans[i]  = globals->trans_info.translations[i];
          scales[i] = globals->trans_info.scales[i];
          shears[i] = globals->trans_info.shears[i];
          quats[i]   = globals->trans_info.quaternions[i];
        }
        quats[3] = globals->trans_info.quaternions[3];

        step =  globals->trans_info.weights[j-1] * simplex_size/ Matlab_num_steps;
        
        for(i=-Matlab_num_steps; i<=Matlab_num_steps; i++) {
          
          switch (j) {
          case  1: trans[0] =start + i*step; break;
          case  2: trans[1] =start + i*step; break;
          case  3: trans[2] =start + i*step; break;
          case  4: quats[0]  =start + i*step; break;
          case  5: quats[1]  =start + i*step; break;
          case  6: quats[2]  =start + i*step; break;
          case  7: quats[3]  =start + i*step; break;
          case  8: scales[0]=start + i*step; break;
          case  9: scales[1]=start + i*step; break;
          case 10: scales[2]=start + i*step; break;
          case 11: shears[0]=start + i*step; break;
          case 12: shears[1]=start + i*step; break;
          case 13: shears[2]=start + i*step; break; 
          }
          
          parameters_to_vector_quater(trans,
                                      quats,
                                      scales,
                                      shears,
                                      p,
                                      globals->trans_info.weights);
    
          (void)fprintf (ofd, "%f %f %f\n",i*step, start+i*step, fit_function_quater(p));
        }

        (void)fprintf (ofd,"];\n"); 
      }
    }
    
    status = close_file(ofd);
    if ( status != VIO_OK ) 
      print_error_and_line_num ("filename `%s' cannot be closed.", 
                   __FILE__, __LINE__, globals->filenames.matlab_file);
    
    
    FREE(p);
  }
  }
  if (globals->obj_function == vr_objective) {
    if (!free_segment_table(segment_table)) {
      (void)fprintf(stderr, "Can't free segment table.\n");
      (void)fprintf(stderr, "Error in line %d, file %s\n",__LINE__, __FILE__);
    }
  } else
  if (globals->obj_function == mutual_information_objective || globals->obj_function == normalized_mutual_information_objective )
                                /* Collignon's mutual information */
    {
      FREE(   prob_fn1 );
      FREE(   prob_fn2 );
      VIO_FREE2D( prob_hash_table);
    }



}
