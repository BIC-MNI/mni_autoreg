/*------------------------------ MNI Header ----------------------------------
@NAME       : def_obj_functions.c
@DESCRIPTION: routines to calculate the objective function used for local
              optimization              
@CREATED    : Nov 4, 1997, Louis Collins
@VERSION    : $Id: def_obj_functions.c,v 1.14 2006-11-30 09:07:32 rotor Exp $
@MODIFIED   : $Log: def_obj_functions.c,v $
@MODIFIED   : Revision 1.14  2006-11-30 09:07:32  rotor
@MODIFIED   :  * many more changes for clean minc 2.0 build
@MODIFIED   :
@MODIFIED   : Revision 1.13  2006/11/29 09:09:34  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 1.12  2005/07/20 20:45:50  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 1.11  2005/06/28 18:56:18  rotor
@MODIFIED   :  * added masking for feature volumes (irina and patricia)
@MODIFIED   :
@MODIFIED   : Revision 1.10  2004/02/12 06:08:19  rotor
@MODIFIED   :  * removed /static defs
@MODIFIED   :
@MODIFIED   : Revision 1.9  2004/02/04 20:44:11  lenezet
@MODIFIED   : *** empty log message ***
@MODIFIED   :
@MODIFIED   : Revision 1.8  2003/02/26 00:56:37  lenezet
@MODIFIED   : for 2D : now computes all 3 coordinates for the "start" (to take into account the slice position).
@MODIFIED   : simplification of build_lattices.
@MODIFIED   : bug correction in amoeba_NL_obj_function.
@MODIFIED   :
@MODIFIED   : Revision 1.7  2002/12/13 21:18:19  lenezet
@MODIFIED   :
@MODIFIED   : A memory leak has been repaired
@MODIFIED   :
@MODIFIED   : Revision 1.6  2002/03/26 14:15:43  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 1.5  2000/03/15 08:42:45  stever
@MODIFIED   : Code cleanup: all functions prototyped (except ParseArgs.c), no useless declarations, etc
@MODIFIED   :
@MODIFIED   : Revision 1.4  1999/06/10 12:51:23  louis
@MODIFIED   : update with optical flow working in addition to xcorr, label, and diff
@MODIFIED   : sub-lattice computed only if needed
@MODIFIED   :
 * Revision 1.3  1999/06/09  13:10:51  louis
 * clean up
 *
-----------------------------------------------------------------------------*/

#include <config.h>                
#include <volume_io.h>        
#include "constants.h"
#include <minctracc_arg_data.h>           /* definition of the global data struct      */
#include <Proglib.h>


/* GLOBALS used within these functions: */

extern int 
  number_dimensions;            /* from do_nonlinear.c */
extern Arg_Data 
  *Gglobals;                    /* from do_nonlinear.c */
extern float
  *Gsqrt_features,
  **Ga1_features,
  *TX, *TY, *TZ;                /* from do_nonlinear.c */
extern VIO_BOOL 
  **masked_samples_in_source;  /* from do_nonlinear.c */
extern int 
  Glen;                         /* from do_nonlinear.c */
extern VIO_Real                     /* from do_nonlinear.c */
  Gtarget_vox_x, Gtarget_vox_y, Gtarget_vox_z,
  Gproj_d1,  Gproj_d1x,  Gproj_d1y,  Gproj_d1z, 
  Gproj_d2,  Gproj_d2x,  Gproj_d2y,  Gproj_d2z;
extern double
  similarity_cost_ratio;
extern VIO_Real     
  Gcost_radius;                 /* from do_nonlinear.c */
int 
  nearest_neighbour_interpolant(VIO_Volume volume, 
                                PointR *coord, double *result);
int target_sample_count=0;

void from_param_to_grid_weights(
   VIO_Real p[],
   VIO_Real grid[]);


float go_get_samples_with_offset(VIO_Volume data, VIO_Volume model_mask,
                                        float *x, float *y, float *z,
                                        VIO_Real  dx, VIO_Real  dy, VIO_Real dz,
                                        int obj_func,
                                        int len, int *sample_count,
                                        float sqrt_s1, float *a1, VIO_BOOL *m1,
                                        VIO_BOOL use_nearest_neighbour);


/* This is the COST FUNCTION TO BE MINIMIZED.
   so that very large displacements are impossible */

static VIO_Real cost_fn(float x, float y, float z, VIO_Real max_length)
{
  VIO_Real v2,v,d;

  v2 = x*x + y*y + z*z;
  v = sqrt(v2);

  v *= v2;

  if (v<max_length)
    d = 0.2 * v / (max_length - v);
  else
    d = 1e38;

  return(d);
}

/* This is the SIMILARITY FUNCTION TO BE MAXIMIZED.

      it is maximum when source and target data are most similar.
      The value is normalized by the weights associated with each feature,
      so that -1.0 <= similarity_fn() <= 0.0
      where similarity_fn() = 0.0 is BEST

      note 0.0 < go_get_samples_with_offset() < 1.0; with 0.0 BEST

      the input parameter D stores the displacement offsets for the
      target lattice:
         D[1] stores the xdisp in voxels
         D[2] stores the ydisp
         D[3] stores the zdisp
*/

static VIO_Real similarity_fn(float *d)
{
  int i;
  VIO_Real
    norm,
    s, func_sim;
   
  /* note: here the displacement order for go_get_samples_with_offset
     is 3,2,1 (=Z,Y,X) since the source and target volumes are stored in
     Z,Y,X order.
     
     This backward ordering is necessary since the same function is
     used for 2D and 3D optimization, and simplex will only modify d[1]
     and d[2] in 2D (d[3] stays const=0). */
  
  s = norm = 0.0;
    
  for(i=0; i<Gglobals->features.number_of_features; i++)  {

                                /* ignore OPTICAL FLOW objective functions, since it is
                                   computed directly and _not_ optimized */

    if (Gglobals->features.obj_func[i] != NONLIN_OPTICALFLOW) {
      func_sim = 
        (VIO_Real)go_get_samples_with_offset(Gglobals->features.model[i],
                Gglobals->features.model_mask[i],
                                         TX,TY,TZ,
                                         d[3], d[2], d[1],
                                         Gglobals->features.obj_func[i],
                                         Glen, &target_sample_count,
                                         Gsqrt_features[i], Ga1_features[i],
                masked_samples_in_source[i],
                                         Gglobals->interpolant==nearest_neighbour_interpolant);
      

      norm += fabs(Gglobals->features.weight[i]);
      s += Gglobals->features.weight[i] * func_sim;
      
      /*
        if ((Gglobals->features.obj_func[i]==NONLIN_CHAMFER) && (func_sim > 1.5))
        do nothing, do not add the chamfer distance info 
      */
      
    }
    
  }
  
  if (norm > 0.0) 
    s = s / norm;
  else
    print_error_and_line_num("The feature weights are null.", 
                             __FILE__, __LINE__);

  return( s );
}

/* 
   this is the objective function that needs to be minimized 
   to give a local deformation
*/
VIO_Real local_objective_function(float *d)
     
{
  VIO_Real
    similarity,
    cost, 
    r;
  
  similarity = (VIO_Real)similarity_fn( d );
  cost       = (VIO_Real)cost_fn( d[1], d[2], d[3], Gcost_radius );
  
  r = 1.0 - 
      similarity * similarity_cost_ratio + 
      cost       * (1.0-similarity_cost_ratio);

  return(r);
}


/*  
    amoeba_NL_obj_function() is minimized in the amoeba() optimization function
*/
VIO_Real amoeba_NL_obj_function(void * dummy, float d[])
{
  int i;
  float p[4];
  VIO_Real
    real_d[VIO_N_DIMENSIONS],
    grid_weights[VIO_N_DIMENSIONS],
    obj_func_val;
  
  for(i=0; i<number_dimensions; i++)
    real_d[i] = d[i];


  from_param_to_grid_weights( real_d, grid_weights);


  for(i=0; i<VIO_N_DIMENSIONS; i++)
    p[i+1] = (float)grid_weights[i];


  obj_func_val =  local_objective_function(p);


  return ( obj_func_val );
}

