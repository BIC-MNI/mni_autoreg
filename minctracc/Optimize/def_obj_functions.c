/*------------------------------ MNI Header ----------------------------------
@NAME       : def_obj_functions.c
@DESCRIPTION: routines to calculate the objective function used for local
              optimization              
@CREATED    : Nov 4, 1997, Louis Collins
@VERSION    : $Id: def_obj_functions.c,v 1.8 2003-02-26 00:56:37 lenezet Exp $
@MODIFIED   : $Log: def_obj_functions.c,v $
@MODIFIED   : Revision 1.8  2003-02-26 00:56:37  lenezet
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
#include <volume_io/internal_volume_io.h>	
#include "constants.h"
#include <arg_data.h>           /* definition of the global data struct      */
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
extern int 
  Glen;                         /* from do_nonlinear.c */
extern Real                     /* from do_nonlinear.c */
  Gtarget_vox_x, Gtarget_vox_y, Gtarget_vox_z,
  Gproj_d1,  Gproj_d1x,  Gproj_d1y,  Gproj_d1z, 
  Gproj_d2,  Gproj_d2x,  Gproj_d2y,  Gproj_d2z;
extern double
  similarity_cost_ratio;
extern Real     
  Gcost_radius;                 /* from do_nonlinear.c */
public int 
  nearest_neighbour_interpolant(Volume volume, 
                                PointR *coord, double *result);

public void from_param_to_grid_weights(
   Real p[],
   Real grid[]);


public float go_get_samples_with_offset(Volume data,
					float *x, float *y, float *z,
					Real  dx, Real  dy, Real dz,
					int obj_func,
					int len, 
					float sqrt_s1, float *a1,
					BOOLEAN use_nearest_neighbour);


/* This is the COST FUNCTION TO BE MINIMIZED.
   so that very large displacements are impossible */

private Real cost_fn(float x, float y, float z, Real max_length)
{
  Real v2,v,d;

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

private Real similarity_fn(float *d)
{
  int i;
  Real
    norm,
    s, func_sim;
   
  /* note: here the displacement order for go_get_samples_with_offset
     is 3,2,1 (=Z,Y,X) since the source and target volumes are stored in
     Z,Y,X order.
     
     This backward ordering is necessary since the same function is
     used for 2D and 3D optimization, and simplex will only modify d[1]
     and d[2] in 2D (d[3] stays const=0). */
  
  s = norm = 0.0;
    
  for_less(i,0,Gglobals->features.number_of_features)  {

				/* ignore OPTICAL FLOW objective functions, since it is
				   computed directly and _not_ optimized */

    if (Gglobals->features.obj_func[i] != NONLIN_OPTICALFLOW) {
      func_sim = 
	(Real)go_get_samples_with_offset(Gglobals->features.model[i],
					 TX,TY,TZ,
					 d[3], d[2], d[1],
					 Gglobals->features.obj_func[i],
					 Glen, 
					 Gsqrt_features[i], Ga1_features[i],
					 Gglobals->interpolant==nearest_neighbour_interpolant);
      

      norm += ABS(Gglobals->features.weight[i]);
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
public Real local_objective_function(float *d)
     
{
  Real
    similarity,
    cost, 
    r;
  
  similarity = (Real)similarity_fn( d );
  cost       = (Real)cost_fn( d[1], d[2], d[3], Gcost_radius );
  
  r = 1.0 - 
      similarity * similarity_cost_ratio + 
      cost       * (1.0-similarity_cost_ratio);

  return(r);
}


/*  
    amoeba_NL_obj_function() is minimized in the amoeba() optimization function
*/
public Real amoeba_NL_obj_function(void * dummy, float d[])
{
  int i;
  float p[4];
  Real
    real_d[N_DIMENSIONS],
    grid_weights[N_DIMENSIONS],
    obj_func_val;

  /* for_less(i,0,number_dimensions)
      p[i+1] = d[i]; */
  
  for_less(i,0,number_dimensions)
    real_d[i] = d[i];


  from_param_to_grid_weights( real_d, grid_weights);


  /*  for_less(i,0,number_dimensions)*/

  for_less(i,0,N_DIMENSIONS)
    p[i+1] = (float)grid_weights[i];


  obj_func_val =  local_objective_function(p);


  return ( obj_func_val );
}

