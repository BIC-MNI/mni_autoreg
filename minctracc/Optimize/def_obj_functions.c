/*------------------------------ MNI Header ----------------------------------
@NAME       : def_obj_functions.c
@DESCRIPTION: routines to calculate the objective function used for local
              optimization              
@CREATED    : Nov 4, 1997, Louis Collins
@MODIFIED   : not yet!
@VERSION    : $Id: def_obj_functions.c,v 1.2 1998-03-13 20:16:44 louis Exp $
-----------------------------------------------------------------------------*/

#include <config.h>		
#include <internal_volume_io.h>	
#include "constants.h"
#include <arg_data.h>           /* definition of the global data struct      */


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
      val[MAX_DIMENSIONS],
      voxel[MAX_DIMENSIONS],
      xw,yw,zw, 
      s, s1, s2, func_sim;
   
   /* note: here the displacement order for go_get_samples_with_offset
      is 3,2,1 (=Z,Y,X) since the source and target volumes are stored in
      Z,Y,X order.
      
      This backward ordering is necessary since the same function is
      used for 2D and 3D optimization, and simplex will only modify d[1]
      and d[2] in 2D (d[3] stays const=0). */
   
   if (Gglobals->trans_info.use_magnitude) {
      
      s = norm = 0.0;
      
      for_less(i,0,Gglobals->features.number_of_features)  {
         func_sim = 
            (Real)go_get_samples_with_offset(Gglobals->features.model[i],
                                             TX,TY,TZ,
                                             d[3], d[2], d[1],
                                             Gglobals->features.obj_func[i],
                                             Glen, 
                                             Gsqrt_features[i], Ga1_features[i],
                                             Gglobals->interpolant==nearest_neighbour_interpolant);
         if ((Gglobals->features.obj_func[i]==NONLIN_CHAMFER) && (func_sim > 1.5)) {
                                /* do nothing, do not add the chamfer distance info */
         }
         else {
            norm += ABS(Gglobals->features.weight[i]);
            s += Gglobals->features.weight[i] * func_sim;
         }
/*  print ("s%12.9f ",func_sim); !!! */
         
      }
      
      if (norm != 0.0) 
         s = s / norm;
      else
         print_error_and_line_num("The feature weights are null.", 
                                  __FILE__, __LINE__);
      
      
   } else {
      /* calc correlation based on projection data */
      voxel[0] = Gtarget_vox_x + d[3];
      voxel[1] = Gtarget_vox_y + d[2]; 
      voxel[2] = Gtarget_vox_z + d[1];
      voxel[3] = 0.0;
      voxel[4] = 0.0;
      
      convert_voxel_to_world(Gglobals->features.model[0], voxel, &xw, &yw, &zw);
      
      evaluate_volume_in_world(Gglobals->features.model[0], xw, yw, zw, 
                               0, TRUE, 0.0, val,
                               NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
      Gproj_d2 = val[0];
      evaluate_volume_in_world(Gglobals->features.model[1], xw, yw, zw, 
                               0, TRUE, 0.0, val,
                               NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
      Gproj_d2x = 2.0*val[0];
      evaluate_volume_in_world(Gglobals->features.model[2], xw, yw, zw, 
                               0, TRUE, 0.0, val,
                               NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
      Gproj_d2y = 2.0*val[0];
      evaluate_volume_in_world(Gglobals->features.model[3], xw, yw, zw, 
                               0, TRUE, 0.0, val,
                               NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
      Gproj_d2z = 2.0*val[0];
      
      s  = (Gproj_d1*Gproj_d2   + Gproj_d1x*Gproj_d2x + 
            Gproj_d1y*Gproj_d2y + Gproj_d1z*Gproj_d2z);
      s1 = (Gproj_d1*Gproj_d1   +  Gproj_d1x*Gproj_d1x + 
            Gproj_d1y*Gproj_d1y + Gproj_d1z*Gproj_d1z);
      s2 = (Gproj_d2*Gproj_d2   + Gproj_d2x*Gproj_d2x + 
            Gproj_d2y*Gproj_d2y + Gproj_d2z*Gproj_d2z);
      
      if ((s1!=0.0) && (s2 !=0.0)) {
         s = s / ( sqrt(s1) * sqrt(s2) );
      }
      else {
         s = 0.0;
      }
      
   }
   return( s );
}

/* this is the function that needs to be minimized */
public Real local_objective_function(float *d)
     
{
  Real
    similarity,
    cost, 
    r;
  
  similarity = (Real)similarity_fn( d );
  cost       = (Real)cost_fn( d[1], d[2], d[3], Gcost_radius );

/*  print ("c%12.9f ",cost);  !!! */
  
  r = 1.0 - 
      similarity * similarity_cost_ratio + 
      cost       * (1.0-similarity_cost_ratio);

/*  print ("  -> %12.9f\n",r);  !!! */
  
  return(r);
}


/*  amoeba_NL_obj_function() is minimized in the amoeba() optimization function
*/
public Real amoeba_NL_obj_function(void * dummy, float d[])
{
  int i;
  float p[4];
  Real obj_func_val;


  for_less(i,0,number_dimensions)
    p[i+1] = d[i];
  obj_func_val =  local_objective_function(p);

  return ( obj_func_val );
  
}

