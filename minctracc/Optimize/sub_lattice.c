/*
------------------------------ MNI Header ----------------------------------
@NAME       : sub_lattice.c
@DESCRIPTION: this file contains a number of rountines used to manipulate the 
              local sublattice for non-linear deformation.

  these include:
     build_source_lattice() -       to build the sublattice in the source volume
     go_get_samples_in_source() -   to interpolate values for sublattice positions 
                                    in the source volume
     go_get_samples_with_offset() - to interpolate values for sublattice positions, 
                                    given a vector offset for the lattice.
     build_target_lattice() -       map the source sublattice thrugh the current xform
                                    to create a sublattice defined on the target.
     
@CREATED    : Mon Nov  3, 1997 , Louis Collins
@VERSION    : $Id: sub_lattice.c,v 1.5 2002-03-26 14:15:46 stever Exp $
@MODIFIED   : $Log: sub_lattice.c,v $
@MODIFIED   : Revision 1.5  2002-03-26 14:15:46  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 1.4  2000/03/15 08:42:47  stever
@MODIFIED   : Code cleanup: all functions prototyped (except ParseArgs.c), no useless declarations, etc
@MODIFIED   :
@MODIFIED   : Revision 1.3  1998/03/13 20:22:22  louis
@MODIFIED   : Modified the loops over the nodes of the sublattice so that no interpolation
@MODIFIED   : in the target volume is done in go_get_samples_with_offset() when the
@MODIFIED   : objective function == NONLIN_CHAMFER and there is no sulci (*a1==0)
@MODIFIED   : in the node list for the source sub-lattice.  This yields a considerable
@MODIFIED   : improvement in speed (62% faster) for the evaluation of chamfer feature.
@MODIFIED   : Overall speedup is on the order of 25-30% when using 2 features, with
@MODIFIED   : no measurable cost when using only a single feature.
@MODIFIED   :
#-----------------------------------------------------------------------------
*/


#include <config.h>		
#include <volume_io/internal_volume_io.h>	

#include <Proglib.h>
#include "constants.h"
#include "arg_data.h"		/* definition of the global data struct      */
#include "sub_lattice.h"
#include "init_lattice.h"


extern Arg_Data *Gglobals;      /* defined in do_nonlinear.c */
extern Volume   Gsuper_sampled_vol; /* defined in do_nonlinear.c */
extern General_transform 
                *Glinear_transform;/* defined in do_nonlinear.c */

                                /* prototypes for functions used here: */

extern float
  *SX, *SY, *SZ;

public  void  general_transform_point_in_trans_plane(
    General_transform   *transform,
    Real                x,
    Real                y,
    Real                z,
    Real                *x_transformed,
    Real                *y_transformed,
    Real                *z_transformed );


/*********************************************************************** 
   build a regular (2D) 3D lattice of coordinates to represent the
   local (circular) spherical neighbourhood surrounding the point
   x,y,z.

   - *length coordinates are returned in PX[], PY[], PZ[]. 
   - The radius of the lattice is defined by width_{x,y,z}.
   - The equivalent retangular lattice has (nx)(ny)(nz) samples,
     but the round lattice returned has *length samples.  

   the coordinate coming in, and those going out are all in WORLD
   coordinates 
*/

public void    build_source_lattice(Real x, Real y, Real z,
				    float PX[], float PY[], float PZ[],
				    Real width_x, Real width_y, Real width_z, 
				    int nx, int ny, int nz,
				    int ndim, int *length)
{
  int 
    c, 
    i,j,k;
  float 
    radius_squared,
    tx,ty,tz;

  *length = 0; 
  c = 1;
  radius_squared = 0.55 * 0.55;	/* a bit bigger than .5^2 */
  
  if (ndim==2) {
    for_less(i,0,nx)
      for_less(j,0,ny) {
	tx = -0.5 + (float)(i)/(float)(nx-1);
	ty = -0.5 + (float)(j)/(float)(ny-1);
	if ((tx*tx + ty*ty) <= radius_squared) {
	  PX[c] = (float)(x + width_x * tx);
	  PY[c] = (float)(y + width_y * ty);
	  PZ[c] = (float)z;
	  c++;
	  (*length)++;
	}
      }
  }
  else {
    for_less(i,0,nx)
      for_less(j,0,ny)
	for_less(k,0,nz) {
	  tx = -0.5 + (float)(i)/(float)(nx-1);
	  ty = -0.5 + (float)(j)/(float)(ny-1);
	  tz = -0.5 + (float)(k)/(float)(nz-1);
	  if ((tx*tx + ty*ty + tz*tz) <= radius_squared) {
	    PX[c] = (float)(x + width_x * tx);
	    PY[c] = (float)(y + width_y * ty);
	    PZ[c] = (float)(z + width_z * tz);
	    c++;
	    (*length)++;
	  }
	}
  }
}

/*********************************************************************** 
   use the world coordinates stored in x[],y[],z[] to interpolate len
   samples from the volume 'data' */

public void go_get_samples_in_source(Volume data,
				     float x[], float y[], float z[],
				     float samples[],
				     int len,
				     int inter_type) 
{
  int 
    c;
  Real 
    val[MAX_DIMENSIONS];
  
  
  for_inclusive(c,1,len) {
     val[0] = 0.0;
     evaluate_volume_in_world(data,
                              (Real)x[c], (Real)y[c], (Real)z[c], 
                              inter_type,
                              TRUE,
                              0.0, 
                              val,
                              NULL, NULL, NULL, 
                              NULL, NULL, NULL, 
                              NULL, NULL, NULL );
     
     samples[c] = (float)val[0];
  }
  
}

/*********************************************************************** 
   use the list of voxel coordinates stored in x[], y[], z[] and the
   voxel offset stored in dx, dy, dz to interpolate len samples from
   the volume 'data' 

   return the value of the (pseudo) normalized objective function
   (indicated by obf_func) between the list of values in *a1 and the values
   interpolated from data.  The value returned _should_ (but is not
   garenteed) to be between -1.0 and 1.0 with 1.0 indicating a PERFECT FIT
   (ie, greater values indicate better fits)

   note: the volume is assumed to be in x,y,z order.

   actually, the order does not matter, except that dx corresponds to
   the displacement along the 1st dimension x[], dy corresponds to
   second dimension y[] and dz to z[], the last.  When doing 2D
   processing, the first dimension (x[]) is assumed to be the slowest
   varying, and the deformation in dx=0.

   CAVEAT 1: only nearest neighbour and tri-linear interpolation
             are supported.
   CAVEAT 2: only Volume data types of UNSIGNED_BYTE, SIGNED_SHORT, and
             UNSIGNED_SHORT are supported.

*/

public float go_get_samples_with_offset(
     Volume data,                  /* The volume of data */
     float *x, float *y, float *z, /* the positions of the sub-lattice */
     Real  dx, Real  dy, Real dz,  /* the local displacement to apply  */
     int obj_func,                 /* the type of obj function req'd   */
     int len,                      /* number of sub-lattice nodes      */
     float sqrt_s1,                /* normalization factor for obj func*/
     float *a1,                    /* feature value for (x,y,z) nodes  */
     BOOLEAN use_nearest_neighbour) /* interpolation flag              */
{
  float
    tmp,
    sample, r,
    s1,s3;			/* to store the sums for f1,f2,f3 */
  int 
    sizes[3],
    ind0, ind1, ind2, 
    c,number_of_nonzero_samples;  
  int xs,ys,zs;
  float
    f_trans, f_scale;

  static double v0, v1, v2;
  static double f0, f1, f2, r0, r1, r2, r1r2, r1f2, f1r2, f1f2;
  static double v000, v001, v010, v011, v100, v101, v110, v111;

  unsigned char  ***byte_ptr;
  unsigned short ***ushort_ptr;
           short ***sshort_ptr;

  Data_types 
    data_type;
  

  number_of_nonzero_samples = 0;

  data_type = get_volume_data_type(data);
  get_volume_sizes(data, sizes);  
  xs = sizes[0];  
  ys = sizes[1];  
  zs = sizes[2];

  f_trans = data->real_value_translation;
  f_scale = data->real_value_scale;


  s1 = 0.0;
  s3 = 0.0;

  ++a1;                         /* inc pointer, so that we are pointing to
                                   the first feature value, corresponding
                                   to the first sub-lattice point x,y,z   */


  if (use_nearest_neighbour) {
				/* then do fast NN interpolation */

    dx += 0.5;			/* to achieve `rounding' for ind0, ind1 and */
    dy += 0.5;			/* ind2 below */
    dz += 0.5;
    
    switch( data_type ) {
    case UNSIGNED_BYTE:  
      
      byte_ptr = VOXEL_DATA (data); 

      ++x; ++y; ++z; 
      
      for_inclusive(c,1,len) {
	
	/*
	   this is code to test timing of David's evaluate_volume_in_world()
	   interpolation code.  While it is _VERY_ general, the eval_vol_in_world()'s
	   NN interpolation is approximately 8-9 times slower than the bit of 
	   fast NN code below.
	   
	   Real
	   sampleR,index[5], wx, wy,wz;
	   
	   index[0] = (Real) ( *x++ + dx );
	   index[1] = (Real) ( *y++ + dy );
	   index[2] = (Real) ( *z++ + dz );
	   index[3] = 0.0;
	   index[4] = 0.0;
	   
	   convert_voxel_to_world(data, index, &wx, &wy, &wz );
	   
	   evaluate_volume_in_world(data, wx, wy, wz,
				    0, TRUE, 0.0, &sampleR,
				    NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
	   sample = (float)sampleR;
	*/
	
      if (  !(obj_func==NONLIN_CHAMFER) ||  (*a1>0) ) {

         ind0 = (int) ( *x++ + dx );
         ind1 = (int) ( *y++ + dy );
         ind2 = (int) ( *z++ + dz );
         
         if (ind0>=0 && ind0<xs &&
             ind1>=0 && ind1<ys &&
             ind2>=0 && ind2<zs) {
            sample = (float)(byte_ptr[ind0][ind1][ind2]) * f_scale + f_trans;
         }
         else
            sample = 0.0;
         
         
#include "switch_obj_func.c"
      }
      else {
         a1++;
         x++;
         y++;
         z++;
      }
    }
    break;
    case SIGNED_SHORT:  
      
      sshort_ptr = VOXEL_DATA (data); 

      ++x; ++y; ++z; 
      
      for_inclusive(c,1,len) {
	
         if (  !(obj_func==NONLIN_CHAMFER) ||  (*a1>0) ) {
            ind0 = (int) ( *x++ + dx );
            ind1 = (int) ( *y++ + dy );
            ind2 = (int) ( *z++ + dz );
            
            if (ind0>=0 && ind0<xs &&
                ind1>=0 && ind1<ys &&
                ind2>=0 && ind2<zs) {
               sample = sshort_ptr[ind0][ind1][ind2] * f_scale + f_trans;
            }
            else
               sample = 0.0;
	
#include "switch_obj_func.c"
         }
         else {
            a1++;
            x++;
            y++;
            z++;
         }
      }
      break;
    case UNSIGNED_SHORT:  
      
      ushort_ptr = VOXEL_DATA (data); 


      ushort_ptr = (unsigned short ***)data->array.data;

      ++x; ++y; ++z; 
      
      for_inclusive(c,1,len) {
	
         if (  !(obj_func==NONLIN_CHAMFER) ||  (*a1>0) ) {
            ind0 = (int) ( *x++ + dx );
            ind1 = (int) ( *y++ + dy );
            ind2 = (int) ( *z++ + dz );
            
            if (ind0>=0 && ind0<xs &&
                ind1>=0 && ind1<ys &&
                ind2>=0 && ind2<zs) {
               sample = (float)(ushort_ptr[ind0][ind1][ind2]) * f_scale + f_trans;
            }
            else
               sample = 0.0;
            
#include "switch_obj_func.c"
         }
         else {
            a1++;
            x++;
            y++;
            z++;
         }
      }
      break;
    default:
      print_error_and_line_num("Data type not supported in go_get_samples_with_offset (only signed_byte, signed_short, unsigned_short allowed)",__FILE__, __LINE__);
    }
    
  }
  else {			/* then do fast trilinear interpolation */
    
     switch( data_type ) {
        case UNSIGNED_BYTE:  
        
        byte_ptr = VOXEL_DATA (data);
        
        ++x; ++y; ++z; 
        
        for_inclusive(c,1,len) {
           
           /*  fast tri-linear interplation */
           
           if (  !(obj_func==NONLIN_CHAMFER) ||  (*a1>0) ) {
              v0 = (Real) ( *x++ + dx );
              v1 = (Real) ( *y++ + dy );
              v2 = (Real) ( *z++ + dz );
              
              ind0 = (int)v0;
              ind1 = (int)v1;
              ind2 = (int)v2;
              
              if (ind0>=0 && ind0<(xs-1) &&
                  ind1>=0 && ind1<(ys-1) &&
                  ind2>=0 && ind2<(zs-1)) {
                 
                 /* get the data */
                 v000 = (Real)(byte_ptr[ind0  ][ind1  ][ind2  ]);
                 v001 = (Real)(byte_ptr[ind0  ][ind1  ][ind2+1]);
                 v010 = (Real)(byte_ptr[ind0  ][ind1+1][ind2  ]);
                 v011 = (Real)(byte_ptr[ind0  ][ind1+1][ind2+1]);
                 v100 = (Real)(byte_ptr[ind0+1][ind1  ][ind2  ]);
                 v101 = (Real)(byte_ptr[ind0+1][ind1  ][ind2+1]);
                 v110 = (Real)(byte_ptr[ind0+1][ind1+1][ind2  ]);
                 v111 = (Real)(byte_ptr[ind0+1][ind1+1][ind2+1]);
                 
                 /* Get the fraction parts */
                 f0 = v0 - ind0;
                 f1 = v1 - ind1;
                 f2 = v2 - ind2;
                 r0 = 1.0 - f0;
                 r1 = 1.0 - f1;
                 r2 = 1.0 - f2;
                 
                 /* Do the interpolation */
                 r1r2 = r1 * r2;
                 r1f2 = r1 * f2;
                 f1r2 = f1 * r2;
                 f1f2 = f1 * f2;
                 
                 sample   = 
                    r0 *  (r1r2 * v000 +
                           r1f2 * v001 +
                           f1r2 * v010 +
                           f1f2 * v011);
                 sample  +=
                    f0 *  (r1r2 * v100 +
                           r1f2 * v101 +
                           f1r2 * v110 +
                           f1f2 * v111);
                 
                 sample = sample * f_scale + f_trans;
                 
                 
              }
              else
                 sample = 0.0;
              
              
#include "switch_obj_func.c"
           }
           else {
              a1++;
              x++;
              y++;
              z++;
           }
           
        }
        break;
        case SIGNED_SHORT:  
        
        sshort_ptr = VOXEL_DATA (data); 
        
        ++x; ++y; ++z; 
        
        for_inclusive(c,1,len) {
           
	/*  fast tri-linear interpolation */
	
           if (  !(obj_func==NONLIN_CHAMFER) ||  (*a1>0) ) {
              v0 = (Real) ( *x++ + dx );
              v1 = (Real) ( *y++ + dy );
              v2 = (Real) ( *z++ + dz );
              
              ind0 = (int)v0;
              ind1 = (int)v1;
              ind2 = (int)v2;
              
              if (ind0>=0 && ind0<(xs-1) &&
                  ind1>=0 && ind1<(ys-1) &&
                  ind2>=0 && ind2<(zs-1)) {
                 
                 /* get the data */
                 v000 = (Real)(sshort_ptr[ind0  ][ind1  ][ind2  ]);
                 v001 = (Real)(sshort_ptr[ind0  ][ind1  ][ind2+1]);
                 v010 = (Real)(sshort_ptr[ind0  ][ind1+1][ind2  ]);
                 v011 = (Real)(sshort_ptr[ind0  ][ind1+1][ind2+1]);
                 v100 = (Real)(sshort_ptr[ind0+1][ind1  ][ind2  ]);
                 v101 = (Real)(sshort_ptr[ind0+1][ind1  ][ind2+1]);
                 v110 = (Real)(sshort_ptr[ind0+1][ind1+1][ind2  ]);
                 v111 = (Real)(sshort_ptr[ind0+1][ind1+1][ind2+1]);
                 
                 /* Get the fraction parts */
                 f0 = v0 - ind0;
                 f1 = v1 - ind1;
                 f2 = v2 - ind2;
                 r0 = 1.0 - f0;
                 r1 = 1.0 - f1;
                 r2 = 1.0 - f2;
                 
                 /* Do the interpolation */
                 r1r2 = r1 * r2;
                 r1f2 = r1 * f2;
                 f1r2 = f1 * r2;
                 f1f2 = f1 * f2;
                 
                 sample   = 
                    r0 *  (r1r2 * v000 +
                           r1f2 * v001 +
                           f1r2 * v010 +
                           f1f2 * v011);
                 sample  +=
                    f0 *  (r1r2 * v100 +
                           r1f2 * v101 +
                           f1r2 * v110 +
                           f1f2 * v111);
                 
                 sample = sample * f_scale + f_trans;
                 
              }
              else
                 sample = 0.0;
              
#include "switch_obj_func.c"
           }
           else {
              a1++;
              x++;
              y++;
              z++;
           }
           
        }
        break;
        case UNSIGNED_SHORT:  
        
        ushort_ptr = VOXEL_DATA (data);
        
        ++x; ++y; ++z; 
        
        for_inclusive(c,1,len) {
           
           if (  !(obj_func==NONLIN_CHAMFER) ||  (*a1>0) ) {
              /*  fast tri-linear interplation */
              
              v0 = (Real) ( *x++ + dx );
              v1 = (Real) ( *y++ + dy );
              v2 = (Real) ( *z++ + dz );
              
              ind0 = (int)v0;
              ind1 = (int)v1;
              ind2 = (int)v2;
              
              if (ind0>=0 && ind0<(xs-1) &&
                  ind1>=0 && ind1<(ys-1) &&
                  ind2>=0 && ind2<(zs-1)) {
                 
                 /* get the data */
                 v000 = (Real)(ushort_ptr[ind0  ][ind1  ][ind2  ]);
                 v001 = (Real)(ushort_ptr[ind0  ][ind1  ][ind2+1]);
                 v010 = (Real)(ushort_ptr[ind0  ][ind1+1][ind2  ]);
                 v011 = (Real)(ushort_ptr[ind0  ][ind1+1][ind2+1]);
                 v100 = (Real)(ushort_ptr[ind0+1][ind1  ][ind2  ]);
                 v101 = (Real)(ushort_ptr[ind0+1][ind1  ][ind2+1]);
                 v110 = (Real)(ushort_ptr[ind0+1][ind1+1][ind2  ]);
                 v111 = (Real)(ushort_ptr[ind0+1][ind1+1][ind2+1]);
                 
                 /* Get the fraction parts */
                 f0 = v0 - ind0;
                 f1 = v1 - ind1;
                 f2 = v2 - ind2;
                 r0 = 1.0 - f0;
                 r1 = 1.0 - f1;
                 r2 = 1.0 - f2;
                 
                 /* Do the interpolation */
                 r1r2 = r1 * r2;
                 r1f2 = r1 * f2;
                 f1r2 = f1 * r2;
                 f1f2 = f1 * f2;
                 
                 sample   = 
                    r0 *  (r1r2 * v000 +
                           r1f2 * v001 +
                           f1r2 * v010 +
                           f1f2 * v011);
                 sample  +=
                    f0 *  (r1r2 * v100 +
                           r1f2 * v101 +
                           f1r2 * v110 +
                           f1f2 * v111);
                 
                 sample = sample * f_scale + f_trans;
                 
              }
              else
                 sample = 0.0;
              
#include "switch_obj_func.c"
           }
           else {
              a1++;
              x++;
              y++;
              z++;
           }

        }
        break;
        default:
        print_error_and_line_num("Data type not supported in go_get_samples_with_offset (only signed_byte, signed_short, unsigned_short allowed)",__FILE__, __LINE__);
     }
     
  }
  
  
				/* do the last bits of the similarity function
				   calculation here - normalizing each obj_func
                                   where-ever possible: */
  r = 0.0;
  switch (obj_func) {

  case NONLIN_XCORR:            /* use standard normalized cross-correlation 
                                   where 0.0 < r < 1.0, where 1.0 is best*/
    if ( sqrt_s1 < 0.001 && s3 < 0.00001) {
      r = 1.0;
    }
    else {
      if ( sqrt_s1 < 0.001 || s3 < 0.00001) {
	r = 0.0;
      }
      else {
	r = s1 / (sqrt_s1*sqrt((double)s3));
      }
    }
    /* r = 1.0 - r;                 now, 0 is best                   */
    break;

  case NONLIN_DIFF:             /* sqrt_s1 stores the number of samples in
                                   the sub-lattice 
                                   s1 stores the sum of the magnitude of
                                   the differences*/

    r = -s1 / sqrt_s1;	        /* r = average intensity difference ; with
                                   -max(intensity range) < r < 0,
                                   where 0 is best                  */
    break;
  case NONLIN_LABEL:
    r = s1 / sqrt_s1;           /* r = average label agreement,
                                   s1 stores the number of similar labels
                                   0 < r < 1.0                      
                                   where 1.0 is best                */
    break;
  case NONLIN_CHAMFER:
    if (number_of_nonzero_samples>0) {
       r = 1.0 - (s1 / (20.0*number_of_nonzero_samples));	
                                /* r = 1- average distance / 20mm 
                                       0 < r < ~1.0 
                                   where 1.0 is best     
                                       and where 2.0cm is an arbitrary value to
                                       norm the dist, corresponding to a guess
                                       at the maximum average cortical variability

                                       so the max(r) could be greater than
                                       1.0, but when it is, shouldn't
                                       chamfer have larger weight to drive
                                       the fit? */
    }
    else
       r = 2.0;                 /* this is simply a value > 1.5, used as a
                                   flag to indicate that there were no
                                   samples used for the chamfer */
    break;
  default:
    print_error_and_line_num("Objective function %d not supported in go_get_samples_with_offset",__FILE__, __LINE__,obj_func);
  }
  
  
  return(r);
}

/* Build the target lattice by transforming the source points through the
   current non-linear transformation stored in:

        Gglobals->trans_info.transformation                        

   both input (px,py,pz) and output (tx,ty,tz) coordinate lists are in
   WORLD COORDINATES

*/
public void    build_target_lattice(float px[], float py[], float pz[],
				    float tx[], float ty[], float tz[],
				    int len, int dim)
{
  int i;
  Real x,y,z;


  for_inclusive(i,1,len) {

    if (dim==3)
      general_transform_point(Gglobals->trans_info.transformation, 
			      (Real)px[i],(Real) py[i], (Real)pz[i], 
			      &x, &y, &z);
    else
      general_transform_point_in_trans_plane(Gglobals->trans_info.transformation, 
				 (Real)px[i],(Real) py[i], (Real)pz[i], 
				 &x, &y, &z);
    
    tx[i] = (float)x;
    ty[i] = (float)y;
    tz[i] = (float)z;
  }
}

/* Build the target lattice by transforming the source points through the
   current non-linear transformation stored in:

        Glinear_transform and Gsuper_d{x,y,z} deformation volumes  

   both input (px,py,pz) and output (tx,ty,tz) coordinate lists are in
   WORLD COORDINATES

*/
public void    build_target_lattice_using_super_sampled_def(
                                     float px[], float py[], float pz[],
				     float tx[], float ty[], float tz[],
				     int len, int dim)
{
  int 
    i,j,
    sizes[MAX_DIMENSIONS],
    xyzv[MAX_DIMENSIONS];
  Real 
    def_vector[N_DIMENSIONS],
    voxel[MAX_DIMENSIONS],
    x,y,z;
  long 
    index[MAX_DIMENSIONS];

  get_volume_sizes(Gsuper_sampled_vol,sizes);
  get_volume_XYZV_indices(Gsuper_sampled_vol,xyzv);

  for_inclusive(i,1,len) {

				/* apply linear part of the transformation */

    general_transform_point(Glinear_transform,
			    (Real)px[i], (Real)py[i], (Real)pz[i], 
			    &x, &y, &z);

				/* now get the non-linear part, using
                                   nearest neighbour interpolation in
                                   the super-sampled deformation
                                   volume. */

    convert_world_to_voxel(Gsuper_sampled_vol, 
			   x,y,z, voxel);

    if ((voxel[ xyzv[X] ] >= -0.5) && (voxel[ xyzv[X] ] < sizes[xyzv[X]]-0.5) &&
	(voxel[ xyzv[Y] ] >= -0.5) && (voxel[ xyzv[Y] ] < sizes[xyzv[Y]]-0.5) &&
	(voxel[ xyzv[Z] ] >= -0.5) && (voxel[ xyzv[Z] ] < sizes[xyzv[Z]]-0.5) ) {

      for_less(j,0,3) index[ xyzv[j] ] = (long) (voxel[ xyzv[j] ]+0.5);
      
      for_less(index[ xyzv[Z+1] ],0, sizes[ xyzv[Z+1] ]) 
	GET_VALUE_4D(def_vector[ index[ xyzv[Z+1] ]  ], \
		     Gsuper_sampled_vol, \
		     index[0], index[1], index[2], index[3]);


      x += def_vector[X];
      y += def_vector[Y];
      z += def_vector[Z];
    }

    tx[i] = (float)x;
    ty[i] = (float)y;
    tz[i] = (float)z;

  }
}
