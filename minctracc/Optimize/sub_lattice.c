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
@VERSION    : $Id: sub_lattice.c,v 1.16 2009-04-03 18:36:59 louis Exp $
@MODIFIED   : $Log: sub_lattice.c,v $
@MODIFIED   : Revision 1.16  2009-04-03 18:36:59  louis
@MODIFIED   : made changes to use only DOUBLES for input source and model volumes, and for all estimation of deformation fields
@MODIFIED   :
@MODIFIED   : Revision 1.15  2006/11/30 09:07:32  rotor
@MODIFIED   :  * many more changes for clean minc 2.0 build
@MODIFIED   :
@MODIFIED   : Revision 1.14  2006/11/29 09:09:34  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 1.13  2005/07/20 20:45:51  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 1.12  2005/07/18 19:14:02  rotor
@MODIFIED   :  * Optimisations to code resulting in 30% speed increase for nonlinear fitting
@MODIFIED   :
@MODIFIED   : Revision 1.11  2005/06/28 18:56:18  rotor
@MODIFIED   :  * added masking for feature volumes (irina and patricia)
@MODIFIED   :
@MODIFIED   : Revision 1.10  2004/02/12 06:08:21  rotor
@MODIFIED   :  * removed /static defs
@MODIFIED   :
@MODIFIED   : Revision 1.9  2003/02/26 00:56:38  lenezet
@MODIFIED   : for 2D : now computes all 3 coordinates for the "start" (to take into account the slice position).
@MODIFIED   : simplification of build_lattices.
@MODIFIED   : bug correction in amoeba_NL_obj_function.
@MODIFIED   :
@MODIFIED   : Revision 1.8  2003/02/04 06:08:46  stever
@MODIFIED   : Add support for correlation coefficient and sum-of-squared difference.
@MODIFIED   :
@MODIFIED   : Revision 1.7  2002/12/13 21:18:20  lenezet
@MODIFIED   :
@MODIFIED   : A memory leak has been repaired
@MODIFIED   :
@MODIFIED   : Revision 1.6  2002/08/09 18:31:33  stever
@MODIFIED   : Fix error messages: it is UNsigned bytes that are allowed.
@MODIFIED   :
@MODIFIED   : Revision 1.5  2002/03/26 14:15:46  stever
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
#include <volume_io.h>        

#include <Proglib.h>
#include "constants.h"
#include "arg_data.h"                /* definition of the global data struct      */
#include "sub_lattice.h"
#include "init_lattice.h"


extern Arg_Data *Gglobals;      /* defined in do_nonlinear.c */
extern VIO_Volume   Gsuper_sampled_vol; /* defined in do_nonlinear.c */
extern VIO_General_transform 
                *Glinear_transform;/* defined in do_nonlinear.c */

                                /* prototypes for functions used here: */

extern float
  *SX, *SY, *SZ;

 void  general_transform_point_in_trans_plane(
    VIO_General_transform   *transform,
    VIO_Real                x,
    VIO_Real                y,
    VIO_Real                z,
    VIO_Real                *x_transformed,
    VIO_Real                *y_transformed,
    VIO_Real                *z_transformed );

int point_not_masked(VIO_Volume volume, 
                            VIO_Real wx, VIO_Real wy, VIO_Real wz);

int voxel_point_not_masked(VIO_Volume volume, 
                                  VIO_Real vx, VIO_Real vy, VIO_Real vz);


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

/* make sure that the lattice is defined on the axis of the 2nd data
   volume (and hence the grid transform volume) and NOT the world
   coordinate axis! */

void    build_source_lattice(VIO_Real x, VIO_Real y, VIO_Real z,
                                    float PX[], float PY[], float PZ[],
                                    VIO_Real width_x, VIO_Real width_y, VIO_Real width_z, 
                                    int nx, int ny, int nz,
                                    int ndim, int *length)
{
  int 
    c, 
    tnx, tny, tnz,
    i,j,k;
  float 
    radius_squared,
    tx,ty,tz;
  float
    abs_step,
    dir[3][3];

  *length = 0; 
  c = 1;
  radius_squared = 0.55 * 0.55;        /* a bit bigger than .5^2 */
  

  for(i=0; i<3; i++) {
    
    abs_step = fabs(Gglobals->step[i]);

    dir[i][0] = Point_x(Gglobals->directions[i]) / abs_step;
    dir[i][1] = Point_y(Gglobals->directions[i]) / abs_step;
    dir[i][2] = Point_z(Gglobals->directions[i]) / abs_step;
  }


  if (Gglobals->count[0] > 1) { tnx = nx; }  else {    tnx = 1;  }
  if (Gglobals->count[1] > 1) { tny = ny; }  else {    tny = 1;  }
  if (Gglobals->count[2] > 1) { tnz = nz; }  else {    tnz = 1;  }


  for(i=0; i<tnx; i++) {
    for(j=0; j<tny; j++) {
      for(k=0; k<tnz; k++) {

        if (tnx>1) { tx = -0.5 + (float)(i)/(float)(tnx-1); } else tx = 0.0;
        if (tny>1) { ty = -0.5 + (float)(j)/(float)(tny-1); } else ty = 0.0;
        if (tnz>1) { tz = -0.5 + (float)(k)/(float)(tnz-1); } else tz = 0.0;
      

        if ((tx*tx + ty*ty + tz*tz) <= radius_squared) {

          tx *= width_x;
          ty *= width_y;
          tz *= width_z;

          PX[c] = (float)x + tx*dir[VIO_X][VIO_X] + ty*dir[VIO_Y][VIO_X] + tz*dir[VIO_Z][VIO_X];
          PY[c] = (float)y + tx*dir[VIO_X][VIO_Y] + ty*dir[VIO_Y][VIO_Y] + tz*dir[VIO_Z][VIO_Y];
          PZ[c] = (float)z + tx*dir[VIO_X][VIO_Z] + ty*dir[VIO_Y][VIO_Z] + tz*dir[VIO_Z][VIO_Z];

          c++;
          (*length)++;
        }
      }
    }
  }


  /*
  if (ndim==2) {
    for(i=0; i<nx; i++)
      for(j=0; j<ny; j++) {
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
    for(i=0; i<nx; i++)
      for(j=0; j<ny; j++)
        for(k=0; k<nz; k++) {
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
  */

}

/*********************************************************************** 
   use the world coordinates stored in x[],y[],z[] to interpolate len
   samples from the volume 'data' */

void go_get_samples_in_source(VIO_Volume data, VIO_Volume mask,
                                     float x[], float y[], float z[],
                                     float samples[], VIO_BOOL masked_samples[],
                                     int len,
                                     int inter_type) 
{
  int 
    c;
  VIO_Real 
    val[VIO_MAX_DIMENSIONS];
  
  for(c=1; c<=len; c++) {  
    if (point_not_masked(mask, (VIO_Real)x[c], (VIO_Real)y[c], (VIO_Real)z[c])){
      masked_samples[c] = FALSE;
      val[0] = 0.0;
      evaluate_volume_in_world(data,
                               (VIO_Real)x[c], (VIO_Real)y[c], (VIO_Real)z[c], 
                               inter_type,
                               TRUE,
                               0.0, 
                               val,
                               NULL, NULL, NULL, 
                               NULL, NULL, NULL, 
                               NULL, NULL, NULL );
      
      samples[c] = (float)val[0];
    }
    else
      {
        masked_samples[c] = TRUE;
        samples[c] = 0.0;
      }
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

   note: the volume is assumed to be in MIzspace, MIyspace, MIxspace
   order. (since the data is loaded in with default_dim_names as an
   option to input_volume).

   actually, the order does not matter, except that dx corresponds to
   the displacement along the 1st dimension x[], dy corresponds to
   second dimension y[] and dz to z[], the last.  When doing 2D
   processing, the first dimension (x[]) is assumed to be the slowest
   varying, and the deformation in dx=0.

   CAVEAT 1: only nearest neighbour and tri-linear interpolation
             are supported.

	     *** AJ + LC: 3/17/2009:  only DOUBLE data now supported.
   CAVEAT 2: only VIO_Volume data types of UNSIGNED_BYTE, SIGNED_SHORT, and
             UNSIGNED_SHORT are supported.

*/

float go_get_samples_with_offset(
				 VIO_Volume data,                  /* The volume of data */
				 VIO_Volume mask,                  /* The target mask */  
				 float *x, float *y, float *z,     /* the positions of the sub-lattice */
				 VIO_Real  dx, VIO_Real  dy, VIO_Real dz,  /* the local displacement to apply  */
				 int obj_func,                     /* the type of obj function req'd   */
				 int len,                          /* number of sub-lattice nodes      */
				 int *sample_target_count,         /* number of nonmasked  sub-lattice nodes */
				 float normalization,              /* normalization factor for obj func*/
				 float *a1,                        /* feature value for (x,y,z) nodes  */
				 VIO_BOOL *m1,                     /* mask flag for (x,y,z) nodes in source */ 
				 VIO_BOOL use_nearest_neighbour)   /* interpolation flag              */
{
  double
    sample, r,
    s1,s2,s3,s4,s5,tmp;                   /* accumulators for inner loop */
  int 
    sizes[3],
    ind0, ind1, ind2, 
    offset0, offset1, offset2,
    c,number_of_nonzero_samples;  
  int xs,ys,zs;
  float
    f_trans, f_scale;

  static double v0, v1, v2;
  static double f0, f1, f2, r0, r1, r2, r1r2, r1f2, f1r2, f1f2;
  static double v000, v001, v010, v011, v100, v101, v110, v111;

  double ***double_ptr;
  
  double mean_s = 0.0;		/* init variables for stats */
  double mean_t = 0.0;
  double var_s = 0.0;
  double var_t = 0.0;
  double covariance = 0.0;
  
  number_of_nonzero_samples = 0;

  get_volume_sizes(data, sizes);  
  xs = sizes[0];  
  ys = sizes[1];  
  zs = sizes[2];


  s1 = s2 = s3 = s4 = s5 = 0.0;
  ++a1;   ++m1;                   /* inc pointer, so that we are pointing to
                                   the first feature value, corresponding
                                   to the first sub-lattice point x,y,z   */



  if (use_nearest_neighbour) {
                                /* then do fast NN interpolation */

    dx += 0.;                        /* to achieve `rounding' for ind0, ind1 and */
    dy += 0.;                        /* ind2 below */
    dz += 0.;
    
    double_ptr = VOXEL_DATA (data); 
      
    /* increment by one as we index from 1...n (but the array is ALLOC'd from 0..n) */
    x++; y++; z++; 
      
    /* for each sub-lattice node */
    for(c=len; c--;) {
        
        /*
         *  this is code to test timing of David's evaluate_volume_in_world()
         *  interpolation code.  While it is _VERY_ general, the eval_vol_in_world()'s
         *  NN interpolation is approximately 8-9 times slower than the bit of 
         *  fast NN code below.
         *  
         *  VIO_Real sampleR, index[5], wx, wy, wz;
         *  index[0] = (VIO_Real) ( *x + dx );
         *  index[1] = (VIO_Real) ( *y + dy );
         *  index[2] = (VIO_Real) ( *z + dz );
         *  index[3] = index[4] = 0.0;
         *  convert_voxel_to_world(data, index, &wx, &wy, &wz );
         *  evaluate_volume_in_world(data, wx, wy, wz,
         *                            0, TRUE, 0.0, &sampleR,
         *                            NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
         *   sample = (double)sampleR;
         */
        
        if (!(*m1) && (voxel_point_not_masked(mask, (VIO_Real)*x, (VIO_Real)*y, (VIO_Real)*z)) && (!(obj_func==NONLIN_CHAMFER) ||  (*a1>0)) ) {

         ind0 = (int) ( *x + dx );
         ind1 = (int) ( *y + dy );
         ind2 = (int) ( *z + dz );
         
         if (ind0>=0 && ind0<xs &&
             ind1>=0 && ind1<ys &&
             ind2>=0 && ind2<zs) {
	   sample = (double)(double_ptr[ind0][ind1][ind2]);
         }
         else{
            sample = 0.0;
         }

#include "switch_obj_func.c"	/* contains a case statement to do the sample-to-sample computations required for each possible non-lin objective function */
      }
        
       /* increment the coordinate, value and mask array pointers */  
       x++;			/* x,y,z are the coords of the sublattice in the fixed (source) image */
       y++;
       z++;
       a1++;			/* a1 is from the fixed image */
       m1++;			/* m1 is rom the mask on the fixed image */
     
    }
  }
  else {                        /* then do fast trilinear interpolation */
    
    /* set up offsets */
    offset0 = (Gglobals->count[VIO_Z] > 1) ? 1 : 0;
    offset1 = (Gglobals->count[VIO_Y] > 1) ? 1 : 0;
    offset2 = (Gglobals->count[VIO_X] > 1) ? 1 : 0;
    
    double_ptr = VOXEL_DATA (data);
        
    /* increment by one as we index from 1...n */
    ++x; ++y; ++z; 
        
    /* for each sub-lattice node */
    for(c=len; c--;) {
           
           /*  fast tri-linear interplation */
           
      if  (  !(*m1) && (voxel_point_not_masked(mask, (VIO_Real)*x, (VIO_Real)*y, (VIO_Real)*z)) && (!(obj_func==NONLIN_CHAMFER) ||  (*a1>0)) ) {
              v0 = (VIO_Real) ( *x + dx );
              v1 = (VIO_Real) ( *y + dy );
              v2 = (VIO_Real) ( *z + dz );
              
              ind0 = (int)v0;
              ind1 = (int)v1;
              ind2 = (int)v2;
              
              if (ind0>=0 && ind0<(xs-offset0) &&
                  ind1>=0 && ind1<(ys-offset1) &&
                  ind2>=0 && ind2<(zs-offset2)) {
                 
                 /* get the data */
                 v000 = (VIO_Real)(double_ptr[ind0        ][ind1        ][ind2        ]);
                 v001 = (VIO_Real)(double_ptr[ind0        ][ind1        ][ind2+offset2]);
                 v010 = (VIO_Real)(double_ptr[ind0        ][ind1+offset1][ind2        ]);
                 v011 = (VIO_Real)(double_ptr[ind0        ][ind1+offset1][ind2+offset2]);
                 v100 = (VIO_Real)(double_ptr[ind0+offset0][ind1        ][ind2        ]);
                 v101 = (VIO_Real)(double_ptr[ind0+offset0][ind1        ][ind2+offset2]);
                 v110 = (VIO_Real)(double_ptr[ind0+offset0][ind1+offset1][ind2        ]);
                 v111 = (VIO_Real)(double_ptr[ind0+offset0][ind1+offset1][ind2+offset2]);
                 
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
                 
                 sample = sample;
                 
                 
              }
              else{
		sample = 0.0;
	      }
              
              
#include "switch_obj_func.c"

      }
      
      /* increment the coordinate, value and mask array pointers */  
      x++;			/* x,y,z are the coords of the sublattice in the fixed (source) image */
      y++;
      z++;
      a1++;			/* a1 is from the fixed image */
      m1++;			/* m1 is rom the mask on the fixed image */
      
    } 
  


                                /* do the last bits of the similarity function
                                   calculation here - normalizing each obj_func
                                   where-ever possible: */
  r = 0.0;

  switch (obj_func) {

  case NONLIN_XCORR:            /* use standard normalized cross-correlation 
                                   where 0.0 < r < 1.0, where 1.0 is best*/
    if ( normalization < 0.001 && s3 < 0.00001) {
      r = 1.0;
    }
    else {
      if ( normalization < 0.001 || s3 < 0.00001) {
        r = 0.0;
      }
      else {
        r = s1 / ((sqrt((double)s2))*(sqrt((double)s3)));
      }
    }
    /* r = 1.0 - r;                 now, 0 is best                   */
    break;

  case NONLIN_DIFF:             /* normalization stores the number of samples in
                                   the sub-lattice 
                                   s1 stores the sum of the magnitude of
                                   the differences*/

     r = -s1 /number_of_nonzero_samples;        /* r = average intensity difference ; with
                                   -max(intensity range) < r < 0,
                                   where 0 is best                  */
    break;
  case NONLIN_LABEL:
     r = s1 /number_of_nonzero_samples;           /* r = average label agreement,
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
  case NONLIN_CORRCOEFF:
      {
          /* Accumulators:
           * s1 = sum of source image values
           * s2 = sum of target image values
           * s3 = sum of squared source image values
           * s4 = sum of squared target image values
           * s5 = sum of source*target values
           *
           * normalization = #values considered
           */
          if (number_of_nonzero_samples>0) {
            mean_s = s1 / number_of_nonzero_samples;
            mean_t = s2 / number_of_nonzero_samples;
            var_s = s3 / number_of_nonzero_samples - mean_s*mean_s;
            var_t = s4 / number_of_nonzero_samples - mean_t*mean_t;
            covariance = s5 / number_of_nonzero_samples - mean_s*mean_t;
          }
          else {
            mean_s = 0.0;
            mean_t = 0.0;
            var_s = 0.0;
            var_t = 0.0;
            covariance = 0.0;
          }


          if ((var_s < 0.00001) || (var_t < 0.00001) ) {
            r = 0.0;
          }
          else {
            r = covariance / sqrt( var_s*var_t );            
          }
      }
      break;
          
  case NONLIN_SQDIFF:           /* normalization stores the number of samples 
                                   in the sub-lattice.
                                   s1 stores the sum of the squared intensity
                                   differences */
    r = -s1 /number_of_nonzero_samples;
    break;

  default:
    print_error_and_line_num("Objective function %d not supported in go_get_samples_with_offset",__FILE__, __LINE__,obj_func);
  }
  
  
  return(r);
  }

}




/* Build the target lattice by transforming the source points through the
   current non-linear transformation stored in:

        Gglobals->trans_info.transformation                        

   both input (px,py,pz) and output (tx,ty,tz) coordinate lists are in
   WORLD COORDINATES

*/
void    build_target_lattice(float px[], float py[], float pz[],
			     float tx[], float ty[], float tz[],
			     int len, int dim)
{
  int i;
  VIO_Real x,y,z;


  for(i=1; i<=len; i++) {

    general_transform_point(Gglobals->trans_info.transformation, 
                            (VIO_Real)px[i],(VIO_Real) py[i], (VIO_Real)pz[i], 
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
void    build_target_lattice_using_super_sampled_def(
                                     float px[], float py[], float pz[],
                                     float tx[], float ty[], float tz[],
                                     int len, int dim)
{
  int 
    i,j,
    sizes[VIO_MAX_DIMENSIONS],
    xyzv[VIO_MAX_DIMENSIONS];
  VIO_Real 
    def_vector[VIO_N_DIMENSIONS],
    voxel[VIO_MAX_DIMENSIONS],
    x,y,z;
  long 
    index[VIO_MAX_DIMENSIONS];

  get_volume_sizes(Gsuper_sampled_vol,sizes);
  get_volume_XYZV_indices(Gsuper_sampled_vol,xyzv);

  for(i=1; i<=len; i++) {

                                /* apply linear part of the transformation */

    general_transform_point(Glinear_transform,
                            (VIO_Real)px[i], (VIO_Real)py[i], (VIO_Real)pz[i], 
                            &x, &y, &z);

                                /* now get the non-linear part, using
                                   nearest neighbour interpolation in
                                   the super-sampled deformation
                                   volume. */

    convert_world_to_voxel(Gsuper_sampled_vol, 
                           x,y,z, voxel);

    if ((voxel[ xyzv[VIO_X] ] >= -0.5) && (voxel[ xyzv[VIO_X] ] < sizes[xyzv[VIO_X]]-0.5) &&
        (voxel[ xyzv[VIO_Y] ] >= -0.5) && (voxel[ xyzv[VIO_Y] ] < sizes[xyzv[VIO_Y]]-0.5) &&
        (voxel[ xyzv[VIO_Z] ] >= -0.5) && (voxel[ xyzv[VIO_Z] ] < sizes[xyzv[VIO_Z]]-0.5) ) {

      for(j=0; j<3; j++) index[ xyzv[j] ] = (long) (voxel[ xyzv[j] ]+0.5);
      
      for(index[xyzv[VIO_Z+1]]=0; index[xyzv[VIO_Z+1]]<sizes[xyzv[VIO_Z+1]]; index[xyzv[VIO_Z+1]]++) 
        GET_VALUE_4D(def_vector[ index[ xyzv[VIO_Z+1] ]  ], \
                     Gsuper_sampled_vol, \
                     index[0], index[1], index[2], index[3]);


      x += def_vector[VIO_X];
      y += def_vector[VIO_Y];
      z += def_vector[VIO_Z];
    }

    tx[i] = (float)x;
    ty[i] = (float)y;
    tz[i] = (float)z;

  }
}
