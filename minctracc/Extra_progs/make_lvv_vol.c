/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_lvv_vol
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : status
@DESCRIPTION: Program to calculate an Lvv volume from a blurred input
              volume
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Thur Oct 5 08:45:43 MET 1995
@MODIFIED   : $Log: make_lvv_vol.c,v $
@MODIFIED   : Revision 1.7  2006-11-30 09:07:31  rotor
@MODIFIED   :  * many more changes for clean minc 2.0 build
@MODIFIED   :
@MODIFIED   : Revision 1.6  2006/11/29 09:09:31  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 1.5  2005/07/20 20:45:46  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 1.4  2004/02/12 05:54:06  rotor
@MODIFIED   :  * removed public/private defs
@MODIFIED   :
@MODIFIED   : Revision 1.3  2002/12/13 21:09:12  lenezet
@MODIFIED   : tests added for direction cosines and for 2d nonlinear
@MODIFIED   :
@MODIFIED   : Revision 1.2  2002/03/26 14:15:30  stever
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
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Extra_progs/make_lvv_vol.c,v 1.7 2006-11-30 09:07:31 rotor Exp $";
#endif

#include <stdio.h>
#include <volume_io.h>
#include <Proglib.h>
#include <config.h>
#include <float.h>

#include <ParseArgv.h>

/*#include <globals.h>*/

/* Constants */
#ifndef TRUE
#  define TRUE 1
#  define FALSE 0
#endif

#define VERY_SMALL_EPS 0.0001        /* this is data dependent! */
#define sqr(a) (a)*(a)

typedef struct {
  VIO_Real 
    u,v,w,
    uu,vv,ww,
    uv,uw,vw;
} deriv_3D_struct;

static char *default_dim_names[VIO_N_DIMENSIONS] = { MIxspace, MIyspace, MIzspace };


static VIO_Real return_Lvv(VIO_Real r[3][3][3],
                       VIO_Real eps);

static void init_the_volume_to_zero(VIO_Volume volume);

static void get_volume_XYZV_indices(VIO_Volume data, int xyzv[]);

char *prog_name;
int stat_quad_total;
int stat_quad_zero;
int stat_quad_two;
int stat_quad_plus;
int stat_quad_minus;
int stat_quad_semi;

double  ftol                     = 0.005;
double  simplex_size             = 20.0;
int     iteration_limit          = 4;
double  iteration_weight         = 0.6;
double  smoothing_weight         = 0.5;
double  similarity_cost_ratio    = 0.5;
int     number_dimensions        = 3;
int     Matlab_num_steps         = 15;
int     Diameter_of_local_lattice= 5;

int     invert_mapping_flag      = FALSE;
int     clobber_flag             = FALSE;



int main(int argc, char *argv[])
{

  VIO_progress_struct
    progress;

  VIO_Status 
    stat;

  VIO_Real 
    tmp, max_val, min_val, intensity_threshold,Lvv,
    val[3][3][3];
  
  VIO_Volume 
    data, lvv;

  float ***float_vol,*f_ptr;

  int
    count,
    index[VIO_MAX_DIMENSIONS],
    data_xyzv[VIO_MAX_DIMENSIONS],
    sizes[VIO_MAX_DIMENSIONS],
    m,n,o,i,j,k,
    state;

  char 
    *history,
    output_filename[1024];

  prog_name = argv[0];
  history = time_stamp(argc, argv);

  if (argc!=3) {
    print("usage: %s input.mnc output_basename\n", prog_name);
    exit(EXIT_FAILURE);
  }

  stat = input_volume(argv[1],3,default_dim_names, NC_UNSPECIFIED, FALSE, 
                      0.0,0.0,TRUE, &data, (minc_input_options *)NULL);

  if (stat != VIO_OK) {
    print ("Error: cannot read %s.\n",argv[1]);
    exit(EXIT_FAILURE);
  }

  lvv   = copy_volume_definition(data, NC_UNSPECIFIED, FALSE, 0.0, 0.0);
  
  get_volume_sizes(data,sizes);
  get_volume_XYZV_indices(data,data_xyzv);


  max_val = -DBL_MAX;
  min_val =  DBL_MAX;
  
  VIO_ALLOC3D(float_vol  , sizes[0], sizes[1], sizes[2]);
  
  for(i=0; i<sizes[0]; i++)
    for(j=0; j<sizes[1]; j++)
      for(k=0; k<sizes[2]; k++) {
        float_vol  [i][j][k] = 0.0;
        tmp = get_volume_real_value(data, i,j,k,0,0);
        if (tmp>max_val) max_val = tmp;
        if (tmp<min_val) min_val = tmp;
      }
  intensity_threshold = 0.01 * max_val;
  print("\nintensity_threshold = %f\n",intensity_threshold);

  initialize_progress_report(&progress, FALSE, sizes[0]*sizes[1]*sizes[2]+1,
                             "Building Lvv:");
  count = 0;
  max_val = -1000000.0;
  min_val =  1000000.0;
 
  for(index[data_xyzv[VIO_X]]=1; index[data_xyzv[VIO_X]]<sizes[data_xyzv[VIO_X]]-1; index[data_xyzv[VIO_X]]++)
    for(index[data_xyzv[VIO_Y]]=1; index[data_xyzv[VIO_Y]]<sizes[data_xyzv[VIO_Y]]-1; index[data_xyzv[VIO_Y]]++)
      for(index[data_xyzv[VIO_Z]]=1; index[data_xyzv[VIO_Z]]<sizes[data_xyzv[VIO_Z]]-1; index[data_xyzv[VIO_Z]]++) {

        tmp = get_volume_real_value(data, 
                                    index[data_xyzv[VIO_X]], 
                                    index[data_xyzv[VIO_Y]], 
                                    index[data_xyzv[VIO_Z]], 0,0);

        Lvv=0.0;
     
        if (tmp > intensity_threshold) {
          for(m=-1; m<=1; m++)
            for(n=-1; n<=1; n++)
              for(o=-1; o<=1; o++)
                val[m+1][n+1][o+1] =  
                  get_volume_real_value(data, 
                                        index[data_xyzv[VIO_X]]+m, 
                                        index[data_xyzv[VIO_Y]]+n, 
                                        index[data_xyzv[VIO_Z]]+o, 0,0);
          Lvv = return_Lvv(val, (VIO_Real)VERY_SMALL_EPS);
        }
       
        if (max_val < Lvv) max_val = Lvv;
        if (min_val > Lvv) min_val = Lvv;
         
        float_vol  [index[data_xyzv[VIO_X]]][index[data_xyzv[VIO_Y]]][index[data_xyzv[VIO_Z]]]=Lvv;
        
        count++;
        update_progress_report( &progress, count );
        
  }
  terminate_progress_report(&progress);
 
  min_val *= 0.9;
  max_val *= 0.9;  
  

  
  set_volume_real_range(lvv  , min_val, max_val);



  for(i=0; i<sizes[0]; i++)
    for(j=0; j<sizes[1]; j++)
      for(k=0; k<sizes[2]; k++) 
       {
        Lvv = float_vol[i][j][k];
        if (Lvv < min_val) Lvv = min_val;
        if (Lvv > max_val) Lvv = max_val;
          

        set_volume_real_value(lvv  , i,j,k, 0, 0, Lvv);
      }

  VIO_FREE3D(float_vol);

  print ("Saving data (%f,%f)...\n",max_val, min_val);
  

  sprintf(output_filename,"%s_Lvv.mnc",argv[2]);
  stat = output_modified_volume(output_filename, NC_UNSPECIFIED, FALSE, 
                                0.0, 0.0,  lvv, argv[1], history, NULL);
  if (stat != VIO_OK) {
    print ("Error: cannot write Lvv %s.\n",output_filename);
    exit(EXIT_FAILURE);
  }

  print ("done.\n");

  exit(EXIT_SUCCESS);
}



static void init_the_volume_to_zero(VIO_Volume volume)
{
    int             v0, v1, v2, v3, v4;
    VIO_Real            zero;
  
    zero = CONVERT_VALUE_TO_VOXEL(volume, 0.0);

    BEGIN_ALL_VOXELS( volume, v0, v1, v2, v3, v4 )
      
      set_volume_voxel_value( volume, v0, v1, v2, v3, v4, zero );
    
    END_ALL_VOXELS

}

static void get_volume_XYZV_indices(VIO_Volume data, int xyzv[])
{
  
  int 
    axis, i, vol_dims;
  char 
    **data_dim_names;
  
  vol_dims       = get_volume_n_dimensions(data);
  data_dim_names = get_volume_dimension_names(data);
  
  for(i=0; i<VIO_N_DIMENSIONS+1; i++) xyzv[i] = -1;
  for(i=0; i<vol_dims; i++) {
    if (convert_dim_name_to_spatial_axis(data_dim_names[i], &axis )) {
      xyzv[axis] = i; 
    } 
    else {     /* not a spatial axis */
      xyzv[VIO_Z+1] = i;
    }
  }
  delete_dimension_names(data,data_dim_names);
  
}


static void estimate_3D_derivatives(VIO_Real r[3][3][3], 
                                    deriv_3D_struct *d) 

{

  VIO_Real        *p11, *p12, *p13;
  VIO_Real        *p21, *p22, *p23;
  VIO_Real        *p31, *p32, *p33;

  VIO_Real        slice_u1, slice_u2, slice_u3;
  VIO_Real        slice_v1, slice_v2, slice_v3;
  VIO_Real        slice_w1, slice_w2, slice_w3;
  
  VIO_Real        edge_u1_v1, /* edge_u1_v2,*/ edge_u1_v3;
/*  VIO_Real        edge_u2_v1, edge_u2_v2, edge_u2_v3; */
  VIO_Real        edge_u3_v1,  /* edge_u3_v2,*/ edge_u3_v3;
  VIO_Real        edge_u1_w1, edge_u1_w2, edge_u1_w3;
  VIO_Real        edge_u2_w1, edge_u2_w2, edge_u2_w3;
  VIO_Real        edge_u3_w1, edge_u3_w2, edge_u3_w3;
  VIO_Real        edge_v1_w1, edge_v1_w2, edge_v1_w3;
  VIO_Real        edge_v2_w1, edge_v2_w2, edge_v2_w3;
  VIO_Real        edge_v3_w1, edge_v3_w2, edge_v3_w3;
  
  /* --- 3x3x3 [u][v][w] --- */
  
  p11 = r[0][0]; 
  p12 = r[0][1]; 
  p13 = r[0][2]; 
  p21 = r[1][0]; 
  p22 = r[1][1]; 
  p23 = r[1][2]; 
  p31 = r[2][0]; 
  p32 = r[2][1]; 
  p33 = r[2][2]; 
  
                                /* lines varying along w */
  edge_u1_v1 = ( *p11     + *(p11+1) + *(p11+2));

/*  edge_u1_v2 = ( *p12     + *(p12+1) + *(p12+2)); */
  edge_u1_v3 = ( *p13     + *(p13+1) + *(p13+2));
/*  edge_u2_v1 = ( *p21     + *(p21+1) + *(p21+2)); */
/*  edge_u2_v2 = ( *p22     + *(p22+1) + *(p22+2)); */
/*  edge_u2_v3 = ( *p23     + *(p23+1) + *(p23+2)); */
  edge_u3_v1 = ( *p31     + *(p31+1) + *(p31+2));
/*  edge_u3_v2 = ( *p32     + *(p32+1) + *(p32+2)); */
  edge_u3_v3 = ( *p33     + *(p33+1) + *(p33+2));
  
                                /* lines varying along v */
  edge_u1_w1 = (  *p11    +  *p12    +  *p13   );
  edge_u1_w2 = ( *(p11+1) + *(p12+1) + *(p13+1));
  edge_u1_w3 = ( *(p11+2) + *(p12+2) + *(p13+2));
  edge_u2_w1 = (  *p21    +  *p22    +  *p23   );
  edge_u2_w2 = ( *(p21+1) + *(p22+1) + *(p23+1));
  edge_u2_w3 = ( *(p21+2) + *(p22+2) + *(p23+2));
  edge_u3_w1 = (  *p31    +  *p32    +  *p33   );
  edge_u3_w2 = ( *(p31+1) + *(p32+1) + *(p33+1)); 
  edge_u3_w3 = ( *(p31+2) + *(p32+2) + *(p33+2));
  
                                /* lines varying along u */
  edge_v1_w1 = (  *p11    +  *p21    +  *p31   );
  edge_v1_w2 = ( *(p11+1) + *(p21+1) + *(p31+1));
  edge_v1_w3 = ( *(p11+2) + *(p21+2) + *(p31+2));
  edge_v2_w1 = (  *p12    +  *p22    +  *p32   );
  edge_v2_w2 = ( *(p12+1) + *(p22+1) + *(p32+1));
  edge_v2_w3 = ( *(p12+2) + *(p22+2) + *(p32+2));
  edge_v3_w1 = (  *p13    +  *p23    +  *p33   );
  edge_v3_w2 = ( *(p13+1) + *(p23+1) + *(p33+1));
  edge_v3_w3 = ( *(p13+2) + *(p23+2) + *(p33+2));
  
  slice_u1 =  (edge_u1_w1 + edge_u1_w2 + edge_u1_w3);
  slice_u2 =  (edge_u2_w1 + edge_u2_w2 + edge_u2_w3);
  slice_u3 =  (edge_u3_w1 + edge_u3_w2 + edge_u3_w3);
  slice_v1 =  (edge_v1_w1 + edge_v1_w2 + edge_v1_w3);
  slice_v2 =  (edge_v2_w1 + edge_v2_w2 + edge_v2_w3);
  slice_v3 =  (edge_v3_w1 + edge_v3_w2 + edge_v3_w3);
  slice_w1 =  (edge_u1_w1 + edge_u2_w1 + edge_u3_w1);
  slice_w2 =  (edge_u1_w2 + edge_u2_w2 + edge_u3_w2);
  slice_w3 =  (edge_u1_w3 + edge_u2_w3 + edge_u3_w3);
  
  d->u  = (slice_u3 - slice_u1) / 18.0;                          
  d->v  = (slice_v3 - slice_v1) / 18.0;                          
  d->w  = (slice_w3 - slice_w1) / 18.0;                          
  d->uu = (slice_u3 + slice_u1 - 2*slice_u2) / 9.0;                   
  d->vv = (slice_v3 + slice_v1 - 2*slice_v2) / 9.0;                   
  d->ww = (slice_w3 + slice_w1 - 2*slice_w2) / 9.0;                   
  d->uv = (edge_u3_v3 + edge_u1_v1 - edge_u3_v1 - edge_u1_v3) / 12.0; 
  d->uw = (edge_u3_w3 + edge_u1_w1 - edge_u3_w1 - edge_u1_w3) / 12.0;  
  d->vw = (edge_v3_w3 + edge_v1_w1 - edge_v3_w1 - edge_v1_w3) / 12.0;  
  
} /* estimate_3D_derivatives */



static VIO_Real return_Lvv(VIO_Real r[3][3][3],
                       VIO_Real eps)
     
{
  deriv_3D_struct 
    d;                        

  VIO_Real
    S,                                /* mean curvature          */
    Lvv,
    sq_mag_grad,                /* square of magnitude of gradient   */
    x,y,z,                        /* first order derivatives */
    xx,yy,zz,xy,xz,yz;                /* second order derivative */


  d.u  = (r[2][1][1] - r[0][1][1] ) / 2.0;
  d.v  = (r[1][2][1] - r[1][0][1] ) / 2.0;
  d.w  = (r[1][1][2] - r[1][1][0] ) / 2.0;
  d.uu = (r[2][1][1] + r[0][1][1] - 2*r[1][1][1]);
  d.vv = (r[1][2][1] + r[1][0][1] - 2*r[1][1][1]);
  d.ww = (r[1][1][2] + r[1][1][0] - 2*r[1][1][1]);
  d.uv = (r[2][2][1] + r[0][0][1] - r[0][2][1] - r[2][0][1]) / 4.0;
  d.uw = (r[2][1][2] + r[0][1][0] - r[0][1][2] - r[2][1][0]) / 4.0;
  d.vw = (r[1][2][2] + r[1][0][0] - r[1][0][2] - r[1][2][0]) / 4.0;

  x  = d.u;  y  = d.v;  z  = d.w;
  xx = d.uu; yy = d.vv; zz = d.ww;
  xy = d.uv; yz = d.vw; xz = d.uw;

  Lvv = 0.0;
  sq_mag_grad = x*x + y*y + z*z;

  if ( fabs(sq_mag_grad) > eps )  {
                                /* Mean curvature */
    S = (
         x*x*(yy + zz) - 2*y*z*yz +
         y*y*(xx + zz) - 2*x*z*xz +
         z*z*(xx + yy) - 2*x*y*xy
         )
          / (2 * sqrt(sq_mag_grad*sq_mag_grad*sq_mag_grad));

    Lvv = sq_mag_grad * S;
  }

  return(Lvv);
}
