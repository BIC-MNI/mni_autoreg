/* ----------------------------- MNI Header -----------------------------------
@NAME       : test_deriv
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : status
@DESCRIPTION: Program to deriv fitting routines.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Mon May  8 11:29:40 MET DST 1995  LC
@MODIFIED   : $Log: test_deriv.c,v $
@MODIFIED   : Revision 1.1  1999-10-25 19:52:11  louis
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
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Extra_progs/test_deriv.c,v 1.1 1999-10-25 19:52:11 louis Exp $";
#endif

#include <stdio.h>
#include <internal_volume_io.h>

/* Constants */
#ifndef TRUE
#  define TRUE 1
#  define FALSE 0
#endif
#ifndef public
#  define public
#  define private static
#endif


typedef struct {
  Real 
    u,v,w,
    uu,vv,ww,
    uv,uw,vw;
} deriv_3D_struct;

static char *default_dim_names[N_DIMENSIONS] =
   { MIzspace, MIyspace, MIxspace };

public  void    estimate_3D_derivatives(Real r[3][3][3], 
					deriv_3D_struct *c);	     


int main(int argc, char *argv[])
{
  float tmp[3][3][3];
  Real local_corr[3][3][3];
  Real 
    r[3][3][3];
  int 
    sizes[MAX_DIMENSIONS],
    i,j,k,
    u,v,w,flag;
  float test_val;
  Status status;

  deriv_3D_struct d;
  Volume
    data,
    v_du, v_dv, v_dw,
    v_uu, v_vv, v_ww,
    v_uv, v_uw, v_vw;

  if (argc!=2) {
    print ("usage:  test_deriv orig.mnc \n");
    exit(EXIT_FAILURE);
  }

  if (input_volume(argv[1], 3, default_dim_names, NC_UNSPECIFIED, FALSE, 0.0,0.0,
		   TRUE, &data, (minc_input_options *)NULL) != OK) { 
    print ("Error reading %s.\n",argv[1]);
    exit(EXIT_FAILURE);
  }

  v_du = copy_volume_definition(data, NC_UNSPECIFIED, FALSE, 0.0,0.0);
  set_volume_real_range(v_du, -5.0, 5.0);
  v_dv = copy_volume_definition(v_du, NC_UNSPECIFIED, FALSE, 0.0,0.0);
  v_dw = copy_volume_definition(v_du, NC_UNSPECIFIED, FALSE, 0.0,0.0);
  v_uu = copy_volume_definition(v_du, NC_UNSPECIFIED, FALSE, 0.0,0.0);
  v_vv = copy_volume_definition(v_du, NC_UNSPECIFIED, FALSE, 0.0,0.0);
  v_ww = copy_volume_definition(v_du, NC_UNSPECIFIED, FALSE, 0.0,0.0);
  v_uv = copy_volume_definition(v_du, NC_UNSPECIFIED, FALSE, 0.0,0.0);
  v_uw = copy_volume_definition(v_du, NC_UNSPECIFIED, FALSE, 0.0,0.0);
  v_vw = copy_volume_definition(v_du, NC_UNSPECIFIED, FALSE, 0.0,0.0);
  v_du = copy_volume_definition(v_du, NC_UNSPECIFIED, FALSE, 0.0,0.0);

  get_volume_sizes(data, sizes);

  for_less(i, 1, sizes[0]-1) {
    print ("slice %3d:%3d\n",i,sizes[0]-1);
    for_less(j, 1, sizes[1]-1) 
      for_less(k, 1, sizes[2]-1) {
	
	
	for_inclusive(u, -1, 1)
	  for_inclusive(v, -1, 1)
	    for_inclusive(w, -1, 1) {
	      GET_VALUE_3D( r[u+1][v+1][w+1], data, i+w, j+v, k+u );
	    }
	     
	estimate_3D_derivatives(r, &d);
	      
	set_volume_real_value(v_du, i, j, k, 0, 0, d.u  );
	set_volume_real_value(v_dv, i, j, k, 0, 0, d.v  );
	set_volume_real_value(v_dw, i, j, k, 0, 0, d.w  );
	set_volume_real_value(v_uu, i, j, k, 0, 0, d.uu );
	set_volume_real_value(v_vv, i, j, k, 0, 0, d.vv );
	set_volume_real_value(v_ww, i, j, k, 0, 0, d.ww );
	set_volume_real_value(v_uv, i, j, k, 0, 0, d.uv );
	set_volume_real_value(v_uw, i, j, k, 0, 0, d.uw );
	set_volume_real_value(v_vw, i, j, k, 0, 0, d.vw );
	
      }
  }  

  status = output_volume("v_du.mnc", NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			 v_du, NULL, (minc_output_options *)NULL);
  if (status == OK)
    status = output_volume("v_dv.mnc", NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			   v_dv, NULL, (minc_output_options *)NULL);
  if (status == OK)
    status = output_volume("v_dw.mnc", NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			   v_dw, NULL, (minc_output_options *)NULL);
  if (status == OK)
    status = output_volume("v_uu.mnc", NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			   v_uu, NULL, (minc_output_options *)NULL);
  if (status == OK)
    status = output_volume("v_vv.mnc", NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			   v_vv, NULL, (minc_output_options *)NULL);
  if (status == OK)
    status = output_volume("v_ww.mnc", NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			   v_ww, NULL, (minc_output_options *)NULL);
  if (status == OK)
    status = output_volume("v_uv.mnc", NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			   v_uv, NULL, (minc_output_options *)NULL);
  if (status == OK)
    status = output_volume("v_uw.mnc", NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			   v_uw, NULL, (minc_output_options *)NULL);
  if (status == OK)
    status = output_volume("v_vw.mnc", NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			   v_vw, NULL, (minc_output_options *)NULL);
  if (status != OK)
    print("error in output\n");

  
}
