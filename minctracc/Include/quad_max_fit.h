/*------------------------------ MNI Header ----------------------------------
@NAME       : quad_max_fit.h
@DESCRIPTION: prototypes for Numerical/quad_max_fit.c
@CREATED    : Mon Nov  3, 1997 , Louis Collins
@MODIFIED   : not yet!
@VERSION    : $Id: quad_max_fit.h,v 1.3 2004-02-12 05:54:16 rotor Exp $
-----------------------------------------------------------------------------*/
    /* local structures */

typedef struct {
  Real 
    u,v,w,
    uu,vv,ww,
    uv,uw,vw;
} deriv_3D_struct;


typedef struct {
  Real 
    u,v,
    uu,vv,
    uv;
} deriv_2D_struct;


  void    estimate_3D_derivatives(Real r[3][3][3], 
					deriv_3D_struct *c);	     
  void    estimate_3D_derivatives_new(Real r[3][3][3], 
					     deriv_3D_struct *c);	     
  void    estimate_3D_derivatives_weighted(Real r[3][3][3], 
						  deriv_3D_struct *c);	     
BOOLEAN return_2D_disp_from_quad_fit(Real r[3][3], 
					    Real *dispu, 
					    Real *dispv);

BOOLEAN eigen(double **inputMat, 
		     int    ndim, 
		     double *eigen_val, 
		     double **eigen_vec, 
		     int    *iters);

BOOLEAN return_3D_disp_from_quad_fit(Real r[3][3][3], 
					    Real *dispu, 
					    Real *dispv, 
					    Real *dispw);

BOOLEAN return_3D_disp_from_min_quad_fit(Real r[3][3][3], 
						Real *dispu, 
						Real *dispv, 
						Real *dispw);

BOOLEAN return_2D_disp_from_quad_fit(Real r[3][3], 
					    Real *dispu, 
					    Real *dispv);

void estimate_3D_derivatives(Real r[3][3][3], 
				    deriv_3D_struct *d);

void estimate_3D_derivatives_weighted(Real r[3][3][3], 
					     deriv_3D_struct *d);

void estimate_3D_derivatives_new(Real r[3][3][3], 
					deriv_3D_struct *d);

BOOLEAN return_principal_directions(Real r[3][3][3],
					   Real dir_1[3],
					   Real dir_2[3],
					   Real *r_K,
					   Real *r_S,
					   Real *r_k1,
					   Real *r_k2,
					   Real *r_norm,
					   Real *r_Lvv,
					   Real eps);

BOOLEAN return_2D_principal_directions(Real r[3][3],
					      Real norm[3],
					      Real tang[3],
					      Real *K,
					      Real eps);

Real return_Lvv(Real r[3][3][3],
		       Real eps);

BOOLEAN return_local_eigen(Real r[3][3][3],
				  Real dir_1[3],
				  Real dir_2[3],
				  Real dir_3[3],
				  Real val[3]);

BOOLEAN return_local_eigen_from_hessian(Real r[3][3][3],
					       Real dir_1[3],
					       Real dir_2[3],
					       Real dir_3[3],
					       Real val[3]);
