/*------------------------------ MNI Header ----------------------------------
@NAME       : quad_max_fit.h
@DESCRIPTION: prototypes for Numerical/quad_max_fit.c
@CREATED    : Mon Nov  3, 1997 , Louis Collins
@MODIFIED   : not yet!
@VERSION    : $Id: quad_max_fit.h,v 1.4 2006-11-29 09:09:32 rotor Exp $
-----------------------------------------------------------------------------*/
    /* local structures */

typedef struct {
  VIO_Real 
    u,v,w,
    uu,vv,ww,
    uv,uw,vw;
} deriv_3D_struct;


typedef struct {
  VIO_Real 
    u,v,
    uu,vv,
    uv;
} deriv_2D_struct;


  void    estimate_3D_derivatives(VIO_Real r[3][3][3], 
                                        deriv_3D_struct *c);             
  void    estimate_3D_derivatives_new(VIO_Real r[3][3][3], 
                                             deriv_3D_struct *c);             
  void    estimate_3D_derivatives_weighted(VIO_Real r[3][3][3], 
                                                  deriv_3D_struct *c);             
VIO_BOOL return_2D_disp_from_quad_fit(VIO_Real r[3][3], 
                                            VIO_Real *dispu, 
                                            VIO_Real *dispv);

VIO_BOOL eigen(double **inputMat, 
                     int    ndim, 
                     double *eigen_val, 
                     double **eigen_vec, 
                     int    *iters);

VIO_BOOL return_3D_disp_from_quad_fit(VIO_Real r[3][3][3], 
                                            VIO_Real *dispu, 
                                            VIO_Real *dispv, 
                                            VIO_Real *dispw);

VIO_BOOL return_3D_disp_from_min_quad_fit(VIO_Real r[3][3][3], 
                                                VIO_Real *dispu, 
                                                VIO_Real *dispv, 
                                                VIO_Real *dispw);

VIO_BOOL return_2D_disp_from_quad_fit(VIO_Real r[3][3], 
                                            VIO_Real *dispu, 
                                            VIO_Real *dispv);

void estimate_3D_derivatives(VIO_Real r[3][3][3], 
                                    deriv_3D_struct *d);

void estimate_3D_derivatives_weighted(VIO_Real r[3][3][3], 
                                             deriv_3D_struct *d);

void estimate_3D_derivatives_new(VIO_Real r[3][3][3], 
                                        deriv_3D_struct *d);

VIO_BOOL return_principal_directions(VIO_Real r[3][3][3],
                                           VIO_Real dir_1[3],
                                           VIO_Real dir_2[3],
                                           VIO_Real *r_K,
                                           VIO_Real *r_S,
                                           VIO_Real *r_k1,
                                           VIO_Real *r_k2,
                                           VIO_Real *r_norm,
                                           VIO_Real *r_Lvv,
                                           VIO_Real eps);

VIO_BOOL return_2D_principal_directions(VIO_Real r[3][3],
                                              VIO_Real norm[3],
                                              VIO_Real tang[3],
                                              VIO_Real *K,
                                              VIO_Real eps);

VIO_Real return_Lvv(VIO_Real r[3][3][3],
                       VIO_Real eps);

VIO_BOOL return_local_eigen(VIO_Real r[3][3][3],
                                  VIO_Real dir_1[3],
                                  VIO_Real dir_2[3],
                                  VIO_Real dir_3[3],
                                  VIO_Real val[3]);

VIO_BOOL return_local_eigen_from_hessian(VIO_Real r[3][3][3],
                                               VIO_Real dir_1[3],
                                               VIO_Real dir_2[3],
                                               VIO_Real dir_3[3],
                                               VIO_Real val[3]);
