
#ifndef  DEF_SPLINES
#define  DEF_SPLINES

#include  <volume_io.h>

#define  LINEAR_COEF_00     1.0
#define  LINEAR_COEF_01     0.0
#define  LINEAR_COEF_10    -1.0
#define  LINEAR_COEF_11     1.0

#define  LINEAR_COEF_0( v0, v1 ) \
            (v0)
#define  LINEAR_COEF_1( v0, v1 ) \
            ((v1) - (v0))

#define LINEAR_UNIVAR( v0, v1, u ) \
            ( (v0) + (u) * ((v1) - (v0)) )

/* ------------------ QUADRATIC INTERPOLATING SPLINES ----------------------- */

#define  QUADRATIC_COEF_00    0.5
#define  QUADRATIC_COEF_01    0.5
#define  QUADRATIC_COEF_02    0.0
#define  QUADRATIC_COEF_10  (-1.0)
#define  QUADRATIC_COEF_11    1.0
#define  QUADRATIC_COEF_12    0.0
#define  QUADRATIC_COEF_20    0.5
#define  QUADRATIC_COEF_21  (-1.0)
#define  QUADRATIC_COEF_22    0.5

#define  QUADRATIC_COEF_OF_V0( u ) \
          ( QUADRATIC_COEF_00 + (u) * ( \
            QUADRATIC_COEF_10 + (u) * \
            QUADRATIC_COEF_20) )
#define  QUADRATIC_COEF_OF_V1( u ) \
          ( QUADRATIC_COEF_01 + (u) * ( \
            QUADRATIC_COEF_11 + (u) * \
            QUADRATIC_COEF_21) )
#define  QUADRATIC_COEF_OF_V2( u ) \
          ( QUADRATIC_COEF_02 + (u) * ( \
            QUADRATIC_COEF_12 + (u) * \
            QUADRATIC_COEF_22) )

#define  COMPUTE_QUADRATIC_COEFFS( u, c0, c1, c2 ) \
          { \
              (c0) = QUADRATIC_COEF_OF_V0( u ); \
              (c1) = QUADRATIC_COEF_OF_V1( u ); \
              (c2) = QUADRATIC_COEF_OF_V2( u ); \
          }

#define  QUADRATIC_UNIVAR( v0, v1, v2, u )                               \
     ( (v0) * QUADRATIC_COEF_OF_V0(u) + \
       (v1) * QUADRATIC_COEF_OF_V1(u) + \
       (v2) * QUADRATIC_COEF_OF_V2(u) )

#define  QUADRATIC_DERIV_COEF_OF_V0( u ) \
                ( QUADRATIC_COEF_10 + 2.0 * (u) * \
                  QUADRATIC_COEF_20 )
#define  QUADRATIC_DERIV_COEF_OF_V1( u ) \
                ( QUADRATIC_COEF_11 + 2.0 * (u) * \
                  QUADRATIC_COEF_21 )
#define  QUADRATIC_DERIV_COEF_OF_V2( u ) \
                ( QUADRATIC_COEF_12 + 2.0 * (u) * \
                  QUADRATIC_COEF_22 )

#define  COMPUTE_QUADRATIC_DERIV_COEFFS( u, c0, c1, c2 ) \
          { \
              (c0) = QUADRATIC_DERIV_COEF_OF_V0( u ); \
              (c1) = QUADRATIC_DERIV_COEF_OF_V1( u ); \
              (c2) = QUADRATIC_DERIV_COEF_OF_V2( u ); \
          }

#define  QUADRATIC_UNIVAR_DERIV( v0, v1, v2, u )                              \
     ( (v0) * QUADRATIC_DERIV_COEF_OF_V0(u) + \
       (v1) * QUADRATIC_DERIV_COEF_OF_V1(u) + \
       (v2) * QUADRATIC_DERIV_COEF_OF_V2(u) )

#define  QUADRATIC_DERIV2_COEF_OF_V0( u ) \
           ( 2.0 * QUADRATIC_COEF_20 )
#define  QUADRATIC_DERIV2_COEF_OF_V1( u ) \
           ( 2.0 * QUADRATIC_COEF_21 )
#define  QUADRATIC_DERIV2_COEF_OF_V2( u ) \
           ( 2.0 * QUADRATIC_COEF_22 )

#define  COMPUTE_QUADRATIC_DERIV2_COEFFS( u, c0, c1, c2 ) \
          { \
              (c0) = QUADRATIC_DERIV2_COEF_OF_V0( u ); \
              (c1) = QUADRATIC_DERIV2_COEF_OF_V1( u ); \
              (c2) = QUADRATIC_DERIV2_COEF_OF_V2( u ); \
          }

#define  QUADRATIC_UNIVAR_DERIV2( v0, v1, v2, u )                             \
     ( (v0) * QUADRATIC_DERIV2_COEF_OF_V0(u) + \
       (v1) * QUADRATIC_DERIV2_COEF_OF_V1(u) + \
       (v2) * QUADRATIC_DERIV2_COEF_OF_V2(u) )

#define  DOT3( a1, b1, c1, a2, b2, c2 ) \
         ( (a1) * (a2) + (b1) * (b2) + (c1) * (c2) )

#define  MULT3( r1, r2, r3, v1, v2, v3, v00, v01, v02, v10, v11, v12, v20, v21, v22 ) \
         { (r1) = DOT3( v1, v2, v3, v00, v01, v02 ); \
           (r2) = DOT3( v1, v2, v3, v10, v11, v12 ); \
           (r3) = DOT3( v1, v2, v3, v20, v21, v22 ); \
         }

#define  QUADRATIC_BIVAR( v00, v01, v02, v10, v11, v12, v20, v21, v22, u_parm, v_parm, val ) \
     { \
         Real  wv0, wv1, wv2; \
         Real  v0, v1, v2; \
 \
         COMPUTE_QUADRATIC_COEFFS( v_parm, wv0, wv1, wv2 ); \
 \
         MULT3( v0, v1, v2, wv0, wv1, wv2, \
                 v00, v01, v02, v10, v11, v12, v20, v21, v22 ); \
 \
         (val) = QUADRATIC_UNIVAR( v0, v1, v2, u_parm ); \
     }

#define  QUADRATIC_BIVAR_DERIV( v00, v01, v02, v10, v11, v12, v20, v21, v22, u, v, val, du, dv ) \
     { \
         Real  wu0, wu1, wu2; \
         Real  wv0, wv1, wv2; \
         Real  wdv0, wdv1, wdv2; \
         Real  v0, v1, v2; \
         Real  dv0, dv1, dv2; \
 \
         COMPUTE_QUADRATIC_COEFFS( u_parm, wu0, wu1, wu2 ); \
         COMPUTE_QUADRATIC_COEFFS( v_parm, wv0, wv1, wv2 ); \
         COMPUTE_QUADRATIC_DERIV_COEFFS( v_parm, wdv0, wdv1, wdv2 ); \
 \
         MULT3( v0, v1, v2, wv0, wv1, wv2, \
                 v00, v01, v02, v10, v11, v12, v20, v21, v22 ); \
\
         MULT3( dv0, dv1, dv2, wdv0, wdv1, wdv2, \
                 v00, v01, v02, v10, v11, v12, v20, v21, v22 ); \
 \
         (dv) = DOT3( wu0, wu1, wu2, dv0, dv1, dv2 ); \
         (du) = QUADRATIC_UNIVAR_DERIV( v0, v1, v2, u_parm ); \
         (val) = DOT3( wu0, wu1, wu2, v0, v1, v2 ); \
     }

#define  QUADRATIC_BIVAR_DERIV2( v00, v01, v02, v10, v11, v12, v20, v21, v22, u_parm, v_parm, val, du, dv, duu, duv, dvv ) \
     { \
         Real  wu0, wu1, wu2; \
         Real  wdu0, wdu1, wdu2; \
         Real  wv0, wv1, wv2; \
         Real  wdv0, wdv1, wdv2; \
         Real  wdvv0, wdvv1, wdvv2; \
         Real  v0, v1, v2; \
         Real  dv0, dv1, dv2; \
         Real  dvv0, dvv1, dvv2; \
 \
         COMPUTE_QUADRATIC_COEFFS( u_parm, wu0, wu1, wu2 ); \
         COMPUTE_QUADRATIC_DERIV_COEFFS( u_parm, wdu0, wdu1, wdu2 ); \
         COMPUTE_QUADRATIC_COEFFS( v_parm, wv0, wv1, wv2 ); \
         COMPUTE_QUADRATIC_DERIV_COEFFS( v_parm, wdv0, wdv1, wdv2 ); \
         COMPUTE_QUADRATIC_DERIV2_COEFFS( v_parm, wdvv0, wdvv1, wdvv2 ); \
 \
         MULT3( v0, v1, v2, wv0, wv1, wv2, \
                 v00, v01, v02, v10, v11, v12, v20, v21, v22 ); \
         MULT3( dv0, dv1, dv2, wdv0, wdv1, wdv2, \
                 v00, v01, v02, v10, v11, v12, v20, v21, v22 ); \
         MULT3( dvv0, dvv1, dvv2, wdvv0, wdvv1, wdvv2, \
                 v00, v01, v02, v10, v11, v12, v20, v21, v22 ); \
 \
         (val) = DOT3( wu0, wu1, wu2,   v0,   v1,   v2 ); \
         (dv)  = DOT3( wu0, wu1, wu2,  dv0,  dv1,  dv2 ); \
         (dvv) = DOT3( wu0, wu1, wu2, dvv0, dvv1, dvv2 ); \
         (du)  = DOT3( wdu0, wdu1, wdu2, v0, v1, v2 ); \
         (duv) = DOT3( wdu0, wdu1, wdu2, dv0, dv1, dv2 ); \
         (duu) = QUADRATIC_UNIVAR_DERIV2( v0, v1, v2, u_parm ); \
     }

#define  QUADRATIC_TRIVAR( c, u_parm, v_parm, w_parm, val ) \
    { \
         Real  wv0, wv1, wv2; \
         Real  ww0, ww1, ww2; \
         Real  v00, v01, v02, v10, v11, v12, v20, v21, v22; \
         Real  v0, v1, v2; \
 \
         COMPUTE_QUADRATIC_COEFFS( v_parm, wv0, wv1, wv2 ); \
         COMPUTE_QUADRATIC_COEFFS( w_parm, ww0, ww1, ww2 ); \
 \
         MULT3( v00, v01, v02, ww0, ww1, ww2, \
                 GLUE(c,000), GLUE(c,001), GLUE(c,002), \
                 GLUE(c,010), GLUE(c,011), GLUE(c,012), \
                 GLUE(c,020), GLUE(c,021), GLUE(c,022) ); \
         MULT3( v10, v11, v12, ww0, ww1, ww2, \
                 GLUE(c,100), GLUE(c,101), GLUE(c,102), \
                 GLUE(c,110), GLUE(c,111), GLUE(c,112), \
                 GLUE(c,120), GLUE(c,121), GLUE(c,122) ); \
         MULT3( v20, v21, v22, ww0, ww1, ww2, \
                 GLUE(c,200), GLUE(c,201), GLUE(c,202), \
                 GLUE(c,210), GLUE(c,211), GLUE(c,212), \
                 GLUE(c,220), GLUE(c,221), GLUE(c,222) ); \
 \
         MULT3( v0, v1, v2, wv0, wv1, wv2, \
                 v00, v01, v02, v10, v11, v12, v20, v21, v22 ); \
 \
         (val) = QUADRATIC_UNIVAR( v0, v1, v2, u_parm ); \
    }

#define  QUADRATIC_TRIVAR_DERIV( c, u_parm, v_parm, w_parm, val, du, dv, dw ) \
    { \
         Real  wu0, wu1, wu2; \
         Real  wv0, wv1, wv2; \
         Real  wdv0, wdv1, wdv2; \
         Real  ww0, ww1, ww2; \
         Real  wdw0, wdw1, wdw2; \
         Real  v00, v01, v02, v10, v11, v12, v20, v21, v22; \
         Real  dw00, dw01, dw02, dw10, dw11, dw12, dw20, dw21, dw22; \
         Real  v0, v1, v2; \
         Real  dv0, dv1, dv2; \
         Real  dw0, dw1, dw2; \
 \
         COMPUTE_QUADRATIC_COEFFS( u_parm, wu0, wu1, wu2 ); \
         COMPUTE_QUADRATIC_COEFFS( v_parm, wv0, wv1, wv2 ); \
         COMPUTE_QUADRATIC_COEFFS( w_parm, ww0, ww1, ww2 ); \
         COMPUTE_QUADRATIC_DERIV_COEFFS( v_parm, wdv0, wdv1, wdv2 ); \
         COMPUTE_QUADRATIC_DERIV_COEFFS( w_parm, wdw0, wdw1, wdw2 ); \
 \
         MULT3( v00, v01, v02, ww0, ww1, ww2, \
                 GLUE(c,000), GLUE(c,001), GLUE(c,002), \
                 GLUE(c,010), GLUE(c,011), GLUE(c,012), \
                 GLUE(c,020), GLUE(c,021), GLUE(c,022) ); \
         MULT3( v10, v11, v12, ww0, ww1, ww2, \
                 GLUE(c,100), GLUE(c,101), GLUE(c,102), \
                 GLUE(c,110), GLUE(c,111), GLUE(c,112), \
                 GLUE(c,120), GLUE(c,121), GLUE(c,122) ); \
         MULT3( v20, v21, v22, ww0, ww1, ww2, \
                 GLUE(c,200), GLUE(c,201), GLUE(c,202), \
                 GLUE(c,210), GLUE(c,211), GLUE(c,212), \
                 GLUE(c,220), GLUE(c,221), GLUE(c,222) ); \
 \
         MULT3( dw00, dw01, dw02, wdw0, wdw1, wdw2, \
                 GLUE(c,000), GLUE(c,001), GLUE(c,002), \
                 GLUE(c,010), GLUE(c,011), GLUE(c,012), \
                 GLUE(c,020), GLUE(c,021), GLUE(c,022) ); \
         MULT3( dw10, dw11, dw12, wdw0, wdw1, wdw2, \
                 GLUE(c,100), GLUE(c,101), GLUE(c,102), \
                 GLUE(c,110), GLUE(c,111), GLUE(c,112), \
                 GLUE(c,120), GLUE(c,121), GLUE(c,122) ); \
         MULT3( dw20, dw21, dw22, wdw0, wdw1, wdw2, \
                 GLUE(c,200), GLUE(c,201), GLUE(c,202), \
                 GLUE(c,210), GLUE(c,211), GLUE(c,212), \
                 GLUE(c,220), GLUE(c,221), GLUE(c,222) ); \
 \
         MULT3( v0, v1, v2, wv0, wv1, wv2, \
                 v00, v01, v02, v10, v11, v12, v20, v21, v22 ); \
         MULT3( dv0, dv1, dv2, wdv0, wdv1, wdv2, \
                 v00, v01, v02, v10, v11, v12, v20, v21, v22 ); \
         MULT3( dw0, dw1, dw2, wv0, wv1, wv2, \
                 dw00, dw01, dw02, dw10, dw11, dw12, dw20, dw21, dw22 ); \
 \
         (val) = DOT3( v0, v1, v2, wu0, wu1, wu2 ); \
         (du) = QUADRATIC_UNIVAR_DERIV( v0, v1, v2, u_parm ); \
         (dv) = DOT3( dv0, dv1, dv2, wu0, wu1, wu2 ); \
         (dw) = DOT3( dw0, dw1, dw2, wu0, wu1, wu2 ); \
    }

#define  QUADRATIC_TRIVAR_DERIV2( c, u_parm, v_parm, w_parm, val, du, dv, dw, duu, duv, duw, dvv, dvw, dww ) \
    { \
         Real  wu0, wu1, wu2; \
         Real  wdu0, wdu1, wdu2; \
         Real  wv0, wv1, wv2; \
         Real  wdv0, wdv1, wdv2; \
         Real  wdvv0, wdvv1, wdvv2; \
         Real  ww0, ww1, ww2; \
         Real  wdw0, wdw1, wdw2; \
         Real  wdww0, wdww1, wdww2; \
         Real  v00, v01, v02, v10, v11, v12, v20, v21, v22; \
         Real  dw00, dw01, dw02, dw10, dw11, dw12, dw20, dw21, dw22; \
         Real  dww00, dww01, dww02, dww10, dww11, dww12, dww20, dww21,dww22;\
         Real  v0, v1, v2; \
         Real  dv0, dv1, dv2; \
         Real  dvv0, dvv1, dvv2; \
         Real  dw0, dw1, dw2; \
         Real  dww0, dww1, dww2; \
         Real  dvw0, dvw1, dvw2; \
 \
         COMPUTE_QUADRATIC_COEFFS( u_parm, wu0, wu1, wu2 ); \
         COMPUTE_QUADRATIC_COEFFS( v_parm, wv0, wv1, wv2 ); \
         COMPUTE_QUADRATIC_COEFFS( w_parm, ww0, ww1, ww2 ); \
         COMPUTE_QUADRATIC_DERIV_COEFFS( u_parm, wdu0, wdu1, wdu2 ); \
         COMPUTE_QUADRATIC_DERIV_COEFFS( v_parm, wdv0, wdv1, wdv2 ); \
         COMPUTE_QUADRATIC_DERIV_COEFFS( w_parm, wdw0, wdw1, wdw2 ); \
         COMPUTE_QUADRATIC_DERIV2_COEFFS( v_parm, wdvv0, wdvv1, wdvv2 ); \
         COMPUTE_QUADRATIC_DERIV2_COEFFS( w_parm, wdww0, wdww1, wdww2 ); \
 \
         MULT3( v00, v01, v02, ww0, ww1, ww2, \
                 GLUE(c,000), GLUE(c,001), GLUE(c,002), \
                 GLUE(c,010), GLUE(c,011), GLUE(c,012), \
                 GLUE(c,020), GLUE(c,021), GLUE(c,022) ); \
         MULT3( v10, v11, v12, ww0, ww1, ww2, \
                 GLUE(c,100), GLUE(c,101), GLUE(c,102), \
                 GLUE(c,110), GLUE(c,111), GLUE(c,112), \
                 GLUE(c,120), GLUE(c,121), GLUE(c,122) ); \
         MULT3( v20, v21, v22, ww0, ww1, ww2, \
                 GLUE(c,200), GLUE(c,201), GLUE(c,202), \
                 GLUE(c,210), GLUE(c,211), GLUE(c,212), \
                 GLUE(c,220), GLUE(c,221), GLUE(c,222) ); \
 \
         MULT3( dw00, dw01, dw02, wdw0, wdw1, wdw2, \
                 GLUE(c,000), GLUE(c,001), GLUE(c,002), \
                 GLUE(c,010), GLUE(c,011), GLUE(c,012), \
                 GLUE(c,020), GLUE(c,021), GLUE(c,022) ); \
         MULT3( dw10, dw11, dw12, wdw0, wdw1, wdw2, \
                 GLUE(c,100), GLUE(c,101), GLUE(c,102), \
                 GLUE(c,110), GLUE(c,111), GLUE(c,112), \
                 GLUE(c,120), GLUE(c,121), GLUE(c,122) ); \
         MULT3( dw20, dw21, dw22, wdw0, wdw1, wdw2, \
                 GLUE(c,200), GLUE(c,201), GLUE(c,202), \
                 GLUE(c,210), GLUE(c,211), GLUE(c,212), \
                 GLUE(c,220), GLUE(c,221), GLUE(c,222) ); \
 \
         MULT3( dww00, dww01, dww02, wdww0, wdww1, wdww2, \
                 GLUE(c,000), GLUE(c,001), GLUE(c,002), \
                 GLUE(c,010), GLUE(c,011), GLUE(c,012), \
                 GLUE(c,020), GLUE(c,021), GLUE(c,022) ); \
         MULT3( dww10, dww11, dww12, wdww0, wdww1, wdww2, \
                 GLUE(c,100), GLUE(c,101), GLUE(c,102), \
                 GLUE(c,110), GLUE(c,111), GLUE(c,112), \
                 GLUE(c,120), GLUE(c,121), GLUE(c,122) ); \
         MULT3( dww20, dww21, dww22, wdww0, wdww1, wdww2, \
                 GLUE(c,200), GLUE(c,201), GLUE(c,202), \
                 GLUE(c,210), GLUE(c,211), GLUE(c,212), \
                 GLUE(c,220), GLUE(c,221), GLUE(c,222) ); \
 \
         MULT3( v0, v1, v2, wv0, wv1, wv2, \
                v00, v01, v02, v10, v11, v12, v20, v21, v22 ); \
         MULT3( dv0, dv1, dv2, wdv0, wdv1, wdv2, \
                v00, v01, v02, v10, v11, v12, v20, v21, v22 ); \
         MULT3( dvv0, dvv1, dvv2, wdvv0, wdvv1, wdvv2, \
                v00, v01, v02, v10, v11, v12, v20, v21, v22 ); \
         MULT3( dw0, dw1, dw2, wv0, wv1, wv2, \
                dw00, dw01, dw02, dw10, dw11, dw12, dw20, dw21, dw22 ); \
         MULT3( dww0, dww1, dww2, wv0, wv1, wv2, \
                dww00, dww01, dww02, dww10, dww11, dww12, dww20, dww21, dww22);\
         MULT3( dvw0, dvw1, dvw2, wdv0, wdv1, wdv2, \
                dw00, dw01, dw02, dw10, dw11, dw12, dw20, dw21, dw22 ); \
 \
         (val) = DOT3( v0, v1, v2,       wu0, wu1, wu2 ); \
         (dv)  = DOT3( dv0, dv1, dv2,    wu0, wu1, wu2 ); \
         (dw)  = DOT3( dw0, dw1, dw2,    wu0, wu1, wu2 ); \
         (dvv) = DOT3( dvv0, dvv1, dvv2, wu0, wu1, wu2 ); \
         (dvw) = DOT3( dvw0, dvw1, dvw2, wu0, wu1, wu2 ); \
         (dww) = DOT3( dww0, dww1, dww2, wu0, wu1, wu2 ); \
         (du)  = DOT3( v0, v1, v2,       wdu0, wdu1, wdu2 ); \
         (duv) = DOT3( dv0, dv1, dv2,    wdu0, wdu1, wdu2 ); \
         (duw) = DOT3( dw0, dw1, dw2,    wdu0, wdu1, wdu2 ); \
         (duu) = QUADRATIC_UNIVAR_DERIV2( v0, v1, v2, u_parm ); \
    }

/* -------------------- CUBIC INTERPOLATING SPLINES ----------------------- */

#define  CUBIC_COEF_00    0.0
#define  CUBIC_COEF_01    1.0
#define  CUBIC_COEF_02    0.0
#define  CUBIC_COEF_03    0.0

#define  CUBIC_COEF_10  (-0.5)
#define  CUBIC_COEF_11    0.0
#define  CUBIC_COEF_12    0.5
#define  CUBIC_COEF_13    0.0

#define  CUBIC_COEF_20    1.0
#define  CUBIC_COEF_21  (-2.5)
#define  CUBIC_COEF_22    2.0
#define  CUBIC_COEF_23  (-0.5)

#define  CUBIC_COEF_30  (-0.5)
#define  CUBIC_COEF_31    1.5
#define  CUBIC_COEF_32  (-1.5)
#define  CUBIC_COEF_33    0.5

#define  CUBIC_COEF_OF_V0( u ) \
          ( CUBIC_COEF_00 + (u) * ( \
            CUBIC_COEF_10 + (u) * (\
            CUBIC_COEF_20 + (u) * \
            CUBIC_COEF_30)) )
#define  CUBIC_COEF_OF_V1( u ) \
          ( CUBIC_COEF_01 + (u) * ( \
            CUBIC_COEF_11 + (u) * (\
            CUBIC_COEF_21 + (u) * \
            CUBIC_COEF_31)) )
#define  CUBIC_COEF_OF_V2( u ) \
          ( CUBIC_COEF_02 + (u) * ( \
            CUBIC_COEF_12 + (u) * (\
            CUBIC_COEF_22 + (u) * \
            CUBIC_COEF_32)) )
#define  CUBIC_COEF_OF_V3( u ) \
          ( CUBIC_COEF_03 + (u) * ( \
            CUBIC_COEF_13 + (u) * (\
            CUBIC_COEF_23 + (u) * \
            CUBIC_COEF_33)) )

#define  COMPUTE_CUBIC_COEFFS( u, c0, c1, c2, c3 ) \
          { \
              (c0) = CUBIC_COEF_OF_V0( u ); \
              (c1) = CUBIC_COEF_OF_V1( u ); \
              (c2) = CUBIC_COEF_OF_V2( u ); \
              (c3) = CUBIC_COEF_OF_V3( u ); \
          }

#define  CUBIC_UNIVAR( v0, v1, v2, v3, u )                               \
     ( (v0) * CUBIC_COEF_OF_V0(u) + \
       (v1) * CUBIC_COEF_OF_V1(u) + \
       (v2) * CUBIC_COEF_OF_V2(u) + \
       (v3) * CUBIC_COEF_OF_V3(u) )

#define  CUBIC_DERIV_COEF_OF_V0( u ) \
         (        CUBIC_COEF_10 + (u) * \
           (2.0 * CUBIC_COEF_20 + (u) * \
            3.0 * CUBIC_COEF_30 ) )
#define  CUBIC_DERIV_COEF_OF_V1( u ) \
         (        CUBIC_COEF_11 + (u) * \
           (2.0 * CUBIC_COEF_21 + (u) * \
            3.0 * CUBIC_COEF_31 ) )
#define  CUBIC_DERIV_COEF_OF_V2( u ) \
         (        CUBIC_COEF_12 + (u) * \
           (2.0 * CUBIC_COEF_22 + (u) * \
            3.0 * CUBIC_COEF_32 ) )
#define  CUBIC_DERIV_COEF_OF_V3( u ) \
         (        CUBIC_COEF_13 + (u) * \
           (2.0 * CUBIC_COEF_23 + (u) * \
            3.0 * CUBIC_COEF_33 ) )

#define  COMPUTE_CUBIC_DERIV_COEFFS( u, c0, c1, c2, c3 ) \
          { \
              (c0) = CUBIC_DERIV_COEF_OF_V0( u ); \
              (c1) = CUBIC_DERIV_COEF_OF_V1( u ); \
              (c2) = CUBIC_DERIV_COEF_OF_V2( u ); \
              (c3) = CUBIC_DERIV_COEF_OF_V3( u ); \
          }

#define  CUBIC_UNIVAR_DERIV( v0, v1, v2, v3, u )                              \
     ( (v0) * CUBIC_DERIV_COEF_OF_V0(u) + \
       (v1) * CUBIC_DERIV_COEF_OF_V1(u) + \
       (v2) * CUBIC_DERIV_COEF_OF_V2(u) + \
       (v3) * CUBIC_DERIV_COEF_OF_V3(u) )

#define  CUBIC_DERIV2_COEF_OF_V0( u ) \
           ( 2.0 * CUBIC_COEF_20 + (u) * \
             6.0 * CUBIC_COEF_30 )
#define  CUBIC_DERIV2_COEF_OF_V1( u ) \
           ( 2.0 * CUBIC_COEF_21 + (u) * \
             6.0 * CUBIC_COEF_31 )
#define  CUBIC_DERIV2_COEF_OF_V2( u ) \
           ( 2.0 * CUBIC_COEF_22 + (u) * \
             6.0 * CUBIC_COEF_32 )
#define  CUBIC_DERIV2_COEF_OF_V3( u ) \
           ( 2.0 * CUBIC_COEF_23 + (u) * \
             6.0 * CUBIC_COEF_33 )

#define  COMPUTE_CUBIC_DERIV2_COEFFS( u, c0, c1, c2, c3 ) \
          { \
              (c0) = CUBIC_DERIV2_COEF_OF_V0( u ); \
              (c1) = CUBIC_DERIV2_COEF_OF_V1( u ); \
              (c2) = CUBIC_DERIV2_COEF_OF_V2( u ); \
              (c3) = CUBIC_DERIV2_COEF_OF_V3( u ); \
          }

#define  CUBIC_UNIVAR_DERIV2( v0, v1, v2, v3, u )                             \
     ( (v0) * CUBIC_DERIV2_COEF_OF_V0(u) + \
       (v1) * CUBIC_DERIV2_COEF_OF_V1(u) + \
       (v2) * CUBIC_DERIV2_COEF_OF_V2(u) + \
       (v3) * CUBIC_DERIV2_COEF_OF_V3(u) )

#define  DOT4( a1, b1, c1, d1, a2, b2, c2, d2 ) \
         ( (a1) * (a2) + (b1) * (b2) + (c1) * (c2) + (d1) * (d2) )

#define  MULT4( r1, r2, r3, r4, w1, w2, w3, w4, v ) \
       { (r1) = DOT4(w1,w2,w3,w4,GLUE(v,00),GLUE(v,01),GLUE(v,02),GLUE(v,03));\
         (r2) = DOT4(w1,w2,w3,w4,GLUE(v,10),GLUE(v,11),GLUE(v,12),GLUE(v,13));\
         (r3) = DOT4(w1,w2,w3,w4,GLUE(v,20),GLUE(v,21),GLUE(v,22),GLUE(v,23));\
         (r4) = DOT4(w1,w2,w3,w4,GLUE(v,30),GLUE(v,31),GLUE(v,32),GLUE(v,33));\
       }

#define  CUBIC_BIVAR(cv,u_parm,v_parm,val) \
     { \
         Real  wv0, wv1, wv2, wv3; \
         Real  v0, v1, v2, v3; \
 \
         COMPUTE_CUBIC_COEFFS( v_parm, wv0, wv1, wv2, wv3 ); \
 \
         MULT4( v0, v1, v2, v3, wv0, wv1, wv2, wv3,cv); \
 \
         (val) = CUBIC_UNIVAR( v0, v1, v2, v3, u_parm ); \
     }

#define  CUBIC_BIVAR_DERIV( cv, u_parm, v_parm, val, du, dv ) \
     { \
         Real  wu0, wu1, wu2, wu3; \
         Real  wv0, wv1, wv2, wv3; \
         Real  wdv0, wdv1, wdv2, wdv3; \
         Real  v0, v1, v2, v3; \
         Real  dv0, dv1, dv2, dv3; \
 \
         COMPUTE_CUBIC_COEFFS( u_parm, wu0, wu1, wu2, wu3 ); \
         COMPUTE_CUBIC_COEFFS( v_parm, wv0, wv1, wv2, wv3 ); \
         COMPUTE_CUBIC_DERIV_COEFFS( v_parm, wdv0, wdv1, wdv2, wdv3 ); \
 \
         MULT4(  v0,  v1,  v2,  v3,  wv0,  wv1,  wv2,  wv3,cv); \
         MULT4( dv0, dv1, dv2, dv3, wdv0, wdv1, wdv2, wdv3,cv); \
 \
         (dv) = DOT4( wu0, wu1, wu2, wu3, dv0, dv1, dv2, dv3 ); \
         (du) = CUBIC_UNIVAR_DERIV( v0, v1, v2, v3, u_parm ); \
         (val) = DOT4( wu0, wu1, wu2, wu3, v0, v1, v2, v3 ); \
     }

#define  CUBIC_BIVAR_DERIV2( cv, u_parm, v_parm, val, du, dv, duu, duv, dvv ) \
     { \
         Real  wu0, wu1, wu2, wu3; \
         Real  wdu0, wdu1, wdu2, wdu3; \
         Real  wv0, wv1, wv2, wv3; \
         Real  wdv0, wdv1, wdv2, wdv3; \
         Real  wdvv0, wdvv1, wdvv2, wdvv3; \
         Real  v0, v1, v2, v3; \
         Real  dv0, dv1, dv2, dv3; \
         Real  dvv0, dvv1, dvv2, dvv3; \
 \
         COMPUTE_CUBIC_COEFFS( u_parm, wu0, wu1, wu2, wu3 ); \
         COMPUTE_CUBIC_DERIV_COEFFS( u_parm, wdu0, wdu1, wdu2, wdu3 ); \
         COMPUTE_CUBIC_COEFFS( v_parm, wv0, wv1, wv2, wv3 ); \
         COMPUTE_CUBIC_DERIV_COEFFS( v_parm, wdv0, wdv1, wdv2, wdv3 ); \
         COMPUTE_CUBIC_DERIV2_COEFFS( v_parm, wdvv0, wdvv1, wdvv2, wdvv3 ); \
 \
         MULT4( v0, v1, v2, v3, wv0, wv1, wv2, wv3,cv); \
         MULT4( dv0, dv1, dv2, dv3, wdv0, wdv1, wdv2, wdv3,cv); \
         MULT4( dvv0, dvv1, dvv2, dvv3, wdvv0, wdvv1, wdvv2, wdvv3,cv); \
 \
         (val) = DOT4( wu0, wu1, wu2, wu3,   v0,   v1,   v2,   v3 ); \
         (dv)  = DOT4( wu0, wu1, wu2, wu3,  dv0,  dv1,  dv2,  dv3 ); \
         (dvv) = DOT4( wu0, wu1, wu2, wu3, dvv0, dvv1, dvv2, dvv3 ); \
         (du)  = DOT4( wdu0, wdu1, wdu2, wdu3, v0, v1, v2, v3 ); \
         (duv) = DOT4( wdu0, wdu1, wdu2, wdu3, dv0, dv1, dv2, dv3 ); \
         (duu) = CUBIC_UNIVAR_DERIV2( v0, v1, v2, v3, u_parm ); \
     }

#define  CUBIC_TRIVAR( c, u_parm, v_parm, w_parm, val ) \
    { \
         Real  wv0, wv1, wv2, wv3; \
         Real  ww0, ww1, ww2, ww3; \
         Real  v00, v01, v02, v03, v10, v11, v12, v13; \
         Real  v20, v21, v22, v23, v30, v31, v32, v33; \
         Real  v0, v1, v2, v3; \
 \
         COMPUTE_CUBIC_COEFFS( v_parm, wv0, wv1, wv2, wv3 ); \
         COMPUTE_CUBIC_COEFFS( w_parm, ww0, ww1, ww2, ww3 ); \
 \
         MULT4( v00, v01, v02, v03, ww0, ww1, ww2, ww3,GLUE(c,0));\
         MULT4( v10, v11, v12, v13, ww0, ww1, ww2, ww3,GLUE(c,1));\
         MULT4( v20, v21, v22, v23, ww0, ww1, ww2, ww3,GLUE(c,2));\
         MULT4( v30, v31, v32, v33, ww0, ww1, ww2, ww3,GLUE(c,3));\
 \
         MULT4( v0, v1, v2, v3, wv0, wv1, wv2, wv3,v); \
 \
         (val) = CUBIC_UNIVAR( v0, v1, v2, v3, u_parm ); \
    }

#define  CUBIC_TRIVAR_DERIV( c, u_parm, v_parm, w_parm, val, deriv_u, deriv_v, deriv_w ) \
    { \
         Real  wu0, wu1, wu2, wu3; \
         Real  wv0, wv1, wv2, wv3; \
         Real  wdv0, wdv1, wdv2, wdv3; \
         Real  ww0, ww1, ww2, ww3; \
         Real  wdw0, wdw1, wdw2, wdw3; \
         Real  v00, v01, v02, v03, v10, v11, v12, v13; \
         Real  v20, v21, v22, v23, v30, v31, v32, v33; \
         Real  dw00, dw01, dw02, dw03, dw10, dw11, dw12, dw13; \
         Real  dw20, dw21, dw22, dw23, dw30, dw31, dw32, dw33; \
         Real  v0, v1, v2, v3; \
         Real  dv0, dv1, dv2, dv3; \
         Real  dw0, dw1, dw2, dw3; \
 \
         COMPUTE_CUBIC_COEFFS( u_parm, wu0, wu1, wu2, wu3 ); \
         COMPUTE_CUBIC_COEFFS( v_parm, wv0, wv1, wv2, wv3 ); \
         COMPUTE_CUBIC_COEFFS( w_parm, ww0, ww1, ww2, ww3 ); \
         COMPUTE_CUBIC_DERIV_COEFFS( v_parm, wdv0, wdv1, wdv2, wdv3 ); \
         COMPUTE_CUBIC_DERIV_COEFFS( w_parm, wdw0, wdw1, wdw2, wdw3 ); \
 \
         MULT4( v00, v01, v02, v03, ww0, ww1, ww2, ww3,GLUE(c,0)); \
         MULT4( v10, v11, v12, v13, ww0, ww1, ww2, ww3,GLUE(c,1)); \
         MULT4( v20, v21, v22, v23, ww0, ww1, ww2, ww3,GLUE(c,2)); \
         MULT4( v30, v31, v32, v33, ww0, ww1, ww2, ww3,GLUE(c,3)); \
 \
         MULT4( dw00, dw01, dw02, dw03, wdw0, wdw1, wdw2, wdw3,GLUE(c,0));\
         MULT4( dw10, dw11, dw12, dw13, wdw0, wdw1, wdw2, wdw3,GLUE(c,1));\
         MULT4( dw20, dw21, dw22, dw23, wdw0, wdw1, wdw2, wdw3,GLUE(c,2));\
         MULT4( dw30, dw31, dw32, dw33, wdw0, wdw1, wdw2, wdw3,GLUE(c,3));\
 \
         MULT4( v0, v1, v2, v3, wv0, wv1, wv2, wv3,v); \
         MULT4( dv0, dv1, dv2, dv3, wdv0, wdv1, wdv2, wdv3,v);\
         MULT4( dw0, dw1, dw2, dw3, wv0, wv1, wv2, wv3,dw);\
 \
         (val) = DOT4( v0, v1, v2, v3, wu0, wu1, wu2, wu3 ); \
         (deriv_u) = CUBIC_UNIVAR_DERIV( v0, v1, v2, v3, u_parm ); \
         (deriv_v) = DOT4( dv0, dv1, dv2, dv3, wu0, wu1, wu2, wu3 ); \
         (deriv_w) = DOT4( dw0, dw1, dw2, dw3, wu0, wu1, wu2, wu3 ); \
    }

#define  CUBIC_TRIVAR_DERIV2( c, u_parm, v_parm, w_parm, val, deriv_u, deriv_v, deriv_w, deriv_uu, deriv_uv, deriv_uw, deriv_vv, deriv_vw, deriv_ww ) \
    { \
         Real  wu0, wu1, wu2, wu3; \
         Real  wdu0, wdu1, wdu2, wdu3; \
         Real  wv0, wv1, wv2, wv3; \
         Real  wdv0, wdv1, wdv2, wdv3; \
         Real  wdvv0, wdvv1, wdvv2, wdvv3; \
         Real  ww0, ww1, ww2, ww3; \
         Real  wdw0, wdw1, wdw2, wdw3; \
         Real  wdww0, wdww1, wdww2, wdww3; \
         Real  v00, v01, v02, v03, v10, v11, v12, v13; \
         Real  v20, v21, v22, v23, v30, v31, v32, v33; \
         Real  dw00, dw01, dw02, dw03, dw10, dw11, dw12, dw13; \
         Real  dw20, dw21, dw22, dw23, dw30, dw31, dw32, dw33; \
         Real  dww00, dww01, dww02, dww03, dww10, dww11, dww12, dww13; \
         Real  dww20, dww21, dww22, dww23, dww30, dww31, dww32, dww33; \
         Real  v0, v1, v2, v3; \
         Real  dv0, dv1, dv2, dv3; \
         Real  dvv0, dvv1, dvv2, dvv3; \
         Real  dw0, dw1, dw2, dw3; \
         Real  dww0, dww1, dww2, dww3; \
         Real  dvw0, dvw1, dvw2, dvw3; \
 \
         COMPUTE_CUBIC_COEFFS( u_parm, wu0, wu1, wu2, wu3 ); \
         COMPUTE_CUBIC_COEFFS( v_parm, wv0, wv1, wv2, wv3 ); \
         COMPUTE_CUBIC_COEFFS( w_parm, ww0, ww1, ww2, ww3 ); \
         COMPUTE_CUBIC_DERIV_COEFFS( u_parm, wdu0, wdu1, wdu2, wdu3 ); \
         COMPUTE_CUBIC_DERIV_COEFFS( v_parm, wdv0, wdv1, wdv2, wdv3 ); \
         COMPUTE_CUBIC_DERIV_COEFFS( w_parm, wdw0, wdw1, wdw2, wdw3 ); \
         COMPUTE_CUBIC_DERIV2_COEFFS( v_parm, wdvv0, wdvv1, wdvv2, wdvv3 ); \
         COMPUTE_CUBIC_DERIV2_COEFFS( w_parm, wdww0, wdww1, wdww2, wdww3 ); \
 \
         MULT4( v00, v01, v02, v03, ww0, ww1, ww2, ww3,GLUE(c,0)); \
         MULT4( v10, v11, v12, v13, ww0, ww1, ww2, ww3,GLUE(c,1)); \
         MULT4( v20, v21, v22, v23, ww0, ww1, ww2, ww3,GLUE(c,2)); \
         MULT4( v30, v31, v32, v33, ww0, ww1, ww2, ww3,GLUE(c,3)); \
 \
         MULT4( dw00, dw01, dw02, dw03, wdw0, wdw1, wdw2, wdw3,GLUE(c,0)); \
         MULT4( dw10, dw11, dw12, dw13, wdw0, wdw1, wdw2, wdw3,GLUE(c,1)); \
         MULT4( dw20, dw21, dw22, dw23, wdw0, wdw1, wdw2, wdw3,GLUE(c,2)); \
         MULT4( dw30, dw31, dw32, dw33, wdw0, wdw1, wdw2, wdw3,GLUE(c,3)); \
 \
         MULT4( dww00, dww01, dww02, dww03, wdww0, wdww1, wdww2, wdww3,GLUE(c,0)); \
         MULT4( dww10, dww11, dww12, dww13, wdww0, wdww1, wdww2, wdww3,GLUE(c,1)); \
         MULT4( dww20, dww21, dww22, dww23, wdww0, wdww1, wdww2, wdww3,GLUE(c,2)); \
         MULT4( dww30, dww31, dww32, dww33, wdww0, wdww1, wdww2, wdww3,GLUE(c,3)); \
 \
         MULT4( v0, v1, v2, v3, wv0, wv1, wv2, wv3,v); \
         MULT4( dv0, dv1, dv2, dv3, wdv0, wdv1, wdv2, wdv3,v); \
         MULT4( dvv0, dvv1, dvv2, dvv3, wdvv0, wdvv1, wdvv2, wdvv3,v); \
         MULT4( dw0, dw1, dw2, dw3, wv0, wv1, wv2, wv3,dw); \
         MULT4( dww0, dww1, dww2, dww3, wv0, wv1, wv2, wv3,dww); \
         MULT4( dvw0, dvw1, dvw2, dvw3, wdv0, wdv1, wdv2, wdv3,dw); \
 \
         (val) = DOT4( v0, v1, v2, v3,         wu0, wu1, wu2, wu3 ); \
         (deriv_v)  = DOT4( dv0, dv1, dv2, dv3,     wu0, wu1, wu2, wu3 ); \
         (deriv_w)  = DOT4( dw0, dw1, dw2, dw3,     wu0, wu1, wu2, wu3 ); \
         (deriv_vv) = DOT4( dvv0, dvv1, dvv2, dvv3, wu0, wu1, wu2, wu3 ); \
         (deriv_vw) = DOT4( dvw0, dvw1, dvw2, dvw3, wu0, wu1, wu2, wu3 ); \
         (deriv_ww) = DOT4( dww0, dww1, dww2, dww3, wu0, wu1, wu2, wu3 ); \
         (deriv_u)  = DOT4( v0, v1, v2, v3,         wdu0, wdu1, wdu2, wdu3 ); \
         (deriv_uv) = DOT4( dv0, dv1, dv2, dv3,     wdu0, wdu1, wdu2, wdu3 ); \
         (deriv_uw) = DOT4( dw0, dw1, dw2, dw3,     wdu0, wdu1, wdu2, wdu3 ); \
         (deriv_uu) = CUBIC_UNIVAR_DERIV2( v0, v1, v2, v3, u_parm ); \
    }

#endif
