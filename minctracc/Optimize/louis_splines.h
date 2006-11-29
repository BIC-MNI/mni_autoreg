
#ifndef  DEF_SPLINES
#define  DEF_SPLINES

#include  <volume_io/basic.h>

#define  LINEAR_COEF_00     1.0
#define  LINEAR_COEF_01     0.0
#define  LINEAR_COEF_10    -1.0
#define  LINEAR_COEF_11     1.0

#define  LINEAR_COEF_0( v0, v1 ) \
            (v0)
#define  LINEAR_COEF_1( v0, v1 ) \
            ((v1) - (v0))

#define LINEAR_UNIVAR( v0, v1, u ) \
            (LINEAR_COEF_1(v0,v1 ) * u + LINEAR_COEF_0(v0,v1))

/* ------------------ QUADRATIC INTERPOLATING SPLINES ----------------------- */

#define  QUADRATIC_COEF_00   0.5
#define  QUADRATIC_COEF_01   0.5
#define  QUADRATIC_COEF_02   0.0
#define  QUADRATIC_COEF_10  -1.5
#define  QUADRATIC_COEF_11   2.0
#define  QUADRATIC_COEF_12  -0.5
#define  QUADRATIC_COEF_20   1.0
#define  QUADRATIC_COEF_21  -2.0
#define  QUADRATIC_COEF_22   1.0

#define  QUADRATIC_COEF_0( v0, v1, v2 ) \
          (QUADRATIC_COEF_00*(v0) + QUADRATIC_COEF_01*(v1))
#define  QUADRATIC_COEF_1( v0, v1, v2 ) \
          (QUADRATIC_COEF_10*(v0)+QUADRATIC_COEF_11*(v1)+QUADRATIC_COEF_12*(v2))
#define  QUADRATIC_COEF_2( v0, v1, v2 ) \
          ((v0)+QUADRATIC_COEF_21*(v1)+(v2))

#define  QUADRATIC_UNIVAR( v0, v1, v2, u )                               \
     ( QUADRATIC_COEF_0(v0,v1,v2) + (u) * \
       (QUADRATIC_COEF_1( v0, v1, v2 ) + (u) * QUADRATIC_COEF_2(v0,v1,v2)) )

#define  QUADRATIC_UNIVAR_DERIV( v0, v1, v2, u ) \
     ( QUADRATIC_COEF_1(v0,v1,v2) + (u) * 2.0 * QUADRATIC_COEF_2(v0,v1,v2) )

#define  QUADRATIC_UNIVAR_DERIV2( v0, v1, v2, u ) \
     ( 2.0 * QUADRATIC_COEF_2(v0,v1,v2) )

#define  QUADRATIC_BIVAR( v00, v01, v02, v10, v11, v12, v20, v21, v22, u, v, val ) \
     { \
         double  cu0, cu1, cu2; \
 \
         cu0 = QUADRATIC_UNIVAR( v00, v01, v02, v ); \
         cu1 = QUADRATIC_UNIVAR( v10, v11, v12, v ); \
         cu2 = QUADRATIC_UNIVAR( v20, v21, v22, v ); \
 \
         val = QUADRATIC_UNIVAR( cu0, cu1, cu2, u ); \
     }

#define  QUADRATIC_BIVAR_DERIV( v00, v01, v02, v10, v11, v12, v20, v21, v22, u, v, du, dv ) \
     { \
         double  cu0, cu1, cu2; \
         double  dv0, dv1, dv2; \
 \
         dv0 = QUADRATIC_UNIVAR_DERIV( v00, v01, v02, v ); \
         dv1 = QUADRATIC_UNIVAR_DERIV( v10, v11, v12, v ); \
         dv2 = QUADRATIC_UNIVAR_DERIV( v20, v21, v22, v ); \
 \
         dv = QUADRATIC_UNIVAR( dv0, dv1, dv2, u ); \
 \
         cu0 = QUADRATIC_UNIVAR( v00, v01, v02, v ); \
         cu1 = QUADRATIC_UNIVAR( v10, v11, v12, v ); \
         cu2 = QUADRATIC_UNIVAR( v20, v21, v22, v ); \
 \
         du = QUADRATIC_UNIVAR_DERIV( cu0, cu1, cu2, u ); \
     }

#define  QUADRATIC_BIVAR_DERIV2( v00, v01, v02, v10, v11, v12, v20, v21, v22, u, v, val, du, dv, duu, duv, dvv ) \
     { \
         double  cu0, cu1, cu2; \
         double  dv0, dv1, dv2; \
         double  dvv0, dvv1, dvv2; \
 \
         dv0 = QUADRATIC_UNIVAR_DERIV( v00, v01, v02, v ); \
         dv1 = QUADRATIC_UNIVAR_DERIV( v10, v11, v12, v ); \
         dv2 = QUADRATIC_UNIVAR_DERIV( v20, v21, v22, v ); \
 \
         dvv0 = QUADRATIC_UNIVAR_DERIV2( v00, v01, v02, v ); \
         dvv1 = QUADRATIC_UNIVAR_DERIV2( v10, v11, v12, v ); \
         dvv2 = QUADRATIC_UNIVAR_DERIV2( v20, v21, v22, v ); \
 \
         dv = QUADRATIC_UNIVAR( dv0, dv1, dv2, u ); \
         dvv = QUADRATIC_UNIVAR( dvv0, dvv1, dvv2, u ); \
 \
         duv = QUADRATIC_UNIVAR_DERIV( dv0, dv1, dv2, u ); \
 \
         cu0 = QUADRATIC_UNIVAR( v00, v01, v02, v ); \
         cu1 = QUADRATIC_UNIVAR( v10, v11, v12, v ); \
         cu2 = QUADRATIC_UNIVAR( v20, v21, v22, v ); \
 \
         val = QUADRATIC_UNIVAR( cu0, cu1, cu2, u ); \
         du = QUADRATIC_UNIVAR_DERIV( cu0, cu1, cu2, u ); \
         duu = QUADRATIC_UNIVAR_DERIV2( cu0, cu1, cu2, u ); \
     }

#define  QUADRATIC_TRIVAR( c, u, v, w, val ) \
    { \
        double  c00, c01, c02, c10, c11, c12, c20, c21, c22; \
        double  c0, c1, c2; \
 \
        c00 = QUADRATIC_UNIVAR(GLUE(c,000),GLUE(c,001),GLUE(c,002), w ); \
        c01 = QUADRATIC_UNIVAR(GLUE(c,010),GLUE(c,011),GLUE(c,012), w ); \
        c02 = QUADRATIC_UNIVAR(GLUE(c,020),GLUE(c,021),GLUE(c,022), w ); \
 \
        c10 = QUADRATIC_UNIVAR(GLUE(c,100),GLUE(c,101),GLUE(c,102), w ); \
        c11 = QUADRATIC_UNIVAR(GLUE(c,110),GLUE(c,111),GLUE(c,112), w ); \
        c12 = QUADRATIC_UNIVAR(GLUE(c,120),GLUE(c,121),GLUE(c,122), w ); \
 \
        c20 = QUADRATIC_UNIVAR(GLUE(c,200),GLUE(c,201),GLUE(c,202), w ); \
        c21 = QUADRATIC_UNIVAR(GLUE(c,210),GLUE(c,211),GLUE(c,212), w ); \
        c22 = QUADRATIC_UNIVAR(GLUE(c,220),GLUE(c,221),GLUE(c,222), w ); \
 \
        c0 = QUADRATIC_UNIVAR( c00, c01, c02, v ); \
        c1 = QUADRATIC_UNIVAR( c10, c11, c12, v ); \
        c2 = QUADRATIC_UNIVAR( c20, c21, c22, v ); \
 \
        (val) = QUADRATIC_UNIVAR( c0, c1, c2, u ); \
    }

#define  QUADRATIC_TRIVAR_DERIV( c, u, v, w, du, dv, dw ) \
    { \
        double  c00, c01, c02, c10, c11, c12, c20, c21, c22; \
        double  dw00, dw01, dw02, dw10, dw11, dw12, dw20, dw21, dw22; \
        double  dv0, dv1, dv2; \
        double  dw0, dw1, dw2; \
        double  c0, c1, c2; \
 \
        c00 = QUADRATIC_UNIVAR(GLUE(c,000),GLUE(c,001),GLUE(c,002), w ); \
        c01 = QUADRATIC_UNIVAR(GLUE(c,010),GLUE(c,011),GLUE(c,012), w ); \
        c02 = QUADRATIC_UNIVAR(GLUE(c,020),GLUE(c,021),GLUE(c,022), w ); \
 \
        c10 = QUADRATIC_UNIVAR(GLUE(c,100),GLUE(c,101),GLUE(c,102), w ); \
        c11 = QUADRATIC_UNIVAR(GLUE(c,110),GLUE(c,111),GLUE(c,112), w ); \
        c12 = QUADRATIC_UNIVAR(GLUE(c,120),GLUE(c,121),GLUE(c,122), w ); \
 \
        c20 = QUADRATIC_UNIVAR(GLUE(c,200),GLUE(c,201),GLUE(c,202), w ); \
        c21 = QUADRATIC_UNIVAR(GLUE(c,210),GLUE(c,211),GLUE(c,212), w ); \
        c22 = QUADRATIC_UNIVAR(GLUE(c,220),GLUE(c,221),GLUE(c,222), w ); \
 \
        dw00 = QUADRATIC_UNIVAR_DERIV(GLUE(c,000),GLUE(c,001),GLUE(c,002), w ); \
        dw01 = QUADRATIC_UNIVAR_DERIV(GLUE(c,010),GLUE(c,011),GLUE(c,012), w ); \
        dw02 = QUADRATIC_UNIVAR_DERIV(GLUE(c,020),GLUE(c,021),GLUE(c,022), w ); \
 \
        dw10 = QUADRATIC_UNIVAR_DERIV(GLUE(c,100),GLUE(c,101),GLUE(c,102), w ); \
        dw11 = QUADRATIC_UNIVAR_DERIV(GLUE(c,110),GLUE(c,111),GLUE(c,112), w ); \
        dw12 = QUADRATIC_UNIVAR_DERIV(GLUE(c,120),GLUE(c,121),GLUE(c,122), w ); \
 \
        dw20 = QUADRATIC_UNIVAR_DERIV(GLUE(c,200),GLUE(c,201),GLUE(c,202), w ); \
        dw21 = QUADRATIC_UNIVAR_DERIV(GLUE(c,210),GLUE(c,211),GLUE(c,212), w ); \
        dw22 = QUADRATIC_UNIVAR_DERIV(GLUE(c,220),GLUE(c,221),GLUE(c,222), w ); \
 \
        c0 = QUADRATIC_UNIVAR( c00, c01, c02, v ); \
        c1 = QUADRATIC_UNIVAR( c10, c11, c12, v ); \
        c2 = QUADRATIC_UNIVAR( c20, c21, c22, v ); \
 \
        dv0 = QUADRATIC_UNIVAR_DERIV( c00, c01, c02, v ); \
        dv1 = QUADRATIC_UNIVAR_DERIV( c10, c11, c12, v ); \
        dv2 = QUADRATIC_UNIVAR_DERIV( c20, c21, c22, v ); \
 \
        dw0 = QUADRATIC_UNIVAR( dw00, dw01, dw02, v ); \
        dw1 = QUADRATIC_UNIVAR( dw10, dw11, dw12, v ); \
        dw2 = QUADRATIC_UNIVAR( dw20, dw21, dw22, v ); \
 \
        (du) = QUADRATIC_UNIVAR_DERIV( c0, c1, c2, u ); \
        (dv) = QUADRATIC_UNIVAR( dv0, dv1, dv2, u ); \
        (dw) = QUADRATIC_UNIVAR( dw0, dw1, dw2, u ); \
    }

#define  QUADRATIC_TRIVAR_DERIV2( c, u, v, w, duu, duv, duw, dvv, dvw, dww ) \
    { \
        double  c00, c01, c02, c10, c11, c12, c20, c21, c22; \
        double  c0, c1, c2; \
        double  dw00, dw01, dw02, dw10, dw11, dw12, dw20, dw21, dw22; \
        double  dww00, dww01, dww02, dww10, dww11, dww12, dww20, dww21, dww22; \
        double  dww0, dww1, dww2; \
        double  dv0, dv1, dv2; \
        double  dvv0, dvv1, dvv2; \
        double  dvw0, dvw1, dvw2; \
        double  dw0, dw1, dw2; \
 \
        c00 = QUADRATIC_UNIVAR(GLUE(c,000),GLUE(c,001),GLUE(c,002),w); \
        c01 = QUADRATIC_UNIVAR(GLUE(c,010),GLUE(c,011),GLUE(c,012),w); \
        c02 = QUADRATIC_UNIVAR(GLUE(c,020),GLUE(c,021),GLUE(c,022),w); \
 \
        c10 = QUADRATIC_UNIVAR(GLUE(c,100),GLUE(c,101),GLUE(c,102),w); \
        c11 = QUADRATIC_UNIVAR(GLUE(c,110),GLUE(c,111),GLUE(c,112),w); \
        c12 = QUADRATIC_UNIVAR(GLUE(c,120),GLUE(c,121),GLUE(c,122),w); \
 \
        c20 = QUADRATIC_UNIVAR(GLUE(c,200),GLUE(c,201),GLUE(c,202),w); \
        c21 = QUADRATIC_UNIVAR(GLUE(c,210),GLUE(c,211),GLUE(c,212),w); \
        c22 = QUADRATIC_UNIVAR(GLUE(c,220),GLUE(c,221),GLUE(c,222),w); \
 \
        dw00 = QUADRATIC_UNIVAR_DERIV(GLUE(c,000),GLUE(c,001),GLUE(c,002),w);\
        dw01 = QUADRATIC_UNIVAR_DERIV(GLUE(c,010),GLUE(c,011),GLUE(c,012),w);\
        dw02 = QUADRATIC_UNIVAR_DERIV(GLUE(c,020),GLUE(c,021),GLUE(c,022),w);\
\
        dw10 = QUADRATIC_UNIVAR_DERIV(GLUE(c,100),GLUE(c,101),GLUE(c,102),w);\
        dw11 = QUADRATIC_UNIVAR_DERIV(GLUE(c,110),GLUE(c,111),GLUE(c,112),w);\
        dw12 = QUADRATIC_UNIVAR_DERIV(GLUE(c,120),GLUE(c,121),GLUE(c,122),w);\
\
        dw20 = QUADRATIC_UNIVAR_DERIV(GLUE(c,200),GLUE(c,201),GLUE(c,202),w);\
        dw21 = QUADRATIC_UNIVAR_DERIV(GLUE(c,210),GLUE(c,211),GLUE(c,212),w);\
        dw22 = QUADRATIC_UNIVAR_DERIV(GLUE(c,220),GLUE(c,221),GLUE(c,222),w);\
 \
 \
        dww00 = QUADRATIC_UNIVAR_DERIV2(GLUE(c,000),GLUE(c,001),GLUE(c,002),w);\
        dww01 = QUADRATIC_UNIVAR_DERIV2(GLUE(c,010),GLUE(c,011),GLUE(c,012),w);\
        dww02 = QUADRATIC_UNIVAR_DERIV2(GLUE(c,020),GLUE(c,021),GLUE(c,022),w);\
\
        dww10 = QUADRATIC_UNIVAR_DERIV2(GLUE(c,100),GLUE(c,101),GLUE(c,102),w);\
        dww11 = QUADRATIC_UNIVAR_DERIV2(GLUE(c,110),GLUE(c,111),GLUE(c,112),w);\
        dww12 = QUADRATIC_UNIVAR_DERIV2(GLUE(c,120),GLUE(c,121),GLUE(c,122),w);\
\
        dww20 = QUADRATIC_UNIVAR_DERIV2(GLUE(c,200),GLUE(c,201),GLUE(c,202),w);\
        dww21 = QUADRATIC_UNIVAR_DERIV2(GLUE(c,210),GLUE(c,211),GLUE(c,212),w);\
        dww22 = QUADRATIC_UNIVAR_DERIV2(GLUE(c,220),GLUE(c,221),GLUE(c,222),w);\
\
        c0 = QUADRATIC_UNIVAR( c00, c01, c02, v ); \
        c1 = QUADRATIC_UNIVAR( c10, c11, c12, v ); \
        c2 = QUADRATIC_UNIVAR( c20, c21, c22, v ); \
 \
        dv0 = QUADRATIC_UNIVAR_DERIV( c00, c01, c02, v ); \
        dv1 = QUADRATIC_UNIVAR_DERIV( c10, c11, c12, v ); \
        dv2 = QUADRATIC_UNIVAR_DERIV( c20, c21, c22, v ); \
 \
        dvv0 = QUADRATIC_UNIVAR_DERIV2( c00, c01, c02, v ); \
        dvv1 = QUADRATIC_UNIVAR_DERIV2( c10, c11, c12, v ); \
        dvv2 = QUADRATIC_UNIVAR_DERIV2( c20, c21, c22, v ); \
 \
        dw0 = QUADRATIC_UNIVAR( dw00, dw01, dw02, v ); \
        dw1 = QUADRATIC_UNIVAR( dw10, dw11, dw12, v ); \
        dw2 = QUADRATIC_UNIVAR( dw20, dw21, dw22, v ); \
 \
        dvw0 = QUADRATIC_UNIVAR_DERIV( dw00, dw01, dw02, v ); \
        dvw1 = QUADRATIC_UNIVAR_DERIV( dw10, dw11, dw12, v ); \
        dvw2 = QUADRATIC_UNIVAR_DERIV( dw20, dw21, dw22, v ); \
 \
        dww0 = QUADRATIC_UNIVAR( dww00, dww01, dww02, v ); \
        dww1 = QUADRATIC_UNIVAR( dww10, dww11, dww12, v ); \
        dww2 = QUADRATIC_UNIVAR( dww20, dww21, dww22, v ); \
 \
        (duu) = QUADRATIC_UNIVAR_DERIV2( c0, c1, c2, u ); \
        (duv) = QUADRATIC_UNIVAR_DERIV( dv0, dv1, dv2, u ); \
        (duw) = QUADRATIC_UNIVAR_DERIV( dw0, dw1, dw2, u ); \
        (dvv) = QUADRATIC_UNIVAR( dvv0, dvv1, dvv2, u ); \
        (dvw) = QUADRATIC_UNIVAR( dvw0, dvw1, dvw2, u ); \
        (dww) = QUADRATIC_UNIVAR( dww0, dww1, dww2, u ); \
    }

/* -------------------- CUBIC INTERPOLATING SPLINES ----------------------- */

#define  CUBIC_COEF_00   0.0
#define  CUBIC_COEF_01   1.0
#define  CUBIC_COEF_02   0.0
#define  CUBIC_COEF_03   0.0

#define  CUBIC_COEF_10  -0.5
#define  CUBIC_COEF_11   0.0
#define  CUBIC_COEF_12   0.5
#define  CUBIC_COEF_13   0.0

#define  CUBIC_COEF_20   1.0
#define  CUBIC_COEF_21  -2.5
#define  CUBIC_COEF_22   2.0
#define  CUBIC_COEF_23  -0.5

#define  CUBIC_COEF_30  -0.5
#define  CUBIC_COEF_31   1.5
#define  CUBIC_COEF_32  -1.5
#define  CUBIC_COEF_33   0.5

#define  CUBIC_COEF_0( v0, v1, v2, v3 ) \
             (v1)
#define  CUBIC_COEF_1( v0, v1, v2, v3 ) \
             (0.5 * ((v2)-(v0)))
#define  CUBIC_COEF_2( v0, v1, v2, v3 ) \
             ((v0) - 2.5 * (v1) + 2.0 * (v2) - 0.5 * (v3))
#define  CUBIC_COEF_3( v0, v1, v2, v3 ) \
             (-0.5 * (v0) + 1.5 * (v1) - 1.5 * (v2) + 0.5 * (v3))

#define  CUBIC_UNIVAR( v0, v1, v2, v3, u ) \
     (  CUBIC_COEF_0(v0,v1,v2,v3) + (u) * ( \
        CUBIC_COEF_1(v0,v1,v2,v3) + (u) * ( \
        CUBIC_COEF_2(v0,v1,v2,v3) + (u) * \
        CUBIC_COEF_3(v0,v1,v2,v3) ) ) \
     )

#define  CUBIC_UNIVAR_DERIV( v0, v1, v2, v3, u ) \
     ( \
        0.5 * ((v2) - (v0)) + (u) * ( \
        2.0 * (v0) - 5.0 * (v1) + 4.0 * (v2) - (v3) + (u) * ( \
        -1.5 * (v0) + 4.5 * (v1) - 4.5 * (v2) + 1.5 * (v3)  ) \
                                    ) \
     )

#define  CUBIC_UNIVAR_DERIV2( v0, v1, v2, v3, u ) \
     ( \
        2.0 * (v0) - 5.0 * (v1) + 4.0 * (v2) - (v3) + (u) * ( \
        -3.0 * (v0) + 9.0 * (v1) - 9.0 * (v2) + 3.0 * (v3)  ) \
     )

#define  CUBIC_UNIVAR_VAL_DERIV1( v0, v1, v2, v3, u, val, deriv )              \
     {                                                                         \
         VIO_Real   _c1, _c2, _c3;                                                 \
                                                                               \
         _c1 = CUBIC_COEF_1(v0,v1,v2,v3);                                      \
         _c2 = CUBIC_COEF_2(v0,v1,v2,v3);                                      \
         _c3 = CUBIC_COEF_3(v0,v1,v2,v3);                                      \
         (val) = CUBIC_COEF_0(v0,v1,v2,v3) + (u) * (_c1+(u)*(_c2+(u)*_c3));    \
         (deriv) = _c1 + (u) * (2.0*_c2 + 3.0 * (u) * _c3);                    \
     }

#define  CUBIC_UNIVAR_VAL_DERIV2( v0, v1, v2, v3, u, val, deriv1, deriv2 )     \
     {                                                                         \
         VIO_Real   _c1, _c2, _c3;                                                 \
                                                                               \
         _c1 = CUBIC_COEF_1(v0,v1,v2,v3);                                      \
         _c2 = CUBIC_COEF_2(v0,v1,v2,v3);                                      \
         _c3 = CUBIC_COEF_3(v0,v1,v2,v3);                                      \
         (val) = CUBIC_COEF_0(v0,v1,v2,v3) + (u) * (_c1 + (u) * (_c2+(u)*_c3));\
         (deriv1) = _c1 + (u) * (2.0*_c2 + 3.0 * (u) * _c3);                   \
         (deriv2) = 2.0 * _c2 + 6.0 * (u) * _c3;                               \
     }

#define  CUBIC_BIVAR( v00, v01, v02, v03, v10, v11, v12, v13, v20, v21, v22, v23, v30, v31, v32, v33, u, v, val ) \
     { \
         double  cu0, cu1, cu2, cu3; \
 \
         cu0 = CUBIC_UNIVAR( v00, v01, v02, v03, v ); \
         cu1 = CUBIC_UNIVAR( v10, v11, v12, v13, v ); \
         cu2 = CUBIC_UNIVAR( v20, v21, v22, v23, v ); \
         cu3 = CUBIC_UNIVAR( v30, v31, v32, v33, v ); \
 \
         val = CUBIC_UNIVAR( cu0, cu1, cu2, cu3, u ); \
     }

#define  CUBIC_BIVAR_DERIV( v00, v01, v02, v03, v10, v11, v12, v13, v20, v21, v22, v23, v30, v31, v32, v33, u, v, du, dv ) \
     { \
         double  cu0, cu1, cu2, cu3; \
         double  dv0, dv1, dv2, dv3; \
 \
         dv0 = CUBIC_UNIVAR_DERIV( v00, v01, v02, v03, v ); \
         dv1 = CUBIC_UNIVAR_DERIV( v10, v11, v12, v13, v ); \
         dv2 = CUBIC_UNIVAR_DERIV( v20, v21, v22, v23, v ); \
         dv3 = CUBIC_UNIVAR_DERIV( v30, v31, v32, v33, v ); \
 \
         dv = CUBIC_UNIVAR( dv0, dv1, dv2, dv3, u ); \
 \
         cu0 = CUBIC_UNIVAR( v00, v01, v02, v03, v ); \
         cu1 = CUBIC_UNIVAR( v10, v11, v12, v13, v ); \
         cu2 = CUBIC_UNIVAR( v20, v21, v22, v23, v ); \
         cu3 = CUBIC_UNIVAR( v30, v31, v32, v33, v ); \
 \
         du = CUBIC_UNIVAR_DERIV( cu0, cu1, cu2, cu3, u ); \
     }

#define  CUBIC_BIVAR_DERIV2( v00, v01, v02, v03, v10, v11, v12, v13, v20, v21, v22, v23, v30, v31, v32, v33, u, v, duu, duv, dvv ) \
     { \
         double  cu0, cu1, cu2, cu3; \
         double  dv0, dv1, dv2, dv3; \
         double  dvv0, dvv1, dvv2, dvv3; \
 \
         dv0 = CUBIC_UNIVAR_DERIV( v00, v01, v02, v03, v ); \
         dv1 = CUBIC_UNIVAR_DERIV( v10, v11, v12, v13, v ); \
         dv2 = CUBIC_UNIVAR_DERIV( v20, v21, v22, v23, v ); \
         dv3 = CUBIC_UNIVAR_DERIV( v30, v31, v32, v33, v ); \
 \
         dvv0 = CUBIC_UNIVAR_DERIV2( v00, v01, v02, v03, v ); \
         dvv1 = CUBIC_UNIVAR_DERIV2( v10, v11, v12, v13, v ); \
         dvv2 = CUBIC_UNIVAR_DERIV2( v20, v21, v22, v23, v ); \
         dvv3 = CUBIC_UNIVAR_DERIV2( v30, v31, v32, v33, v ); \
 \
         dvv = CUBIC_UNIVAR( dvv0, dvv1, dvv2, dvv3, u ); \
 \
         duv = CUBIC_UNIVAR_DERIV( dv0, dv1, dv2, dv3, u ); \
 \
         cu0 = CUBIC_UNIVAR( v00, v01, v02, v03, v ); \
         cu1 = CUBIC_UNIVAR( v10, v11, v12, v13, v ); \
         cu2 = CUBIC_UNIVAR( v20, v21, v22, v23, v ); \
         cu3 = CUBIC_UNIVAR( v30, v31, v32, v33, v ); \
         duu = CUBIC_UNIVAR_DERIV2( cu0, cu1, cu2, cu3, u ); \
     }

#define  CUBIC_BIVAR_VAL_DERIV2( v00, v01, v02, v03, v10, v11, v12, v13, v20, v21, v22, v23, v30, v31, v32, v33, u, v, val, du, dv, duu, duv, dvv ) \
     {                                                                        \
         double  c0, c1, c2, c3;                                              \
         double  dv0, dv1, dv2, dv3;                                          \
         double  dvv0, dvv1, dvv2, dvv3;                                      \
                                                                              \
         CUBIC_UNIVAR_VAL_DERIV2( v00, v01, v02, v03, v, c0, dv0, dvv0 );     \
         CUBIC_UNIVAR_VAL_DERIV2( v10, v11, v12, v13, v, c1, dv1, dvv1 );     \
         CUBIC_UNIVAR_VAL_DERIV2( v20, v21, v22, v23, v, c2, dv2, dvv2 );     \
         CUBIC_UNIVAR_VAL_DERIV2( v30, v31, v32, v33, v, c3, dv3, dvv3 );     \
                                                                              \
         dvv = CUBIC_UNIVAR( dvv0, dvv1, dvv2, dvv3, u );                     \
         CUBIC_UNIVAR_VAL_DERIV1( dv0, dv1, dv2, dv3, u, dv, duv );           \
         CUBIC_UNIVAR_VAL_DERIV2( c0, c1, c2, c3, u, val, du, duu );          \
     }

#define  CUBIC_TRIVAR( c, u, v, w, val ) \
    { \
        double  c00, c01, c02, c03, c10, c11, c12, c13; \
        double  c20, c21, c22, c23, c30, c31, c32, c33; \
        double  c0, c1, c2, c3; \
 \
        c00 = CUBIC_UNIVAR(GLUE(c,000),GLUE(c,001),GLUE(c,002),GLUE(c,003), w ); \
        c01 = CUBIC_UNIVAR(GLUE(c,010),GLUE(c,011),GLUE(c,012),GLUE(c,013), w ); \
        c02 = CUBIC_UNIVAR(GLUE(c,020),GLUE(c,021),GLUE(c,022),GLUE(c,023), w ); \
        c03 = CUBIC_UNIVAR(GLUE(c,030),GLUE(c,031),GLUE(c,032),GLUE(c,033), w ); \
 \
        c10 = CUBIC_UNIVAR(GLUE(c,100),GLUE(c,101),GLUE(c,102),GLUE(c,103), w ); \
        c11 = CUBIC_UNIVAR(GLUE(c,110),GLUE(c,111),GLUE(c,112),GLUE(c,113), w ); \
        c12 = CUBIC_UNIVAR(GLUE(c,120),GLUE(c,121),GLUE(c,122),GLUE(c,123), w ); \
        c13 = CUBIC_UNIVAR(GLUE(c,130),GLUE(c,131),GLUE(c,132),GLUE(c,133), w ); \
 \
        c20 = CUBIC_UNIVAR(GLUE(c,200),GLUE(c,201),GLUE(c,202),GLUE(c,203), w ); \
        c21 = CUBIC_UNIVAR(GLUE(c,210),GLUE(c,211),GLUE(c,212),GLUE(c,213), w ); \
        c22 = CUBIC_UNIVAR(GLUE(c,220),GLUE(c,221),GLUE(c,222),GLUE(c,223), w ); \
        c23 = CUBIC_UNIVAR(GLUE(c,230),GLUE(c,231),GLUE(c,232),GLUE(c,233), w ); \
 \
        c30 = CUBIC_UNIVAR(GLUE(c,300),GLUE(c,301),GLUE(c,302),GLUE(c,303), w ); \
        c31 = CUBIC_UNIVAR(GLUE(c,310),GLUE(c,311),GLUE(c,312),GLUE(c,313), w ); \
        c32 = CUBIC_UNIVAR(GLUE(c,320),GLUE(c,321),GLUE(c,322),GLUE(c,323), w ); \
        c33 = CUBIC_UNIVAR(GLUE(c,330),GLUE(c,331),GLUE(c,332),GLUE(c,333), w ); \
 \
        c0 = CUBIC_UNIVAR( c00, c01, c02, c03, v ); \
        c1 = CUBIC_UNIVAR( c10, c11, c12, c13, v ); \
        c2 = CUBIC_UNIVAR( c20, c21, c22, c23, v ); \
        c3 = CUBIC_UNIVAR( c30, c31, c32, c33, v ); \
 \
        (val) = CUBIC_UNIVAR( c0, c1, c2, c3, u ); \
    }

#define  CUBIC_TRIVAR_DERIV( c, u, v, w, du, dv, dw ) \
    { \
        double  c00, c01, c02, c03, c10, c11, c12, c13; \
        double  c20, c21, c22, c23, c30, c31, c32, c33; \
        double  dw00, dw01, dw02, dw03, dw10, dw11, dw12, dw13; \
        double  dw20, dw21, dw22, dw23, dw30, dw31, dw32, dw33; \
        double  dv0, dv1, dv2, dv3; \
        double  dw0, dw1, dw2, dw3; \
        double  c0, c1, c2, c3; \
 \
        c00 = CUBIC_UNIVAR(GLUE(c,000),GLUE(c,001),GLUE(c,002),GLUE(c,003), w ); \
        c01 = CUBIC_UNIVAR(GLUE(c,010),GLUE(c,011),GLUE(c,012),GLUE(c,013), w ); \
        c02 = CUBIC_UNIVAR(GLUE(c,020),GLUE(c,021),GLUE(c,022),GLUE(c,023), w ); \
        c03 = CUBIC_UNIVAR(GLUE(c,030),GLUE(c,031),GLUE(c,032),GLUE(c,033), w ); \
 \
        c10 = CUBIC_UNIVAR(GLUE(c,100),GLUE(c,101),GLUE(c,102),GLUE(c,103), w ); \
        c11 = CUBIC_UNIVAR(GLUE(c,110),GLUE(c,111),GLUE(c,112),GLUE(c,113), w ); \
        c12 = CUBIC_UNIVAR(GLUE(c,120),GLUE(c,121),GLUE(c,122),GLUE(c,123), w ); \
        c13 = CUBIC_UNIVAR(GLUE(c,130),GLUE(c,131),GLUE(c,132),GLUE(c,133), w ); \
 \
        c20 = CUBIC_UNIVAR(GLUE(c,200),GLUE(c,201),GLUE(c,202),GLUE(c,203), w ); \
        c21 = CUBIC_UNIVAR(GLUE(c,210),GLUE(c,211),GLUE(c,212),GLUE(c,213), w ); \
        c22 = CUBIC_UNIVAR(GLUE(c,220),GLUE(c,221),GLUE(c,222),GLUE(c,223), w ); \
        c23 = CUBIC_UNIVAR(GLUE(c,230),GLUE(c,231),GLUE(c,232),GLUE(c,233), w ); \
 \
        c30 = CUBIC_UNIVAR(GLUE(c,300),GLUE(c,301),GLUE(c,302),GLUE(c,303), w ); \
        c31 = CUBIC_UNIVAR(GLUE(c,310),GLUE(c,311),GLUE(c,312),GLUE(c,313), w ); \
        c32 = CUBIC_UNIVAR(GLUE(c,320),GLUE(c,321),GLUE(c,322),GLUE(c,323), w ); \
        c33 = CUBIC_UNIVAR(GLUE(c,330),GLUE(c,331),GLUE(c,332),GLUE(c,333), w ); \
 \
        dw00 = CUBIC_UNIVAR_DERIV(GLUE(c,000),GLUE(c,001),GLUE(c,002),GLUE(c,003), w ); \
        dw01 = CUBIC_UNIVAR_DERIV(GLUE(c,010),GLUE(c,011),GLUE(c,012),GLUE(c,013), w ); \
        dw02 = CUBIC_UNIVAR_DERIV(GLUE(c,020),GLUE(c,021),GLUE(c,022),GLUE(c,023), w ); \
        dw03 = CUBIC_UNIVAR_DERIV(GLUE(c,030),GLUE(c,031),GLUE(c,032),GLUE(c,033), w ); \
 \
        dw10 = CUBIC_UNIVAR_DERIV(GLUE(c,100),GLUE(c,101),GLUE(c,102),GLUE(c,103), w ); \
        dw11 = CUBIC_UNIVAR_DERIV(GLUE(c,110),GLUE(c,111),GLUE(c,112),GLUE(c,113), w ); \
        dw12 = CUBIC_UNIVAR_DERIV(GLUE(c,120),GLUE(c,121),GLUE(c,122),GLUE(c,123), w ); \
        dw13 = CUBIC_UNIVAR_DERIV(GLUE(c,130),GLUE(c,131),GLUE(c,132),GLUE(c,133), w ); \
 \
        dw20 = CUBIC_UNIVAR_DERIV(GLUE(c,200),GLUE(c,201),GLUE(c,202),GLUE(c,203), w ); \
        dw21 = CUBIC_UNIVAR_DERIV(GLUE(c,210),GLUE(c,211),GLUE(c,212),GLUE(c,213), w ); \
        dw22 = CUBIC_UNIVAR_DERIV(GLUE(c,220),GLUE(c,221),GLUE(c,222),GLUE(c,223), w ); \
        dw23 = CUBIC_UNIVAR_DERIV(GLUE(c,230),GLUE(c,231),GLUE(c,232),GLUE(c,233), w ); \
 \
        dw30 = CUBIC_UNIVAR_DERIV(GLUE(c,300),GLUE(c,301),GLUE(c,302),GLUE(c,303), w ); \
        dw31 = CUBIC_UNIVAR_DERIV(GLUE(c,310),GLUE(c,311),GLUE(c,312),GLUE(c,313), w ); \
        dw32 = CUBIC_UNIVAR_DERIV(GLUE(c,320),GLUE(c,321),GLUE(c,322),GLUE(c,323), w ); \
        dw33 = CUBIC_UNIVAR_DERIV(GLUE(c,330),GLUE(c,331),GLUE(c,332),GLUE(c,333), w ); \
 \
        c0 = CUBIC_UNIVAR( c00, c01, c02, c03, v ); \
        c1 = CUBIC_UNIVAR( c10, c11, c12, c13, v ); \
        c2 = CUBIC_UNIVAR( c20, c21, c22, c23, v ); \
        c3 = CUBIC_UNIVAR( c30, c31, c32, c33, v ); \
 \
        dv0 = CUBIC_UNIVAR_DERIV( c00, c01, c02, c03, v ); \
        dv1 = CUBIC_UNIVAR_DERIV( c10, c11, c12, c13, v ); \
        dv2 = CUBIC_UNIVAR_DERIV( c20, c21, c22, c23, v ); \
        dv3 = CUBIC_UNIVAR_DERIV( c30, c31, c32, c33, v ); \
 \
        dw0 = CUBIC_UNIVAR( dw00, dw01, dw02, dw03, v ); \
        dw1 = CUBIC_UNIVAR( dw10, dw11, dw12, dw13, v ); \
        dw2 = CUBIC_UNIVAR( dw20, dw21, dw22, dw23, v ); \
        dw3 = CUBIC_UNIVAR( dw30, dw31, dw32, dw33, v ); \
 \
        (du) = CUBIC_UNIVAR_DERIV( c0, c1, c2, c3, u ); \
        (dv) = CUBIC_UNIVAR( dv0, dv1, dv2, dv3, u ); \
        (dw) = CUBIC_UNIVAR( dw0, dw1, dw2, dw3, u ); \
    }

#define  CUBIC_TRIVAR_DERIV2( c, u, v, w, duu, duv, duw, dvv, dvw, dww ) \
    { \
        double  c00, c01, c02, c03, c10, c11, c12, c13; \
        double  c20, c21, c22, c23, c30, c31, c32, c33; \
        double  c0, c1, c2, c3; \
        double  dw00, dw01, dw02, dw03, dw10, dw11, dw12, dw13; \
        double  dw20, dw21, dw22, dw23, dw30, dw31, dw32, dw33; \
        double  dww00, dww01, dww02, dww03, dww10, dww11, dww12, dww13; \
        double  dww20, dww21, dww22, dww23, dww30, dww31, dww32, dww33; \
        double  dww0, dww1, dww2, dww3; \
        double  dv0, dv1, dv2, dv3; \
        double  dvv0, dvv1, dvv2, dvv3; \
        double  dvw0, dvw1, dvw2, dvw3; \
        double  dw0, dw1, dw2, dw3; \
 \
        c00 = CUBIC_UNIVAR(GLUE(c,000),GLUE(c,001),GLUE(c,002),GLUE(c,003),w); \
        c01 = CUBIC_UNIVAR(GLUE(c,010),GLUE(c,011),GLUE(c,012),GLUE(c,013),w); \
        c02 = CUBIC_UNIVAR(GLUE(c,020),GLUE(c,021),GLUE(c,022),GLUE(c,023),w); \
        c03 = CUBIC_UNIVAR(GLUE(c,030),GLUE(c,031),GLUE(c,032),GLUE(c,033),w); \
 \
        c10 = CUBIC_UNIVAR(GLUE(c,100),GLUE(c,101),GLUE(c,102),GLUE(c,103),w); \
        c11 = CUBIC_UNIVAR(GLUE(c,110),GLUE(c,111),GLUE(c,112),GLUE(c,113),w); \
        c12 = CUBIC_UNIVAR(GLUE(c,120),GLUE(c,121),GLUE(c,122),GLUE(c,123),w); \
        c13 = CUBIC_UNIVAR(GLUE(c,130),GLUE(c,131),GLUE(c,132),GLUE(c,133),w); \
 \
        c20 = CUBIC_UNIVAR(GLUE(c,200),GLUE(c,201),GLUE(c,202),GLUE(c,203),w); \
        c21 = CUBIC_UNIVAR(GLUE(c,210),GLUE(c,211),GLUE(c,212),GLUE(c,213),w); \
        c22 = CUBIC_UNIVAR(GLUE(c,220),GLUE(c,221),GLUE(c,222),GLUE(c,223),w); \
        c23 = CUBIC_UNIVAR(GLUE(c,230),GLUE(c,231),GLUE(c,232),GLUE(c,233),w); \
 \
        c30 = CUBIC_UNIVAR(GLUE(c,300),GLUE(c,301),GLUE(c,302),GLUE(c,303),w); \
        c31 = CUBIC_UNIVAR(GLUE(c,310),GLUE(c,311),GLUE(c,312),GLUE(c,313),w); \
        c32 = CUBIC_UNIVAR(GLUE(c,320),GLUE(c,321),GLUE(c,322),GLUE(c,323),w); \
        c33 = CUBIC_UNIVAR(GLUE(c,330),GLUE(c,331),GLUE(c,332),GLUE(c,333),w); \
 \
        dw00 = CUBIC_UNIVAR_DERIV(GLUE(c,000),GLUE(c,001),GLUE(c,002),GLUE(c,003),w);\
        dw01 = CUBIC_UNIVAR_DERIV(GLUE(c,010),GLUE(c,011),GLUE(c,012),GLUE(c,013),w);\
        dw02 = CUBIC_UNIVAR_DERIV(GLUE(c,020),GLUE(c,021),GLUE(c,022),GLUE(c,023),w);\
        dw03 = CUBIC_UNIVAR_DERIV(GLUE(c,030),GLUE(c,031),GLUE(c,032),GLUE(c,033),w);\
\
        dw10 = CUBIC_UNIVAR_DERIV(GLUE(c,100),GLUE(c,101),GLUE(c,102),GLUE(c,103),w);\
        dw11 = CUBIC_UNIVAR_DERIV(GLUE(c,110),GLUE(c,111),GLUE(c,112),GLUE(c,113),w);\
        dw12 = CUBIC_UNIVAR_DERIV(GLUE(c,120),GLUE(c,121),GLUE(c,122),GLUE(c,123),w);\
        dw13 = CUBIC_UNIVAR_DERIV(GLUE(c,130),GLUE(c,131),GLUE(c,132),GLUE(c,133),w);\
\
        dw20 = CUBIC_UNIVAR_DERIV(GLUE(c,200),GLUE(c,201),GLUE(c,202),GLUE(c,203),w);\
        dw21 = CUBIC_UNIVAR_DERIV(GLUE(c,210),GLUE(c,211),GLUE(c,212),GLUE(c,213),w);\
        dw22 = CUBIC_UNIVAR_DERIV(GLUE(c,220),GLUE(c,221),GLUE(c,222),GLUE(c,223),w);\
        dw23 = CUBIC_UNIVAR_DERIV(GLUE(c,230),GLUE(c,231),GLUE(c,232),GLUE(c,233),w);\
\
        dw30 = CUBIC_UNIVAR_DERIV(GLUE(c,300),GLUE(c,301),GLUE(c,302),GLUE(c,303),w);\
        dw31 = CUBIC_UNIVAR_DERIV(GLUE(c,310),GLUE(c,311),GLUE(c,312),GLUE(c,313),w);\
        dw32 = CUBIC_UNIVAR_DERIV(GLUE(c,320),GLUE(c,321),GLUE(c,322),GLUE(c,323),w);\
        dw33 = CUBIC_UNIVAR_DERIV(GLUE(c,330),GLUE(c,331),GLUE(c,332),GLUE(c,333),w);\
 \
 \
        dww00 = CUBIC_UNIVAR_DERIV2(GLUE(c,000),GLUE(c,001),GLUE(c,002),GLUE(c,003),w);\
        dww01 = CUBIC_UNIVAR_DERIV2(GLUE(c,010),GLUE(c,011),GLUE(c,012),GLUE(c,013),w);\
        dww02 = CUBIC_UNIVAR_DERIV2(GLUE(c,020),GLUE(c,021),GLUE(c,022),GLUE(c,023),w);\
        dww03 = CUBIC_UNIVAR_DERIV2(GLUE(c,030),GLUE(c,031),GLUE(c,032),GLUE(c,033),w);\
\
        dww10 = CUBIC_UNIVAR_DERIV2(GLUE(c,100),GLUE(c,101),GLUE(c,102),GLUE(c,103),w);\
        dww11 = CUBIC_UNIVAR_DERIV2(GLUE(c,110),GLUE(c,111),GLUE(c,112),GLUE(c,113),w);\
        dww12 = CUBIC_UNIVAR_DERIV2(GLUE(c,120),GLUE(c,121),GLUE(c,122),GLUE(c,123),w);\
        dww13 = CUBIC_UNIVAR_DERIV2(GLUE(c,130),GLUE(c,131),GLUE(c,132),GLUE(c,133),w);\
\
        dww20 = CUBIC_UNIVAR_DERIV2(GLUE(c,200),GLUE(c,201),GLUE(c,202),GLUE(c,203),w);\
        dww21 = CUBIC_UNIVAR_DERIV2(GLUE(c,210),GLUE(c,211),GLUE(c,212),GLUE(c,213),w);\
        dww22 = CUBIC_UNIVAR_DERIV2(GLUE(c,220),GLUE(c,221),GLUE(c,222),GLUE(c,223),w);\
        dww23 = CUBIC_UNIVAR_DERIV2(GLUE(c,230),GLUE(c,231),GLUE(c,232),GLUE(c,233),w);\
\
        dww30 = CUBIC_UNIVAR_DERIV2(GLUE(c,300),GLUE(c,301),GLUE(c,302),GLUE(c,303),w);\
        dww31 = CUBIC_UNIVAR_DERIV2(GLUE(c,310),GLUE(c,311),GLUE(c,312),GLUE(c,313),w);\
        dww32 = CUBIC_UNIVAR_DERIV2(GLUE(c,320),GLUE(c,321),GLUE(c,322),GLUE(c,323),w);\
        dww33 = CUBIC_UNIVAR_DERIV2(GLUE(c,330),GLUE(c,331),GLUE(c,332),GLUE(c,333),w);\
 \
        c0 = CUBIC_UNIVAR( c00, c01, c02, c03, v ); \
        c1 = CUBIC_UNIVAR( c10, c11, c12, c13, v ); \
        c2 = CUBIC_UNIVAR( c20, c21, c22, c23, v ); \
        c3 = CUBIC_UNIVAR( c30, c31, c32, c33, v ); \
 \
        dv0 = CUBIC_UNIVAR_DERIV( c00, c01, c02, c03, v ); \
        dv1 = CUBIC_UNIVAR_DERIV( c10, c11, c12, c13, v ); \
        dv2 = CUBIC_UNIVAR_DERIV( c20, c21, c22, c23, v ); \
        dv3 = CUBIC_UNIVAR_DERIV( c30, c31, c32, c33, v ); \
 \
        dvv0 = CUBIC_UNIVAR_DERIV2( c00, c01, c02, c03, v ); \
        dvv1 = CUBIC_UNIVAR_DERIV2( c10, c11, c12, c13, v ); \
        dvv2 = CUBIC_UNIVAR_DERIV2( c20, c21, c22, c23, v ); \
        dvv3 = CUBIC_UNIVAR_DERIV2( c30, c31, c32, c33, v ); \
 \
        dw0 = CUBIC_UNIVAR( dw00, dw01, dw02, dw03, v ); \
        dw1 = CUBIC_UNIVAR( dw10, dw11, dw12, dw13, v ); \
        dw2 = CUBIC_UNIVAR( dw20, dw21, dw22, dw23, v ); \
        dw3 = CUBIC_UNIVAR( dw30, dw31, dw32, dw33, v ); \
 \
        dvw0 = CUBIC_UNIVAR_DERIV( dw00, dw01, dw02, dw03, v ); \
        dvw1 = CUBIC_UNIVAR_DERIV( dw10, dw11, dw12, dw13, v ); \
        dvw2 = CUBIC_UNIVAR_DERIV( dw20, dw21, dw22, dw23, v ); \
        dvw3 = CUBIC_UNIVAR_DERIV( dw30, dw31, dw32, dw33, v ); \
 \
        dww0 = CUBIC_UNIVAR( dww00, dww01, dww02, dww03, v ); \
        dww1 = CUBIC_UNIVAR( dww10, dww11, dww12, dww13, v ); \
        dww2 = CUBIC_UNIVAR( dww20, dww21, dww22, dww23, v ); \
        dww3 = CUBIC_UNIVAR( dww30, dww31, dww32, dww33, v ); \
 \
        (duu) = CUBIC_UNIVAR_DERIV2( c0, c1, c2, c3, u ); \
        (duv) = CUBIC_UNIVAR_DERIV( dv0, dv1, dv2, dv3, u ); \
        (duw) = CUBIC_UNIVAR_DERIV( dw0, dw1, dw2, dw3, u ); \
        (dvv) = CUBIC_UNIVAR( dvv0, dvv1, dvv2, dvv3, u ); \
        (dvw) = CUBIC_UNIVAR( dvw0, dvw1, dvw2, dvw3, u ); \
        (dww) = CUBIC_UNIVAR( dww0, dww1, dww2, dww3, u ); \
    }

#define  CUBIC_TRIVAR_VAL_DERIV2( c, u, v, w, val, du, dv, dw, duu, duv, duw, dvv, dvw, dww ) \
    { \
        double  c00, c01, c02, c03, c10, c11, c12, c13; \
        double  c20, c21, c22, c23, c30, c31, c32, c33; \
        double  c0, c1, c2, c3; \
        double  dw00, dw01, dw02, dw03, dw10, dw11, dw12, dw13; \
        double  dw20, dw21, dw22, dw23, dw30, dw31, dw32, dw33; \
        double  dww00, dww01, dww02, dww03, dww10, dww11, dww12, dww13; \
        double  dww20, dww21, dww22, dww23, dww30, dww31, dww32, dww33; \
        double  dww0, dww1, dww2, dww3; \
        double  dv0, dv1, dv2, dv3; \
        double  dvv0, dvv1, dvv2, dvv3; \
        double  dvw0, dvw1, dvw2, dvw3; \
        double  dw0, dw1, dw2, dw3; \
 \
        CUBIC_UNIVAR_VAL_DERIV2(GLUE(c,000),GLUE(c,001),GLUE(c,002),GLUE(c,003),w,c00,dw00,dww00); \
        CUBIC_UNIVAR_VAL_DERIV2(GLUE(c,010),GLUE(c,011),GLUE(c,012),GLUE(c,013),w,c01,dw01,dww01); \
        CUBIC_UNIVAR_VAL_DERIV2(GLUE(c,020),GLUE(c,021),GLUE(c,022),GLUE(c,023),w,c02,dw02,dww02); \
        CUBIC_UNIVAR_VAL_DERIV2(GLUE(c,030),GLUE(c,031),GLUE(c,032),GLUE(c,033),w,c03,dw03,dww03); \
 \
        CUBIC_UNIVAR_VAL_DERIV2(GLUE(c,100),GLUE(c,101),GLUE(c,102),GLUE(c,103),w,c10,dw10,dww10); \
        CUBIC_UNIVAR_VAL_DERIV2(GLUE(c,110),GLUE(c,111),GLUE(c,112),GLUE(c,113),w,c11,dw11,dww11); \
        CUBIC_UNIVAR_VAL_DERIV2(GLUE(c,120),GLUE(c,121),GLUE(c,122),GLUE(c,123),w,c12,dw12,dww12); \
        CUBIC_UNIVAR_VAL_DERIV2(GLUE(c,130),GLUE(c,131),GLUE(c,132),GLUE(c,133),w,c13,dw13,dww13); \
 \
        CUBIC_UNIVAR_VAL_DERIV2(GLUE(c,200),GLUE(c,201),GLUE(c,202),GLUE(c,203),w,c20,dw20,dww20); \
        CUBIC_UNIVAR_VAL_DERIV2(GLUE(c,210),GLUE(c,211),GLUE(c,212),GLUE(c,213),w,c21,dw21,dww21); \
        CUBIC_UNIVAR_VAL_DERIV2(GLUE(c,220),GLUE(c,221),GLUE(c,222),GLUE(c,223),w,c22,dw22,dww22); \
        CUBIC_UNIVAR_VAL_DERIV2(GLUE(c,230),GLUE(c,231),GLUE(c,232),GLUE(c,233),w,c23,dw23,dww23); \
 \
        CUBIC_UNIVAR_VAL_DERIV2(GLUE(c,300),GLUE(c,301),GLUE(c,302),GLUE(c,303),w,c30,dw30,dww30); \
        CUBIC_UNIVAR_VAL_DERIV2(GLUE(c,310),GLUE(c,311),GLUE(c,312),GLUE(c,313),w,c31,dw31,dww31); \
        CUBIC_UNIVAR_VAL_DERIV2(GLUE(c,320),GLUE(c,321),GLUE(c,322),GLUE(c,323),w,c32,dw32,dww32); \
        CUBIC_UNIVAR_VAL_DERIV2(GLUE(c,330),GLUE(c,331),GLUE(c,332),GLUE(c,333),w,c33,dw33,dww33); \
 \
        CUBIC_UNIVAR_VAL_DERIV2( c00, c01, c02, c03, v, c0, dv0, dvv0 ); \
        CUBIC_UNIVAR_VAL_DERIV2( c10, c11, c12, c13, v, c1, dv1, dvv1 ); \
        CUBIC_UNIVAR_VAL_DERIV2( c20, c21, c22, c23, v, c2, dv2, dvv2 ); \
        CUBIC_UNIVAR_VAL_DERIV2( c30, c31, c32, c33, v, c3, dv3, dvv3 ); \
 \
        CUBIC_UNIVAR_VAL_DERIV1( dw00, dw01, dw02, dw03, v, dw0, dvw0 ); \
        CUBIC_UNIVAR_VAL_DERIV1( dw10, dw11, dw12, dw13, v, dw1, dvw1 ); \
        CUBIC_UNIVAR_VAL_DERIV1( dw20, dw21, dw22, dw23, v, dw2, dvw2 ); \
        CUBIC_UNIVAR_VAL_DERIV1( dw30, dw31, dw32, dw33, v, dw3, dvw3 ); \
 \
        dww0 = CUBIC_UNIVAR( dww00, dww01, dww02, dww03, v ); \
        dww1 = CUBIC_UNIVAR( dww10, dww11, dww12, dww13, v ); \
        dww2 = CUBIC_UNIVAR( dww20, dww21, dww22, dww23, v ); \
        dww3 = CUBIC_UNIVAR( dww30, dww31, dww32, dww33, v ); \
 \
        CUBIC_UNIVAR_VAL_DERIV2( c0, c1, c2, c3, u, val, du, duu ); \
        CUBIC_UNIVAR_VAL_DERIV1( dv0, dv1, dv2, dv3, u, dv, duv ); \
        CUBIC_UNIVAR_VAL_DERIV1( dw0, dw1, dw2, dw3, u, dw, duw ); \
        (dvv) = CUBIC_UNIVAR( dvv0, dvv1, dvv2, dvv3, u ); \
        (dvw) = CUBIC_UNIVAR( dvw0, dvw1, dvw2, dvw3, u ); \
        (dww) = CUBIC_UNIVAR( dww0, dww1, dww2, dww3, u ); \
    }

#endif
