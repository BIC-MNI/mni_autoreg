/* '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   derivatives()

        evaluate first and second order derivatives by fitting a quadratic
        to given small neighborhood of values. 

        haralick equal weighted least squares fitting problem, p.394.
        note that the pure second order coefficients are scaled by 0.5 here 
        code was provided by visa (see correlation notes).

        for example, the following function
        f(x,y,z) = x+2*y+3*x*x+4*y*y+2*z*z+5*x*y+7*x*z
        gives cu=1, cv=2, cw=0, cuu=6, cvv=8, cww=4, cuv=5, cuw=7, cvw=0.
        (the assumed form of the second order polynomial is: cu*x + cv*y + cw*z 
        +0.5*cuu*x*x + 0.5*cvv*y*y + 0.5*cww*z*z + cuv*x*y + cuw*x*z + cvw*y*z).
        

''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' */
void derivatives(VIO_Real r[3][3][3], /* compute quadratic fit to this region     */ 
                 VIO_Real c[9])          /* store derivatives here                       */

{
  VIO_Real        slice_u1, slice_u2, slice_u3;
  VIO_Real        slice_v1, slice_v2, slice_v3;
  VIO_Real        slice_w1, slice_w2, slice_w3;
  
  VIO_Real        edge_u1_v1, edge_u1_v2, edge_u1_v3;
  VIO_Real        edge_u2_v1, edge_u2_v2, edge_u2_v3;
  VIO_Real        edge_u3_v1, edge_u3_v2, edge_u3_v3;
  VIO_Real        edge_u1_w1, edge_u1_w2, edge_u1_w3;
  VIO_Real        edge_u2_w1, edge_u2_w2, edge_u2_w3;
  VIO_Real        edge_u3_w1, edge_u3_w2, edge_u3_w3;
  VIO_Real        edge_v1_w1, edge_v1_w2, edge_v1_w3;
  VIO_Real        edge_v2_w1, edge_v2_w2, edge_v2_w3;
  VIO_Real        edge_v3_w1, edge_v3_w2, edge_v3_w3;
  
  VIO_Real         psum, nsum;
  VIO_Real  *pc;
  VIO_Real        *p11, *p12, *p13, *p14, *p15;
  VIO_Real        *p21, *p22, *p23, *p24, *p25;
  VIO_Real        *p31, *p32, *p33, *p34, *p35;
  VIO_Real        *p41, *p42, *p43, *p44, *p45;
  VIO_Real        *p51, *p52, *p53, *p54, *p55;
  
  
  /* --- 3x3x3 --- */
  
  p11 = r[0][0]; 
  p12 = r[0][1]; 
  p13 = r[0][2]; 
  p21 = r[1][0]; 
  p22 = r[1][1]; 
  p23 = r[1][2]; 
  p31 = r[2][0]; 
  p32 = r[2][1]; 
  p33 = r[2][2]; 
  
  edge_u1_v1 = ( *p11 + *(p11+1) + *(p11+2))/3.0;
  edge_u1_v2 = ( *p12 + *(p12+1) + *(p12+2))/3.0;
  edge_u1_v3 = ( *p13 + *(p13+1) + *(p13+2))/3.0;
  edge_u2_v1 = ( *p21 + *(p21+1) + *(p21+2))/3.0;
  edge_u2_v2 = ( *p22 + *(p22+1) + *(p22+2))/3.0;
  edge_u2_v3 = ( *p23 + *(p23+1) + *(p23+2))/3.0;
  edge_u3_v1 = ( *p31 + *(p31+1) + *(p31+2))/3.0;
  edge_u3_v2 = ( *p32 + *(p32+1) + *(p32+2))/3.0;
  edge_u3_v3 = ( *p33 + *(p33+1) + *(p33+2))/3.0;
  
  edge_u1_w1 = ( *p11 + *p12 + *p13)/3.0;
  edge_u1_w2 = ( *(p11+1) + *(p12+1) + *(p13+1))/3.0;
  edge_u1_w3 = ( *(p11+2) + *(p12+2) + *(p13+2))/3.0;
  edge_u2_w1 = ( *p21 + *p22 + *p23)/3.0;
  edge_u2_w2 = ( *(p21+1) + *(p22+1) + *(p23+1))/3.0;
  edge_u2_w3 = ( *(p21+2) + *(p22+2) + *(p23+2))/3.0;
  edge_u3_w1 = ( *p31 + *p32 + *p33)/3.0;
  edge_u3_w2 = ( *(p31+2) + *(p32+2) + *(p33+2))/3.0;
  edge_u3_w3 = ( *(p31+2) + *(p32+2) + *(p33+2))/3.0;
  
  edge_v1_w1 = ( *p11 + *p21 + *p31)/3.0;
  edge_v1_w2 = ( *(p11+1) + *(p21+1) + *(p31+1))/3.0;
  edge_v1_w3 = ( *(p11+2) + *(p21+2) + *(p31+2))/3.0;
  edge_v2_w1 = ( *p12 + *p22 + *p32)/3.0;
  edge_v2_w2 = ( *(p12+1) + *(p22+1) + *(p32+1))/3.0;
  edge_v2_w3 = ( *(p12+2) + *(p22+2) + *(p32+2))/3.0;
  edge_v3_w1 = ( *p13 + *p23 + *p33)/3.0;
  edge_v3_w2 = ( *(p13+1) + *(p23+1) + *(p33+1))/3.0;
  edge_v3_w3 = ( *(p13+2) + *(p23+2) + *(p33+2))/3.0;
  
  slice_u1 =  (edge_u1_w1 + edge_u1_w2 + edge_u1_w3)/3.0;
  slice_u2 =  (edge_u2_w1 + edge_u2_w2 + edge_u2_w3)/3.0;
  slice_u3 =  (edge_u3_w1 + edge_u3_w2 + edge_u3_w3)/3.0;
  slice_v1 =  (edge_v1_w1 + edge_v1_w2 + edge_v1_w3)/3.0;
  slice_v2 =  (edge_v2_w1 + edge_v2_w2 + edge_v2_w3)/3.0;
  slice_v3 =  (edge_v3_w1 + edge_v3_w2 + edge_v3_w3)/3.0;
  slice_w1 =  (edge_u1_w1 + edge_u2_w1 + edge_u3_w1)/3.0;
  slice_w2 =  (edge_u1_w2 + edge_u2_w2 + edge_u3_w2)/3.0;
  slice_w3 =  (edge_u1_w3 + edge_u2_w3 + edge_u3_w3)/3.0;
  
  c[0] = (slice_u3 - slice_u1) / 2.0;                          /* cu */
  c[1] = (slice_v3 - slice_v1) / 2.0;                          /* cv */
  c[2] = (slice_w3 - slice_w1) / 2.0;                          /* cw */
  c[3] = (slice_u3 + slice_u1 - 2*slice_u2);                   /* cuu */
  c[4] = (slice_v3 + slice_v1 - 2*slice_v2);                   /* cvv */
  c[5] = (slice_w3 + slice_w1 - 2*slice_w2);                   /* cww */
  c[6] = (edge_u3_v3 + edge_u1_v1 - edge_u3_v1 - edge_u1_v3)/4.0;  /* cuv */
  c[7] = (edge_u3_w3 + edge_u1_w1 - edge_u3_w1 - edge_u1_w3)/4.0;  /* cuw */
  c[8] = (edge_v3_w3 + edge_v1_w1 - edge_v3_w1 - edge_v1_w3)/4.0;  /* cvw */
  
 
  /* --- if positive definite (or convex quadratic), then discard --- */
  /* --- recall for correlation, need to min negative similarity  --- */
  /* --- so discard similarity if it gives a minimum.                  --- */

  if (test_max(c)==0) {

    for(i=0; i<9; i++)
      c[i] = 0.0;

  }


} /* derivatives */



/* '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   test_max()

        return 1 if quadratic represented by input is negative definite;
        otherwise return 0.

        strang, pp. 325

''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' */
int  test_max(VIO_Real c[9])

{
  VIO_Real  d1, d2, d3;
  int    test;
  VIO_Real cof11, cof12, cof13;
  VIO_Real cuu, cvv, cww, cuv, cuw, cvw;

  cuu = c[3];
  cvv = c[4];
  cww = c[5];
  cuv = c[6];
  cuw = c[7];
  cvw = c[8];

  d1 = cuu;
  d2 = cuu*cvv - cuv*cuv;
  cof11 = cvv*cww - cvw*cvw;
  cof12 = cvw*cuw - cuv*cww;
  cof13 = cuv*cvw - cvv*cuw;
  d3 =  cof11*cuu + cof12*cuv + cof13*cuw;

  test = 0;
  if ((d1 <= -0.0) && (d2 <= -0.0) && (d3 <=-0.0)) 
    test = 1;
    
  return test;

} /* test_max */


