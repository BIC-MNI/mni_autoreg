#include <float.h>

#define MY_FLT_EPSILON  0.001             /* 10*FLT_EPSILON*/

#define LESS( a, b )     ( ((a) + MY_FLT_EPSILON) < (b) )
#define GREATER( a, b )  ( ((a) - MY_FLT_EPSILON) > (b) )

#define LESS_EQ( a, b )     ( ((a) - MY_FLT_EPSILON) < (b) )
#define GREATER_EQ( a, b )  ( ((a) + MY_FLT_EPSILON) > (b) )

#define EQUAL( a, b )    ( ABS( (a) - (b) ) < MY_FLT_EPSILON )

#define INSIDE( a, b, c,  a1, b1, c1, a2, b2, c2 ) \
   (GREATER_EQ ( (a) , (a1) ) && LESS_EQ ( (a) , (a2) ) &&  \
    GREATER_EQ ( (b) , (b1) ) && LESS_EQ ( (b) , (b2) ) && \
    GREATER_EQ ( (c) , (c1) ) && LESS_EQ ( (c) , (c2) )) 

#define INSIDE_W( a, b, c, d ) \
   (GREATER_EQ ( (a) , (d)->Xmin ) && LESS_EQ ( (a) , (d)->Xmax ) &&  \
    GREATER_EQ ( (b) , (d)->Ymin ) && LESS_EQ ( (b) , (d)->Ymax ) && \
    GREATER_EQ ( (c) , (d)->Zmin ) && LESS_EQ ( (c) , (d)->Zmax )) 

#define INSIDE_V( a, b, c, d ) \
   (GREATER_EQ ( (a) , -0.5 ) && LESS ( (a) , (d)->cols-0.5 ) &&  \
    GREATER_EQ ( (b) , -0.5 ) && LESS ( (b) , (d)->rows-0.5 ) && \
    GREATER_EQ ( (c) , -0.5 ) && LESS ( (c) , (d)->slices-0.5 )) 

