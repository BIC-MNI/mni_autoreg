#ifndef MINCTRACC_POINT_VECTOR_H
#define MINCTRACC_POINT_VECTOR_H

				/* redefine Point and Vector */

#define N_DIMENSIONS 3

typedef  struct
{
    Real   coords[N_DIMENSIONS];
} PointR;

typedef  struct
{
    Real   coords[N_DIMENSIONS];
} VectorR;


#endif
