#ifndef MINCTRACC_POINT_VECTOR_H
#define MINCTRACC_POINT_VECTOR_H

                                /* redefine Point and Vector */

#define VIO_N_DIMENSIONS 3

typedef  struct
{
    VIO_Real   coords[VIO_N_DIMENSIONS];
} PointR;

typedef  struct
{
    VIO_Real   coords[VIO_N_DIMENSIONS];
} VectorR;


#endif
