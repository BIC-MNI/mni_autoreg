#ifndef MINCTRACC_INIT_LATTICE_H
#define MINCTRACC_INIT_LATTICE_H

#include "point_vector.h"

void get_volume_XYZV_indices(VIO_Volume data, int xyzv[]);

/* 
   find the largest lattice array, aligned with the existing data volume,
   that will fit inside the volumetric space define by the data volume,
   given the user spacing.
   
   return the lattice info in start, count, step and directions. 

   count is constrained to be positive,
   step is constrained to have the same sign as the data volume.

   if one of the dimensions of the input volume is of length 1, then the
   output lattice will be of the same length.

*/

void set_up_lattice(VIO_Volume data,       /* in: volume  */
                           double *user_step, /* in: user requested spacing for lattice */
                           double *start,     /* out:world starting position of lattice */
                           int    *count,     /* out:number of steps in each direction */
                           double *step,      /* out:step size in each direction */
                           VectorR directions[]); /* out: vector directions for each index*/


/* 
   in this procedure, the smallest (in number of samples) 3D lattice
   is defined that complete covers either d1 or d2 (taking into account
   the mask volumes m1 and m2).  The lattice spacing is defined by the 
   globals->step variable, where the steps are stored in X,Y,Z order.

   The procedure returns the globals->start, globals->count and 
   globals->step that define the lattice (each in X Y Z order).  

   globals->smallest_vol ==1 or ==2, indicating on which volume space
   the lattice is defined.  (This volume always ==2 when NONLIN 
   transformations are optimized.

*/


void init_lattice(VIO_Volume d1,
                         VIO_Volume d2,
                         VIO_Volume m1,
                         VIO_Volume m2, 
                         Arg_Data *globals);





#endif
