#include <volume_io.h>
#include "minctracc_minctracc_arg_data.h"


void initializeArgs(Arg_Data *args);

VIO_General_transform* minctracc( VIO_Volume source, 
                                  VIO_Volume target, 
                                  VIO_Volume sourceMask, 
                                  VIO_Volume targetMask, 
                                  VIO_General_transform *initialXFM, 
                                  int iterations, 
                                  float weight, 
                                  float simplexSize, 
                                  float stiffness, 
                                  float similarity, 
                                  float sub_lattice, 
                                  Arg_Data *args);

