#ifndef GLOBALS_H
#define GLOBALS_H

#include <stdlib.h>
#include <float.h>
#include <config.h>
#include <volume_io.h>
#include <ParseArgv.h>
#include <minctracc.h>


/*  ------------------------ global variables  ------------------------ */

extern char  *prog_name;

extern VIO_Volume  mask_data;
extern VIO_Volume  mask_model;
extern VIO_Volume  model;
extern VIO_Volume  model_dx;
extern VIO_Volume  model_dy;
extern VIO_Volume  model_dz;
extern VIO_Volume  model_dxyz;
extern VIO_Volume  data;
extern VIO_Volume  data_dx;
extern VIO_Volume  data_dy;
extern VIO_Volume  data_dz;
extern VIO_Volume  data_dxyz;

extern double  ftol;
extern double  simplex_size;
extern int     iteration_limit;
extern double  iteration_weight;
extern double  smoothing_weight;
extern double  similarity_cost_ratio;
extern int     number_dimensions;
extern int     Matlab_num_steps;
extern int     Diameter_of_local_lattice;

extern int     invert_mapping_flag;
extern int     clobber_flag;

extern VIO_Real initial_corr, final_corr;


  
/*------------------------    Command line arguments --------------------*/


extern ArgvInfo argTable[];


#endif /* GLOBALS_H */


