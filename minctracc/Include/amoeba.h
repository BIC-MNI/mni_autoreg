#ifndef  _DEF_AMOEBA_H
#define  _DEF_AMOEBA_H

#include  <volume_io.h>

typedef  Real    (*amoeba_function) ( void *, float [] );

typedef  struct
{
    int               n_parameters;
    float             **parameters;
    Real              *values;
    amoeba_function   function;
    void              *function_data;
    Real              tolerance;
    Real              *sum;
    int               n_steps_no_improvement;
} amoeba_struct;

#endif
