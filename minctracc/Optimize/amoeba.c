/* ----------------------------- MNI Header -----------------------------------
@NAME       : amoeba.c
@DESCRIPTION: Simplex Optimization
@METHOD     : Nelder & Mead 
@GLOBALS    : none
@CALLS      : volume_io library routines/macros
@COPYRIGHT  :
              Copyright 1995 David MacDonald, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.

@CREATED    : 
@MODIFIED   : $Log: amoeba.c,v $
@MODIFIED   : Revision 1.1  1995-09-07 10:05:11  louis
@MODIFIED   : Initial revision
@MODIFIED   :

   LC: this version of David's ameoba.c  is slightly modified from that
   found in BIC_PL, in that the function numerically_close is included
   here, and bic_pl is not needed for minctracc.

---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Optimize/amoeba.c,v 1.1 1995-09-07 10:05:11 louis Exp $";
#endif


#include <volume_io.h>
#include <amoeba.h>

#define  FLIP_RATIO      1.0
#define  CONTRACT_RATIO  0.5
#define  STRETCH_RATIO   2.0

/* ----------------------------- MNI Header -----------------------------------
@NAME       : numerically_close
@INPUT      : n1
              n2
              threshold_ratio
@OUTPUT     : 
@RETURNS    : TRUE if the numbers are within the threshold ratio
@DESCRIPTION: Decides if two numbers are close to each other.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

private  BOOLEAN  numerically_close(
    Real  n1,
    Real  n2,
    Real  threshold_ratio )
{

    double  avg, diff;

    diff = n1 - n2;
    if( diff < 0.0 )  diff = -diff;

    if( n1 <= threshold_ratio && n1 >= -threshold_ratio &&
        n2 <= threshold_ratio && n2 >= -threshold_ratio )
    {
        return( diff <= threshold_ratio );
    }

    avg = (n1 + n2) / 2.0;

    if( avg == 0.0 )
        return( diff <= (double) threshold_ratio );

    if( avg < 0.0 )  avg = -avg;

    return( (diff / avg) <= (double) threshold_ratio );
}


private  Real  get_function_value(
    amoeba_struct  *amoeba,
    float          parameters[] )
{
    return( (*amoeba->function) ( amoeba->function_data, parameters ) );
}

public  void  initialize_amoeba(
    amoeba_struct     *amoeba,
    int               n_parameters,
    Real              initial_parameters[],
    Real              parameter_delta,
    amoeba_function   function,
    void              *function_data,
    Real              tolerance )
{
    int    i, j;

    amoeba->n_parameters = n_parameters;
    amoeba->function = function;
    amoeba->function_data = function_data;
    amoeba->tolerance = tolerance;
    amoeba->n_steps_no_improvement = 0;
    ALLOC2D( amoeba->parameters, n_parameters+1, n_parameters );
    ALLOC( amoeba->values, n_parameters+1 );

    ALLOC( amoeba->sum, n_parameters );

    for_less( j, 0, n_parameters )
        amoeba->sum[j] = 0.0;

    for_less( i, 0, n_parameters+1 )
    {
        for_less( j, 0, n_parameters )
        {
            amoeba->parameters[i][j] = (float) initial_parameters[j];
            if( i > 0 && j == i - 1 )
                amoeba->parameters[i][j] += parameter_delta;
            amoeba->sum[j] += amoeba->parameters[i][j];
        }

        amoeba->values[i] = get_function_value( amoeba, amoeba->parameters[i] );
    }
}

public  Real  get_amoeba_parameters(
    amoeba_struct  *amoeba,
    Real           parameters[] )
{
    int   i, j, low;

    low = 0;
    for_less( i, 1, amoeba->n_parameters+1 )
    {
        if( amoeba->values[i] < amoeba->values[low] )
            low = i;
    }

    for_less( j, 0, amoeba->n_parameters )
        parameters[j] = (Real) amoeba->parameters[low][j];

    return( amoeba->values[low] );
}

public  void  terminate_amoeba(
    amoeba_struct  *amoeba )
{
    FREE2D( amoeba->parameters );
    FREE( amoeba->values );
    FREE( amoeba->sum );
}

private  Real  try_amoeba(
    amoeba_struct  *amoeba,
    Real           sum[],
    int            high,
    Real           fac )
{
    int    j;
    Real   y_try, fac1, fac2;
    float  *parameters;

    ALLOC( parameters, amoeba->n_parameters );

    fac1 = (1.0 - fac) / amoeba->n_parameters;
    fac2 = fac - fac1;

    for_less( j, 0, amoeba->n_parameters )
        parameters[j] = sum[j] * fac1 + amoeba->parameters[high][j] * fac2;

    y_try = get_function_value( amoeba, parameters );

    if( y_try < amoeba->values[high] )
    {
        amoeba->values[high] = y_try;
        for_less( j, 0, amoeba->n_parameters )
        {
            sum[j] += parameters[j] - amoeba->parameters[high][j];
            amoeba->parameters[high][j] = parameters[j];
        }
    }

    FREE( parameters );

    return( y_try );
}

#define  N_STEPS_NO_IMPROVEMENT  10

public  BOOLEAN  perform_amoeba(
    amoeba_struct  *amoeba )
{
    int     i, j, low, high, next_high;
    Real    y_try, y_save;

    if( amoeba->values[0] > amoeba->values[1] )
    {
        high = 0;
        next_high = 1;
    }
    else
    {
        high = 1;
        next_high = 0;
    }

    low = next_high;

    for_less( i, 2, amoeba->n_parameters+1 )
    {
        if( amoeba->values[i] < amoeba->values[low] )
            low = i;
        else if( amoeba->values[i] > amoeba->values[high] )
        {
            next_high = high;
            high = i;
        }
        else if( amoeba->values[i] > amoeba->values[next_high] )
            next_high = i;
    }

    if( numerically_close( amoeba->values[low], amoeba->values[high],
                           amoeba->tolerance ) )
    {
        ++amoeba->n_steps_no_improvement;
        if( ++amoeba->n_steps_no_improvement == N_STEPS_NO_IMPROVEMENT )
            return( FALSE );
    }
    else
        amoeba->n_steps_no_improvement = 0;

    y_try = try_amoeba( amoeba, amoeba->sum, high, -FLIP_RATIO );

    if( y_try <= amoeba->values[low] )
        y_try = try_amoeba( amoeba, amoeba->sum, high, STRETCH_RATIO );
    else if( y_try >= amoeba->values[next_high] )
    {
        y_save = amoeba->values[high];
        y_try = try_amoeba( amoeba, amoeba->sum, high, CONTRACT_RATIO );

        if( y_try >= y_save )
        {
            for_less( i, 0, amoeba->n_parameters+1 )
            {
                if( i != low )
                {
                    for_less( j, 0, amoeba->n_parameters )
                    {
                        amoeba->parameters[i][j] = (amoeba->parameters[i][j] +
                                            amoeba->parameters[low][j]) / 2.0;
                    }

                    amoeba->values[i] = get_function_value( amoeba,
                                                  amoeba->parameters[i] );
                }
            }

            for_less( j, 0, amoeba->n_parameters )
            {
                amoeba->sum[j] = 0.0;
                for_less( i, 0, amoeba->n_parameters+1 )
                    amoeba->sum[j] += amoeba->parameters[i][j];
            }
        }
    }

    return( TRUE );
}
