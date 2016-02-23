/* ----------------------------- MNI Header -----------------------------------
@NAME       : amoeba.c
@DESCRIPTION: Simplex Optimization
@METHOD     : Nelder & Mead 
@GLOBALS    : none
@CALLS      : internal_volume_io library routines/macros
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
@MODIFIED   : Revision 96.10  2006-11-30 09:07:32  rotor
@MODIFIED   :  * many more changes for clean minc 2.0 build
@MODIFIED   :
@MODIFIED   : Revision 96.9  2006/11/29 09:09:33  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 96.8  2005/07/20 20:45:49  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 96.7  2004/02/12 06:08:19  rotor
@MODIFIED   :  * removed /static defs
@MODIFIED   :
@MODIFIED   : Revision 96.6  2004/02/04 20:44:11  lenezet
@MODIFIED   : *** empty log message ***
@MODIFIED   :
@MODIFIED   : Revision 96.5  2002/03/26 14:15:43  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 96.4  2000/03/15 08:42:45  stever
@MODIFIED   : Code cleanup: all functions prototyped (except ParseArgs.c), no useless declarations, etc
@MODIFIED   :
@MODIFIED   : Revision 96.3  1999/10/25 19:59:06  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :
 * Revision 96.2  1997/11/03  20:05:41  louis
 * no changes
 *
 * Revision 96.2  1997/11/03  20:05:41  louis
 * no changes
 *
 * Revision 96.1  1997/11/03  15:04:23  louis
 * working version, before creation of mni_animal package, and before inserting
 * distance transforms
 *
 * Revision 96.1  1997/11/03  15:04:23  louis
 * working version, before creation of mni_animal package, and before inserting
 * distance transforms
 *
 * Revision 96.0  1996/08/21  18:22:10  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:22:01  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:16:03  louis
 * Never released version 0.95
 *
 * Revision 1.4  1996/08/12  14:15:52  louis
 * Pre-release
 *
 * Revision 1.3  1996/03/25  10:33:15  collins
 * changed perform_amoeba() so that it counts the actual number of
 * objective function evaluations instead of the number of calls to
 * perform_amoeba().  (the previous version was underestimating the
 * number of function calls by a factor of 1.8-2.0!)
 *
 * Revision 1.2  1995/09/14  09:56:02  collins
 * removed the call to numerically_close and replaced it with the same
 * terminating condition as used in the previous version of amoeba, so
 * that the same termination tolerance can be used.
 *
 * Also removed the warning message "No improvement after X steps..."
 *
 * Revision 1.2  1995/09/14  09:56:02  collins
 * removed the call to numerically_close and replaced it with the same
 * terminating condition as used in the previous version of amoeba, so
 * that the same termination tolerance can be used.
 *
 * Also removed the warning message "No improvement after X steps..."
 *
 * Revision 1.1  1995/09/07  10:05:11  collins
 * Initial revision
 *

   LC: this version of David's ameoba.c  is slightly modified from that
   found in BIC_PL, in that the function numerically_close is included
   here, and bic_pl is not needed for minctracc.   This way, only 
   volume_io.h and amoeba.h are needed for inclusion.

   there are four  functions that the user has to call (in order):
   initialize_amoeba(), 
   perform_amoeba(),
   get_amoeba_parameters(),
   terminate_amoeba()


   1: initialize_amoeba is called first to set up the simplex optimization
      process.  It requires the follwing parameters:
   initialize_amoeba(
      amoeba_struct     *amoeba,                 the structure, allocated by
                                                 the calling program
      int               n_parameters,            the number of free parameter
                                                 to be optimized 
      VIO_Real              initial_parameters[],    the array containing the
                                                 starting guess
      VIO_Real              parameter_delta,         the radius of the simplex
      amoeba_function   function,                the name of the objective
                                                 function
      void              *function_data,          pointer to any data you
                                                 want to pass to the objective
                                                 function
      VIO_Real              tolerance )              the stopping tolerance
    
   2: perform_amoeba(amoeba) is then called to actually do the optimization.
      The function returns TRUE if the optimization is successful, otherwise
      it returns FALSE.

   3: get_amoeba_parameters(amoeba, parameters) must be called to extract
      the optimized parameters.  'parameters' is a VIO_Real array, allocated by
      the calling program.

   4: terminate_amoeba(amoeba) is called in the end to free up internal data
      structures used by the procedure.  This procedure should be called after
      get_amoeba_parameters, since the optimized parameters will be freed (and
      lost) in the call to terminate_amoeba().

---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/minctracc/Optimize/amoeba.c,v 96.10 2006-11-30 09:07:32 rotor Exp $";
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif //HAVE_CONFIG_H


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

#if 0
/** call to this function is commented out, so we don't need to compile it **/
static  VIO_BOOL  numerically_close(
    VIO_Real  n1,
    VIO_Real  n2,
    VIO_Real  threshold_ratio )
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
#endif


/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_function_value
@INPUT      : amoeba
              parameters
@OUTPUT     : 
@RETURNS    : function value
@DESCRIPTION: Evaluates the function being minimized by amoeba, by calling
              the user function.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */
static  VIO_Real  get_function_value(
    amoeba_struct  *amoeba,
    float          parameters[] )
{
    return( (*amoeba->function) ( amoeba->function_data, parameters ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_amoeba
@INPUT      : n_parameters
              initial_parameters
              parameter_deltas
              function
              function_data
              tolerance
@OUTPUT     : amoeba
@RETURNS    : 
@DESCRIPTION: Initializes the amoeba structure to minimize the function.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */
 void  initialize_amoeba(
    amoeba_struct     *amoeba,
    int               n_parameters,
    VIO_Real              initial_parameters[],
    VIO_Real              parameter_delta,
    amoeba_function   function,
    void              *function_data,
    VIO_Real              tolerance )
{
    int    i, j;
    
#ifndef MNI_AUTOREG_OLD_AMOEBA_INIT
    float *parameter_fwd;
    float *parameter_bwd;
    VIO_Real cost_fwd;
    VIO_Real cost_bwd;
#endif 

    amoeba->n_parameters = n_parameters;
    amoeba->function = function;
    amoeba->function_data = function_data;

    amoeba->tolerance = tolerance;
    amoeba->n_steps_no_improvement = 0;
    VIO_ALLOC2D( amoeba->parameters, n_parameters+1, n_parameters );
    ALLOC( amoeba->values, n_parameters+1 );

    ALLOC( amoeba->sum, n_parameters );
    
#ifndef MNI_AUTOREG_OLD_AMOEBA_INIT
    ALLOC( parameter_fwd, n_parameters);
    ALLOC( parameter_bwd, n_parameters);
#endif

    for(j=0; j<n_parameters; j++)
        amoeba->sum[j] = 0.0;
    
    for(i=0; i<n_parameters+1; i++)
    {
        for(j=0; j<n_parameters; j++)
        {
#ifdef MNI_AUTOREG_OLD_AMOEBA_INIT
            amoeba->parameters[i][j] = (float) initial_parameters[j];
#else
            parameter_fwd[j] = (float) initial_parameters[j];
            parameter_bwd[j] = (float) initial_parameters[j];
#endif
            if( i > 0 && j == i - 1 )
#ifdef MNI_AUTOREG_OLD_AMOEBA_INIT
                amoeba->parameters[i][j] += parameter_delta;
            amoeba->sum[j] += amoeba->parameters[i][j];
#else
            {
              parameter_fwd[j] += parameter_delta;
              parameter_bwd[j] -= parameter_delta;
            }
#endif
        }

#ifdef MNI_AUTOREG_OLD_AMOEBA_INIT
      amoeba->values[i] = get_function_value( amoeba, amoeba->parameters[i] );
#else
      // Test the parameters 
      cost_fwd = get_function_value( amoeba, parameter_fwd);
      cost_bwd = get_function_value( amoeba, parameter_bwd);

      // Use the one with the lowest value
      for(j=0; j<n_parameters; j++)
      {
          amoeba->parameters[i][j] = ( cost_fwd < cost_bwd ) ? parameter_fwd[j] : parameter_bwd[j];
          amoeba->sum[j] += amoeba->parameters[i][j];
      }
      amoeba->values[i] = ( cost_fwd < cost_bwd ) ? cost_fwd : cost_bwd;
#endif
    }
    
#ifndef MNI_AUTOREG_OLD_AMOEBA_INIT
    FREE( parameter_fwd );
    FREE( parameter_bwd );
#endif    
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_amoeba_parameters
@INPUT      : amoeba
@OUTPUT     : parameters
@RETURNS    : function value
@DESCRIPTION: Passes back the current position of the amoeba (best value),
              and returns the function value at that point.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */
 VIO_Real  get_amoeba_parameters(
    amoeba_struct  *amoeba,
    VIO_Real           parameters[] )
{
    int   i, j, low;

    low = 0;
    for(i=1; i<amoeba->n_parameters+1; i++)
    {
      
        if( amoeba->values[i] < amoeba->values[low] )
            low = i;
    }
    for(j=0; j<amoeba->n_parameters; j++)
        parameters[j] = (VIO_Real) amoeba->parameters[low][j];

    return( amoeba->values[low] );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : terminate_amoeba
@INPUT      : amoeba
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Frees the amoeba.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

 void  terminate_amoeba(
    amoeba_struct  *amoeba )
{
    VIO_FREE2D( amoeba->parameters );
    FREE( amoeba->values );
    FREE( amoeba->sum );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : try_amoeba
@INPUT      : amoeba
              sum
              high
              fac
@OUTPUT     : 
@RETURNS    : value
@DESCRIPTION: Does a modification to the high vertex of the amoeba and
              returns the value of the new point.  If the new point is
              better (smaller value), it replaces the high vertex of the
              amoeba.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */
static  VIO_Real  try_amoeba(
    amoeba_struct  *amoeba,
    VIO_Real       sum[],
    int            high,
    VIO_Real       fac )
{
    int    j;
    VIO_Real   y_try, fac1, fac2;
    float  *parameters;

    ALLOC( parameters, amoeba->n_parameters );

    fac1 = (1.0 - fac) / amoeba->n_parameters;
    fac2 = fac - fac1;

    for(j=0; j<amoeba->n_parameters; j++)
        parameters[j] = sum[j] * fac1 + amoeba->parameters[high][j] * fac2;


    y_try = get_function_value( amoeba, parameters );

    if( y_try < amoeba->values[high] )
    {
        amoeba->values[high] = y_try;
        for(j=0; j<amoeba->n_parameters; j++)
        {
            sum[j] += parameters[j] - amoeba->parameters[high][j];
            amoeba->parameters[high][j] = parameters[j];
        }
    }

    FREE( parameters );

    return( y_try );
}

#define  N_STEPS_NO_IMPROVEMENT  6

/* ----------------------------- MNI Header -----------------------------------
@NAME       : perform_amoeba
@INPUT      : amoeba
@OUTPUT     : 
@RETURNS    : TRUE if numerically significant improvement
@DESCRIPTION: Performs one iteration of an amoeba, returning true if a
              numerically significant improvement has been found recently.
              Even if it returns FALSE, you can keep calling this function,
              since it may be contracting with no improvement, but will 
              eventually shrink small enough to get an improvment.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */
 VIO_BOOL  perform_amoeba(
    amoeba_struct  *amoeba, int *num_funks )
{
    int     i, j, low, high, next_high;
    VIO_Real    y_try, y_save;
    VIO_BOOL  improvement_found;
    VIO_Real tol;

    improvement_found = TRUE;

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

    for(i=2; i<amoeba->n_parameters+1; i++)
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

/*    if( numerically_close( amoeba->values[low], amoeba->values[high],
                           amoeba->tolerance ) )
 */

    tol = 2.0 * fabs(amoeba->values[high]-amoeba->values[low]) /
      (fabs(amoeba->values[high]) + fabs(amoeba->values[low]));

    if (tol < amoeba->tolerance)
   {
        ++amoeba->n_steps_no_improvement;
        if( ++amoeba->n_steps_no_improvement == N_STEPS_NO_IMPROVEMENT ) {
          return( FALSE );
        }
    }
    else
        amoeba->n_steps_no_improvement = 0;

    y_try = try_amoeba( amoeba, amoeba->sum, high, -FLIP_RATIO );
    (*num_funks)++;

    if( y_try <= amoeba->values[low] ) {
        y_try = try_amoeba( amoeba, amoeba->sum, high, STRETCH_RATIO );
        (*num_funks)++;
      }
    else if( y_try >= amoeba->values[next_high] )
    {
        y_save = amoeba->values[high];
        y_try = try_amoeba( amoeba, amoeba->sum, high, CONTRACT_RATIO );
        (*num_funks)++;
        
        if( y_try >= y_save )
        {
            for(i=0; i<amoeba->n_parameters+1; i++)
            {
                if( i != low )
                {
                    for(j=0; j<amoeba->n_parameters; j++)
                    {
                        amoeba->parameters[i][j] = (amoeba->parameters[i][j] +
                                            amoeba->parameters[low][j]) / 2.0;
                    }

                    amoeba->values[i] = get_function_value( amoeba,
                                                  amoeba->parameters[i] );
                    (*num_funks)++;
                }
            }

            for(j=0; j<amoeba->n_parameters; j++)
            {
                amoeba->sum[j] = 0.0;
                for(i=0; i<amoeba->n_parameters+1; i++)
                    amoeba->sum[j] += amoeba->parameters[i][j];
            }
        }
    }

    return( improvement_found );
}


