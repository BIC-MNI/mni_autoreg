/* ----------------------------- MNI Header -----------------------------------
@NAME       : optimize.c
@DESCRIPTION: collection of routines for user-specified optimization.
@METHOD     : now, only simplex method is used.
@MODIFIED   : $Log: optimize.c,v $
@MODIFIED   : Revision 1.5  1993-11-15 16:27:08  louis
@MODIFIED   : working version, with new library, with RCS revision stuff,
@MODIFIED   : before deformations included
@MODIFIED   :
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Optimize/optimize.c,v 1.5 1993-11-15 16:27:08 louis Exp $";
#endif

#include <mni.h>
#include <recipes.h>
#include <limits.h>

#include "constants.h"
#include "arg_data.h"
#include "objectives.h"
#include "make_rots.h"
#include "segment_table.h"

extern Arg_Data main_args;

/* external calls: */

public void make_zscore_volume(Volume d1, Volume m1, 
			       float threshold); 

public void add_speckle_to_volume(Volume d1, 
				  float speckle,
				  double  *start, int *count, VectorR directions[]);

         Volume   Gdata1, Gdata2, Gmask1, Gmask2;
         int      Ginverse_mapping_flag, Gndim;

extern   double   ftol ;        
extern   double   simplex_size ;

         Segment_Table  *segment_table;

/* ----------------------------- MNI Header -----------------------------------
@NAME       : fit_function
@INPUT      : params - a variable length array of floats
@OUTPUT     :               
@RETURNS    : a float value of the user requested objective function,
              measuring the similarity between two data sets.
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */

public float fit_function(float *params) 
{

  Transform *mat;
  int i;
  float r;

  double trans[3];
  double cent[3];
  double rots[3];
  double scale[3];
  double shear[3];

  for_less( i, 0, 3 )		/* set translations */
    trans[i] = params[i+1]; 

  for_less( i, 0, 3 )		/* set rotations */
    rots[i] = params[i+4]; 

  if (Gndim >= 7)  {
    scale[0] = params[7];
    if (Gndim >7) {		/* three scales */
      scale[1] = params[8];
      scale[2] = params[9];
    }
    else {			/* only one scale */
      scale[1] = params[7];
      scale[2] = params[7];
    }
  }
  else {			/* fixed scale */
    scale[0] = 1.0;
    scale[1] = 1.0;
    scale[2] = 1.0;
  }

  for_less( i, 0, 3 )		/* set shears */
    shear[i] = 0.0;
  
  if (Gndim==10) {
      shear[0] = params[10]; 
  }
  else if (Gndim==12) {
    for_less( i, 0, 3 )		
      shear[i] = params[10+i]; 
  }

  for_less( i, 0, 3 )
    cent[i] = main_args.trans_info.center[i]; /* GLOBAL MAIN_ARGS USED HERE */
  

  mat = get_linear_transform_ptr(main_args.trans_info.transformation);

  if (Ginverse_mapping_flag)
    build_inverse_transformation_matrix(mat, cent, trans, scale, shear, rots);
  else
    build_transformation_matrix(mat, cent, trans, scale, shear, rots);

				/* call the needed objective function */
  r = (main_args.obj_function)(Gdata1,Gdata2,Gmask1,Gmask2,&main_args);

  return(r);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : optimize_simplex
                get the parameters necessary to map volume 1 to volume 2
		using the simplex optimizaqtion algorithm and a user specified
		objective function.
@INPUT      : d1,d2:
                two volumes of data (already in memory).
	      m1,m2:
                two mask volumes for data (already in memory).
	      globals:
	        a global data structure containing info from the command line,
		including the input parameters to be optimized, the input matrix,
		and a plethora of flags!
@OUTPUT     : 
@RETURNS    : TRUE if ok, FALSE if error.
@DESCRIPTION: 
@METHOD     : uses the ameoba simplex algorithm from numerical recipes. 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Fri Jun 11 11:16:25 EST 1993 LC
@MODIFIED   : 
---------------------------------------------------------------------------- */
public BOOLEAN optimize_simplex(Volume d1,
				Volume d2,
				Volume m1,
				Volume m2, 
				Arg_Data *globals)
{
  BOOLEAN 
    stat;
  float 
    local_ftol,
    **p, *y;
  int 
    i,j, 
    ndim, nfunk;

  Transform
    *mat;

  double trans[3];
  double cent[3];
  double rots[3];
  double scale[3];
  double shear[3];

  
  stat = TRUE;
  local_ftol = ftol;

  switch (globals->trans_info.transform_type) {
  case TRANS_PROCRUSTES: 
    ndim = 7;
    break;
  case TRANS_LSQ: 
    ndim = 12;
    break;
  case TRANS_LSQ6: 
    ndim = 6;
    break;
  case TRANS_LSQ7: 
    ndim = 7;
    break;
  case TRANS_LSQ9: 
    ndim = 9;
    break;
  case TRANS_LSQ12: 
    ndim = 12;
    break;
  default:
    (void)fprintf(stderr, "Unknown type of transformation requested (%d)\n",
		   globals->trans_info.transform_type);
    (void)fprintf(stderr, "Error in line %d, file %s\n",__LINE__, __FILE__);
    stat = FALSE;
  }


				/* set GLOBALS to communicate with the
				   function to be fitted!              */
  if (stat) {
    Gndim = ndim;
    if (globals->smallest_vol == 1) {
      Gdata1 = d1;
      Gdata2 = d2;
      Gmask1 = m1;
      Gmask2 = m2;
      Ginverse_mapping_flag = FALSE;
    }
    else {
      Gdata1 = d2;
      Gdata2 = d1;
      Gmask1 = m2;
      Gmask2 = m1;
      Ginverse_mapping_flag = TRUE;
    }
  }


  if (stat) {

    p = matrix(1,ndim+1,1,ndim); /* simplex */
    y = vector(1,ndim+1);        /* value of correlation at simplex vertices */
    
    p[1][1]=globals->trans_info.translations[0];
    p[1][2]=globals->trans_info.translations[1];
    p[1][3]=globals->trans_info.translations[2];
    
    p[1][4]=globals->trans_info.rotations[0]; 
    p[1][5]=globals->trans_info.rotations[1]; 
    p[1][6]=globals->trans_info.rotations[2];
    
    if (ndim >= 7) p[1][7]=globals->trans_info.scales[0];
    if (ndim >7) {
      p[1][8]=globals->trans_info.scales[1];
      p[1][9]=globals->trans_info.scales[2];
    }
    
    if (ndim==12) {
      for_less( i, 0, 3 )		/* set shears */
	p[1][10+i] = globals->trans_info.shears[i];
    }


    for (i=2; i<=(ndim+1); ++i)	/* copy initial guess to all points of simplex */
      for (j=1; j<=ndim; ++j)
	p[i][j] = p[1][j];
    
    p[2][1]=p[1][1]+simplex_size;		/* set up all vertices of simplex */
    p[3][2]=p[1][2]+simplex_size;
    p[4][3]=p[1][3]+simplex_size/5;
    
    p[5][4]=p[1][4] + (simplex_size*DEG_TO_RAD);
    p[6][5]=p[1][5] + (simplex_size*DEG_TO_RAD);
    p[7][6]=p[1][6] + (simplex_size*DEG_TO_RAD);
    
    if (ndim >= 7) p[8][7]=p[1][7] + simplex_size/50;	
    if (ndim >7) {
      p[9][8]=p[1][8] + simplex_size/50;
      p[10][9]=p[1][9]+ simplex_size/50;
    }

    if (ndim==10) {
	p[11][10]=p[1][10] + (simplex_size*DEG_TO_RAD);
    } else
      if (ndim==12) {
	for_less( i, 0, 3 )
	  p[11+i][10+i]=p[1][10+i] + (simplex_size*DEG_TO_RAD);
      }


    for (i=1; i<=(ndim+1); ++i)	{   /* set up value of correlation at all points of simplex */
      
      y[i] = fit_function(p[i]);

      (void)print ("corr = %6.4f",y[i]);
      for (j=1; j<=3; ++j)  {
	(void)print (" %7.2f",p[i][j]);
      }
      for (j=4; j<=6; ++j)  {
	(void)print (" %7.3f",p[i][j]*RAD_TO_DEG);
      }
      for (j=7; j<=ndim; ++j)  {
	(void)print (" %6.4f",p[i][j]);
      }
      (void)print ("\n");
      
    }
	
    amoeba(p,y,ndim,local_ftol,fit_function,&nfunk);
    
    (void)print("after %d iterations\n",nfunk);
    
    for (i=1; i<=(ndim+1); ++i)	{   /* print out value of correlation at all points of simplex */
      if (i==1) {
	(void)print ("end corr = %f",y[i]);
	for (j=1; j<=3; ++j)  {
	  (void)print (" %f",p[i][j]);
	}
	for (j=4; j<=6; ++j)  {
	  (void)print (" %f",p[i][j]*RAD_TO_DEG);
	}
	for (j=7; j<=ndim; ++j)  {
	  (void)print (" %f",p[i][j]);
	}
	(void)print ("\n");
      }
    }   
    
				/* copy result into main data structure */
    for_less ( i, 0, 3 ) {
      globals->trans_info.translations[i] = p[1][i+1];
      globals->trans_info.rotations[i]    = p[1][i+4];
    }

    if (ndim==6) {
      for_less ( i, 0, 3 ) 
	globals->trans_info.scales[i]       =  1.0;
    } else
      if (ndim==7) {
	for_less ( i, 0, 3 ) 
	  globals->trans_info.scales[i]       =  p[1][7];
      } else {
	for_less ( i, 0, 3 ) 
	  globals->trans_info.scales[i]       =  p[1][i+7];
      } 
    
    if (ndim==12) {
      for_less( i, 0, 3 )		/* get shears */
	globals->trans_info.shears[i] = p[1][10+i]; 
    }
    else
      for_less( i, 0, 3 )		/* set shears */
	globals->trans_info.shears[i] = 0.0;
    
    

    for_less( i, 0, 3 ) {		/* set translations */
      trans[i] = globals->trans_info.translations[i]; 
      rots[i]  = globals->trans_info.rotations[i];
      scale[i] = globals->trans_info.scales[i];
      shear[i] = globals->trans_info.shears[i];
      cent[i]  = globals->trans_info.center[i];
    }

    mat = get_linear_transform_ptr(globals->trans_info.transformation);
    
    build_transformation_matrix(mat, cent, trans, scale, shear, rots);

    free_matrix(p,1,ndim+1,1,ndim); /* simplex */
    free_vector(y,1,ndim+1);        /* value of correlation at simplex vertices */

  }

  return( stat );
}



/* ----------------------------- MNI Header -----------------------------------
@NAME       : optimize_linear_transformation
                get the parameters necessary to map volume 1 to volume 2
		using a user specified optimization strategy and objective
		function.
@INPUT      : d1,d2:
                two volumes of data (already in memory).
	      m1,m2:
                two mask volumes for data (already in memory).
	      globals:
	        a global data structure containing info from the command line,
		including the input parameters to be optimized, the input matrix,
		and a plethora of flags!
@OUTPUT     : 
@RETURNS    : TRUE if ok, FALSE if error.
@DESCRIPTION: 
@METHOD     :
                1- this routine begins by initializing the volume data structures
		to be used by the objective functions.

		2- optimization function is called
		
		3- the optimized parameters are returned in globals...
@GLOBALS    : 
@CALLS      : 
@CREATED    : Wed Jun  9 12:56:08 EST 1993 LC
@MODIFIED   : 
---------------------------------------------------------------------------- */
public BOOLEAN optimize_linear_transformation(Volume d1,
					      Volume d2,
					      Volume m1,
					      Volume m2, 
					      Arg_Data *globals)
{
  BOOLEAN 
    stat;
  int i;

  stat = TRUE;
  
	     /*----------------- prepare data for optimization ------------ */
  
  if (globals->obj_function == zscore_objective) { /* replace volume d1 and d2 by zscore volume  */
    make_zscore_volume(d1,m1,(float)globals->threshold[0]);
    make_zscore_volume(d2,m2,(float)globals->threshold[1]);
  } 
  else  if (globals->obj_function == ssc_objective) {	/* add speckle to the data set */

    make_zscore_volume(d1,m1,(float)globals->threshold[0]); /* need to make data sets comparable */
    make_zscore_volume(d2,m2,(float)globals->threshold[1]); /* in mean and sd...                 */

    if (globals->smallest_vol == 1)
      add_speckle_to_volume(d1, 
			    globals->speckle,
			    globals->start, globals->count, globals->directions);
    else
      add_speckle_to_volume(d2, 
			    globals->speckle,
			    globals->start, globals->count, globals->directions);    
  } else if (globals->obj_function == vr_objective) {
    if (globals->smallest_vol == 1) {
      if (!build_segment_table(&segment_table, d1, globals->groups))
	print_error("%s\n",__FILE__, __LINE__,"Could not build segment table for source volume");
    }
    else {
      if (!build_segment_table(&segment_table, d2, globals->groups))
	print_error("%s\n",__FILE__, __LINE__,"Could not build segment table for target volume");

    }

    if (globals->flags.debug && globals->flags.verbose>1) {
      print ("groups = %d\n",segment_table->groups);
      for_less(i, segment_table->min, segment_table->max+1) {
	print ("%5d: table = %5d, function = %5d\n",i,segment_table->table[i],
	       (segment_table->segment)(i,segment_table) );
      }
    }

  }
	   /* ---------------- call requested optimization strategy ---------*/

  switch (globals->optimize_type) {
  case OPT_SIMPLEX:
    stat = optimize_simplex(d1, d2, m1, m2, globals);
    break;
  default:
    (void)fprintf(stderr, "Unknown type of optimization requested (%d)\n",globals->optimize_type);
    (void)fprintf(stderr, "Error in line %d, file %s\n",__LINE__, __FILE__);
    stat = FALSE;
  }
  
          /* ----------------finish up parameter/matrix manipulations ------*/

  if (globals->obj_function == vr_objective) {
    stat = free_segment_table(segment_table);
  }

  return(stat);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : measure_fit
                evaluate the objective function comparing d1 to d2 a single
		time.
@INPUT      : d1,d2:
                two volumes of data (already in memory).
	      m1,m2:
                two mask volumes for data (already in memory).
	      globals:
	        a global data structure containing info from the command line,
		including the input parameters, the input matrix,
		and a plethora of flags!
@OUTPUT     : 
@RETURNS    : TRUE if ok, FALSE if error.
@DESCRIPTION: 
@METHOD     :
                1- this routine begins by initializing the volume data structures
		to be used by the objective functions.

		2- the objective function is evaluated
		
		3- its value is returned
@GLOBALS    : 
@CALLS      : 
@CREATED    : Tue Jul 13 11:15:57 EST 1993 LC
@MODIFIED   : 
---------------------------------------------------------------------------- */
public float measure_fit(Volume d1,
			 Volume d2,
			 Volume m1,
			 Volume m2, 
			 Arg_Data *globals)
{
  BOOLEAN 
    stat;
  float 
    **p, y;
  int 
    i, 
    ndim;

/*  Transform
    *mat;
    double trans[3];  double cent[3];  double rots[3];  double scale[3];  double shear[3]; */

  
  stat = TRUE;

	     /*----------------- prepare data for objective function evaluation ------------ */

  
  if (globals->obj_function == zscore_objective) { /* replace volume d1 and d2 by zscore volume  */
    make_zscore_volume(d1,m1,(float)globals->threshold[0]);
    make_zscore_volume(d2,m2,(float)globals->threshold[1]);
  } 
  else  if (globals->obj_function == ssc_objective) {	/* add speckle to the data set */

    make_zscore_volume(d1,m1,(float)globals->threshold[0]); /* need to make data sets comparable */
    make_zscore_volume(d2,m2,(float)globals->threshold[1]); /* in mean and sd...                 */

    if (globals->smallest_vol == 1)
      add_speckle_to_volume(d1, 
			    globals->speckle,
			    globals->start, globals->count, globals->directions);
    else
      add_speckle_to_volume(d2, 
			    globals->speckle,
			    globals->start, globals->count, globals->directions);    
  } else if (globals->obj_function == vr_objective) {

    if (globals->smallest_vol == 1) {
      if (!build_segment_table(&segment_table, d1, globals->groups))
	print_error("%s",__FILE__, __LINE__,"Could not build segment table for source volume\n");
    }
    else {
      if (!build_segment_table(&segment_table, d2, globals->groups))
	print_error("%s",__FILE__, __LINE__,"Could not build segment table for target volume\n");
    }

  }
	   /* ---------------- build the parameters necessary for obj func evaluation ---------*/


  switch (globals->trans_info.transform_type) {
  case TRANS_PROCRUSTES: 
    ndim = 7;
    break;
  case TRANS_LSQ: 
    ndim = 12;
    break;
  case TRANS_LSQ6: 
    ndim = 6;
    break;
  case TRANS_LSQ7: 
    ndim = 7;
    break;
  case TRANS_LSQ9: 
    ndim = 9;
    break;
  case TRANS_LSQ12: 
    ndim = 12;
    break;
  default:
    (void)fprintf(stderr, "Unknown type of transformation requested (%d)\n",
		   globals->trans_info.transform_type);
    (void)fprintf(stderr, "Error in line %d, file %s\n",__LINE__, __FILE__);
    stat = FALSE;
  }

				/* set GLOBALS to communicate with the
				   function to be fitted!              */
  if (stat) {
    Gndim = ndim;
    if (globals->smallest_vol == 1) {
      Gdata1 = d1;
      Gdata2 = d2;
      Gmask1 = m1;
      Gmask2 = m2;
      Ginverse_mapping_flag = FALSE;
    }
    else {
      Gdata1 = d2;
      Gdata2 = d1;
      Gmask1 = m2;
      Gmask2 = m1;
      Ginverse_mapping_flag = TRUE;
    }
  }


  if (stat) {

    p = matrix(1,ndim+1,1,ndim); /* simplex */
    
    p[1][1]=globals->trans_info.translations[0];
    p[1][2]=globals->trans_info.translations[1];
    p[1][3]=globals->trans_info.translations[2];
    
    p[1][4]=globals->trans_info.rotations[0]; 
    p[1][5]=globals->trans_info.rotations[1]; 
    p[1][6]=globals->trans_info.rotations[2];
    
    if (ndim >= 7) p[1][7]=globals->trans_info.scales[0];
    if (ndim >7) {
      p[1][8]=globals->trans_info.scales[1];
      p[1][9]=globals->trans_info.scales[2];
    }
    
    if (ndim==12) {
      for_less( i, 0, 3 )		/* set shears */
	p[1][10+i] = globals->trans_info.shears[i];
    }


      
    y = fit_function(p[1]);	/* evaluate the objective  function */

    free_matrix(p,1,ndim+1,1,ndim); /* simplex */

  }
  
          /* ----------------finish up parameter/matrix manipulations ------*/

  if (globals->obj_function == vr_objective) {
    stat = free_segment_table(segment_table);
  }

  if (stat)
    return(y);
  else
    return(-FLT_MAX);
}

