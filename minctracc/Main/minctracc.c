/* ----------------------------- MNI Header -----------------------------------
   @NAME       : minctracc.c
   @INPUT      : base name of file corresponding to msp volume data, 
                 and output file namefor results
   
   minctracc [options] inputfile.iff modelfile.iff outputfile 
   
   @OUTPUT     : output file will contain the parameters for an affine transformation
                 that best maps the inputfile to the outputfile

   @RETURNS    : TRUE if ok, ERROR if error.

   @DESCRIPTION: This program will read in a volumetric dataset in mnc format
                 This file is assumed to contain rectangular voxels.  


   @METHOD     : This routine works using optimization, the goal of the program
         is to find the affine transformation T that best maps D to M.  

	 'best' is defined by the user-selected objective function to be minimized.

   @GLOBALS    : 
   @CALLS      : 
   @CREATED    : February 3, 1992 - louis collins
   @MODIFIED   : Thu May 20 11:21:22 EST 1993 lc
             complete re-write to use minc library

Wed May 26 13:05:44 EST 1993 lc

        complete re-write of fit_vol to:
	- use minc library (libminc.a) to read and write .mnc files
	- use libmni.a
	- use dave's volume structures.
	- read and write .xfm files, for the transformations
	- allow different interpolation schemes on voxel values
	- allow different objective functions: xcorr, ssc, var ratio
	- allow mask volumes on both source and model
	- calculate the bounding box of both source and model,
	  to map samples from smallest volume into larger one (for speed)

   ---------------------------------------------------------------------------- */

#include <def_mni.h>
#include "minctracc.h"
#include "globals.h"
#include <recipes.h>

/*************************************************************************/
main ( argc, argv )
     int argc;
     char *argv[];
{   /* main */
  
  global_data_struct 
    main_data;

  Status 
    status;
  Linear_Transformation 
    *lt;
  int 
    sizes[3],i,j;
  Real
    step[3];
  float 
    **the_matrix,  **inv_matrix;

  prog_name     = argv[0];	

  /* Call ParseArgv to interpret all command line args */
  if (ParseArgv(&argc, argv, argTable, 0) || (argc!=4)) {
    (void) fprintf(stderr, 
		   "\nUsage: %s [<options>] <sourcefile> <targetfile> <output transfile>\n", 
		   prog_name);
    (void) fprintf(stderr,"       %s [-help]\n\n", prog_name);
    exit(EXIT_FAILURE);
  }

  main_args.filenames.data  = argv[1];	/* set up necessary file names */
  main_args.filenames.model = argv[2];
  main_args.filenames.output_trans = argv[3];
  

				/* set up linear transformation to identity, if
				   not set up in the argument list */

  if (main_args.trans_info.transformation.transform == NULL) { /* then use identity */
    main_args.trans_info.transformation = identity_transformation;
    main_args.trans_info.transformation.trans_data = malloc(sizeof(Linear_Transformation));
    (void) memcpy(main_args.trans_info.transformation.trans_data,
		  identity_transformation.trans_data,
		  sizeof(Linear_Transformation));
  }

				/* set up main data struct, for file writing below */
  main_data.n_tags = 0;
  main_data.n_vols = 2;
  main_data.dim    = 3;
  main_data.tags1  = NULL;
  main_data.tags2  = NULL;
  main_data.labels = NULL;
  main_data.tps    = NULL;
  strcpy(main_data.filename1,main_args.filenames.data);
  strcpy(main_data.filename2,main_args.filenames.model);



  DEBUG_PRINT1 ( "===== Debugging information from %s =====\n", prog_name);
  DEBUG_PRINT1 ( "Data filename       = %s\n", main_args.filenames.data);
  DEBUG_PRINT1 ( "Model filename      = %s\n", main_args.filenames.model);
  DEBUG_PRINT1 ( "Data mask filename  = %s\n", main_args.filenames.mask_data);
  DEBUG_PRINT1 ( "Model mask filename = %s\n", main_args.filenames.mask_model);
  DEBUG_PRINT1 ( "Output filename     = %s\n\n", main_args.filenames.output_trans);
  DEBUG_PRINT  ( "Objective function  = ");
  if (main_args.obj_function == xcorr_objective) {
    DEBUG_PRINT1("cross correlation (threshold = %f)\n",main_args.threshold);
  }
  else
    if (main_args.obj_function == zscore_objective) {
      DEBUG_PRINT1( "zscore (threshold = %f)\n", main_args.threshold );
    }
    else
      if (main_args.obj_function == vr_objective) {
	DEBUG_PRINT2( "ratio of variance (threshold = %f, groups = %d)\n", 
		     main_args.threshold, main_args.groups);
      }
      else
	if (main_args.obj_function == ssc_objective) {
	  DEBUG_PRINT2( "stochastic sign change (threshold = %f, speckle %f %%)\n",
		       main_args.threshold, main_args.speckle);
	}
	else {
	  DEBUG_PRINT("unknown!\n");
	  exit(EXIT_FAILURE);
	}


  DEBUG_PRINT1 ( "Transform linear    = %s\n", (main_args.trans_info.transformation.linear ? "TRUE" : "FALSE") );
  DEBUG_PRINT1 ( "Transform inverted? = %s\n", (main_args.trans_info.invert_mapping_flag ? "TRUE" : "FALSE") );
  DEBUG_PRINT1 ( "Transform type      = %d\n", main_args.trans_info.transform_type );

  if (main_args.trans_info.transformation.linear) {
    lt = (Linear_Transformation *) main_args.trans_info.transformation.trans_data;

    DEBUG_PRINT ( "Transform matrix    = ");
    for_less(i,0,4) DEBUG_PRINT1 ("%9.4f ",lt->mat[0][i]);
    DEBUG_PRINT ( "\n" );
    DEBUG_PRINT ( "                      ");
    for_less(i,0,4) DEBUG_PRINT1 ("%9.4f ",lt->mat[1][i]);
    DEBUG_PRINT ( "\n" );
    DEBUG_PRINT ( "                      ");
    for_less(i,0,4) DEBUG_PRINT1 ("%9.4f ",lt->mat[2][i]);
    DEBUG_PRINT ( "\n" );
  }

  DEBUG_PRINT3 ( "Transform center   = %8.3f %8.3f %8.3f\n", 
		main_args.trans_info.center[0],
		main_args.trans_info.center[1],
		main_args.trans_info.center[2] );
  DEBUG_PRINT3 ( "Transform trans    = %8.3f %8.3f %8.3f\n", 
		main_args.trans_info.translations[0],
		main_args.trans_info.translations[1],
		main_args.trans_info.translations[2] );
  DEBUG_PRINT3 ( "Transform rotation = %8.3f %8.3f %8.3f\n", 
		main_args.trans_info.rotations[0],
		main_args.trans_info.rotations[1],
		main_args.trans_info.rotations[2] );
  DEBUG_PRINT3 ( "Transform scale    = %8.3f %8.3f %8.3f\n\n", 
		main_args.trans_info.scales[0],
		main_args.trans_info.scales[1],
		main_args.trans_info.scales[2] );


  ALLOC( data, 1 );		/* read in source data and target model */
  ALLOC( model, 1 );
  status = input_volume( main_args.filenames.data, &data );
  status = input_volume( main_args.filenames.model, &model );

  get_volume_sizes(data, sizes);
  get_volume_separations(data, step);
  DEBUG_PRINT3 ( "Source volume %3d cols by %3d rows by %d slices\n",
		 sizes[X], sizes[Y], sizes[Z]);
  DEBUG_PRINT3 ( "Source voxel = %8.3f %8.3f %8.3f\n", 
		 step[X], step[Y], step[Z]);
  get_volume_sizes(model, sizes);
  get_volume_separations(model, step);
  DEBUG_PRINT3 ( "Target volume %3d cols by %3d rows by %d slices\n",
		 sizes[X], sizes[Y], sizes[Z]);
  DEBUG_PRINT3 ( "Target voxel = %8.3f %8.3f %8.3f\n\n", 
		 step[X], step[Y], step[Z]);
  


  /* ===========================  translate initial transformation matrix into 
                                  transformation parameters */

  init_params( data, model, mask_data, mask_model, &main_args );

  DEBUG_PRINT  ("AFTER init_params()\n");
  if (main_args.trans_info.transformation.linear) {
    lt = (Linear_Transformation *) main_args.trans_info.transformation.trans_data;

    DEBUG_PRINT ( "Transform matrix    = ");
    for_less(i,0,4) DEBUG_PRINT1 ("%9.4f ",lt->mat[0][i]);
    DEBUG_PRINT ( "\n" );
    DEBUG_PRINT ( "                      ");
    for_less(i,0,4) DEBUG_PRINT1 ("%9.4f ",lt->mat[1][i]);
    DEBUG_PRINT ( "\n" );
    DEBUG_PRINT ( "                      ");
    for_less(i,0,4) DEBUG_PRINT1 ("%9.4f ",lt->mat[2][i]);
    DEBUG_PRINT ( "\n" );
  }
  DEBUG_PRINT ( "\n" );

  DEBUG_PRINT3 ( "Transform center   = %8.3f %8.3f %8.3f\n", 
		main_args.trans_info.center[0],
		main_args.trans_info.center[1],
		main_args.trans_info.center[2] );
  DEBUG_PRINT3 ( "Transform trans    = %8.3f %8.3f %8.3f\n", 
		main_args.trans_info.translations[0],
		main_args.trans_info.translations[1],
		main_args.trans_info.translations[2] );
  DEBUG_PRINT3 ( "Transform rotation = %8.3f %8.3f %8.3f\n", 
		main_args.trans_info.rotations[0],
		main_args.trans_info.rotations[1],
		main_args.trans_info.rotations[2] );
  DEBUG_PRINT3 ( "Transform scale    = %8.3f %8.3f %8.3f\n\n", 
		main_args.trans_info.scales[0],
		main_args.trans_info.scales[1],
		main_args.trans_info.scales[2] );



				/* do not do any optimization if the transformation
				   requested is the Principal Axes Transformation 
				   then:
		                   =======   do linear fitting =============== */

  if (main_args.trans_info.transform_type != TRANS_PAT) {
    
				/* initialize the sampling lattice and figure out
				   which of the two volumes is smaller.           */

    init_lattice( data, model, mask_data, mask_model, &main_args );


    if (main_args.smallest_vol == 1) {
      DEBUG_PRINT("Source volume is smallest\n");
    }
    else {
      DEBUG_PRINT("Target volume is smallest\n");
    }
    DEBUG_PRINT3 ( "Lattice step size  = %8.3f %8.3f %8.3f\n",
		  main_args.step[0],main_args.step[1],main_args.step[2]);
    DEBUG_PRINT3 ( "Lattice start      = %8.3f %8.3f %8.3f\n",
		  main_args.start[0],main_args.start[1],main_args.start[2]);
    DEBUG_PRINT3 ( "Lattice count      = %8d %8d %8d\n\n",
		  main_args.count[0],main_args.count[1],main_args.count[2]);

    if (!optimize_linear_transformation( data, model, mask_data, mask_model, &main_args )) {
      (void) fprintf(stderr, 
		     "\n%s: Error in optimization\n", 
		     prog_name);
      exit(EXIT_FAILURE);
   }
  }



  /* ===========================   write out transformation =============== */

  if (main_args.trans_info.transformation.linear) {
    
    lt = (Linear_Transformation *) main_args.trans_info.transformation.trans_data;
    
    for_less ( i, 0, 3 )       /* copy to main_data transform matrix */
      for_less ( j, 0, 4)
	main_data.xform[i][j] = lt->mat[i][j];
    
  } 
  else {
    (void) fprintf(stderr, 
		   "\n%s: Error - non-linear transformations not yet supported\n", 
		   prog_name);
    exit(EXIT_FAILURE);
  }

  if (main_args.trans_info.invert_mapping_flag && main_args.trans_info.transformation.linear) {
    
    the_matrix= matrix(1,4,1,4);
    inv_matrix= matrix(1,4,1,4);
    
    for (i=1; i<=4; ++i) {      /* copy to main transform matrix */
      for (j=1; j<=3; ++j)
	the_matrix[i][j] = main_data.xform[j-1][i-1];
      if (i!=4)
	the_matrix[i][4] = 0.0;
      else
	the_matrix[i][4] = 1.0;
    }
    
    invertmatrix(4,the_matrix,inv_matrix);
    
    for (i=1; i<=4; ++i)        /* copy to main transform matrix */
      for (j=1; j<=3; ++j)
	main_data.xform[j-1][i-1] = inv_matrix[i][j];
    
    free_matrix(inv_matrix,1,4,1,4);
    free_matrix(the_matrix,1,4,1,4);
  }
  
  if (!save_transform(&main_data, main_args.filenames.output_trans)) {
    (void) fprintf(stderr,"Error saving transformation to %s.\n",main_args.filenames.output_trans);
    return ERROR_STATUS;
  }
  
  return( status );
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_transformation
@INPUT      : dst - Pointer to client data from argument table
              key - argument key
              nextArg - argument following key
@OUTPUT     : (nothing) 
@RETURNS    : TRUE so that ParseArgv will discard nextArg
@DESCRIPTION: Routine called by ParseArgv to read in a transformation file
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : February 15, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public int get_transformation(char *dst, char *key, char *nextArg)
{     /* ARGSUSED */
   Transformation *transformation;
   Linear_Transformation *matrx;
   Thin_plate_spline       *tps;
   FILE *fp;
   
   int dim;
   int num_pts;
   double **warps;
   double **coords;

   /* Get pointer to transformation structure */
   transformation = (Transformation *) dst;

   /* Open the file */
   if (strcmp(nextArg, "-") == 0) {
      fp = stdin;
   }
   else {
      fp = fopen(nextArg, "r");
      if (fp==NULL) {
         (void) fprintf(stderr, "Error opening transformation file.\n");
         exit(EXIT_FAILURE);
      }
   }

   /* figure out what type of transformation file, linear or thin plate spline */
   
   if ( is_tps_file(fp) ) {
     transformation->linear = FALSE;
     transformation->transform = do_non_linear_transformation;
     if (transformation->trans_data != NULL) FREE(transformation->trans_data);
     ALLOC( tps, 1 );
     transformation->trans_data = tps;

     if (!input_tps_transform(fp,
			      &num_pts,
			      &dim,
			      &warps,
			      &coords)) {
       (void) fprintf(stderr, "Cannot input Thin Plate Spline transformation.\n");
       return(FALSE);
     } 
     else {
       tps->num_points = num_pts;
       tps->dim = dim;
       tps->warps = warps;
       tps->coords = coords;
     }
   }
   else {
				/* Must be a linear transformation */
     transformation->linear = TRUE;
     transformation->transform = do_linear_transformation;
     if (transformation->trans_data != NULL) FREE(transformation->trans_data);
     ALLOC( matrx, 1 );
     transformation->trans_data = matrx;
     
     /* Read the file */
     if (!input_transform(fp, matrx->mat)) {
       (void) fprintf(stderr, "Error reading transformation file.\n");
       exit(EXIT_FAILURE);
     }
     
     /* Invert the transformation */

/*
     invert_transformation(transformation, transformation);  only needed for resampling!
*/
   }
   

   /* set a GLOBAL flag, to show that a transformation has been read in */

   main_args.trans_info.use_default = FALSE;	

   return TRUE;
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_mask_file
@INPUT      : dst - Pointer to client data from argument table
              key - argument key
              nextArg - argument following key
@OUTPUT     : (nothing) 
@RETURNS    : TRUE so that ParseArgv will discard nextArg
@DESCRIPTION: Routine called by ParseArgv to read in a binary mask file
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Wed May 26 13:05:44 EST 1993 Collins Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */
public int get_mask_file(char *dst, char *key, char *nextArg)
{     /* ARGSUSED */

  Status status;

  if (strncmp ( "-model_mask", key, 2) == 0) {
    ALLOC( mask_model, 1 );
    status = input_volume( nextArg, &mask_model );
    dst = nextArg;
  }
  else {
    ALLOC( mask_data, 1);
    status = input_volume( nextArg, &mask_data );
    dst = nextArg;
  }

  if (status != OK)
  {
    (void) fprintf(stderr, "Cannot input mask file %s.",nextArg);
    return(FALSE);
  } 

  return TRUE;
  
}


public void print_error(char *s, char * d1, int d2, int d3, int d4, int d5, int d6, int d7)
{
  (void) fprintf(stderr, "Error in %s in file %s, line %d\n",prog_name,d1,d2);
  (void) fprintf(stderr, "   %s\n", s, d3,d4,d5,d6,d7);
  exit(EXIT_FAILURE);

}
