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
   @MODIFIED   : $Log: minctracc.c,v $
   @MODIFIED   : Revision 1.7  1993-11-15 13:12:10  louis
   @MODIFIED   : working version, deform deform installation
   @MODIFIED   :

Thu May 20 11:21:22 EST 1993 lc
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

@MODIFIED   : $Log: minctracc.c,v $
@MODIFIED   : Revision 1.7  1993-11-15 13:12:10  louis
@MODIFIED   : working version, deform deform installation
@MODIFIED   :
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Main/minctracc.c,v 1.7 1993-11-15 13:12:10 louis Exp $";
#endif




#include <mni.h>
#include <recipes.h>
#include <limits.h>
#include <print_error.h>

#include "minctracc.h"
#include "globals.h"

static char *default_dim_names[N_DIMENSIONS] =
   { MIzspace, MIyspace, MIxspace };


/*************************************************************************/
main ( argc, argv )
     int argc;
     char *argv[];
{   /* main */
  

  Status 
    status;
  General_transform
    tmp_invert;
  
  Transform 
    *lt, ident_trans;
  int 
    sizes[3],msizes[3],i;
  Real
    min_value, max_value, step[3];
  float 
    obj_func_val;

  FILE *ofd;

  char 
    *comments;

  Real thickness[3];





  prog_name     = argv[0];	


  comments = time_stamp(argc, argv); /* build comment history line for below */


  /* Call ParseArgv to interpret all command line args */
  if (ParseArgv(&argc, argv, argTable, 0) || 
      ((argc!=4) && (main_args.filenames.measure_file==NULL)) || 
      ((argc!=3) && (main_args.filenames.measure_file!=NULL))
      ) {

    print ("Parameters left:\n");
    for_less(i,0,argc)
      print ("%s ",argv[i]);
    print ("\n");
    
    (void) fprintf(stderr, 
		   "\nUsage: %s [<options>] <sourcefile> <targetfile> <output transfile>\n", 
		   prog_name);
    (void) fprintf(stderr,"       %s [-help]\n\n", prog_name);
    exit(EXIT_FAILURE);
  }

				/* reset the order of step sizes!,
				   since they come in x,y,z and the data is forced
				   to be read in as z,y,x                           */ 
  thickness[COL_IND] = main_args.step[0]; 
  thickness[ROW_IND] = main_args.step[1];
  thickness[SLICE_IND] = main_args.step[2];
  for_less(i,0,3) 
    main_args.step[i]=thickness[i];




  main_args.filenames.data  = argv[1];	/* set up necessary file names */
  main_args.filenames.model = argv[2];
  if (main_args.filenames.measure_file==NULL) 
    main_args.filenames.output_trans = argv[3];

  if ( (main_args.filenames.measure_file==NULL) && !clobber_flag && file_exists(argv[3])) {
    print ("File %s exists.\n",argv[3]);
    print ("Use -clobber to overwrite.\n");
    return ERROR;
  }

  

				/* set up linear transformation to identity, if
				   not set up in the argument list */

  if (main_args.trans_info.transformation == (General_transform *)NULL) { /* then use identity */

				/* build identity General_transform */
    make_identity_transform(&ident_trans);
				/* attach it to global w-to-w transformation */
    ALLOC(main_args.trans_info.transformation,1);
    create_linear_transform(main_args.trans_info.transformation,&ident_trans);

  }

  ALLOC(main_args.trans_info.orig_transformation,1);
  copy_general_transform(main_args.trans_info.transformation,
			 main_args.trans_info.orig_transformation);

  DEBUG_PRINT1 ( "===== Debugging information from %s =====\n", prog_name);
  DEBUG_PRINT1 ( "Data filename       = %s\n", main_args.filenames.data);
  DEBUG_PRINT1 ( "Model filename      = %s\n", main_args.filenames.model);
  DEBUG_PRINT1 ( "Data mask filename  = %s\n", main_args.filenames.mask_data);
  DEBUG_PRINT1 ( "Model mask filename = %s\n", main_args.filenames.mask_model);
  DEBUG_PRINT1 ( "Output filename     = %s\n", main_args.filenames.output_trans);
  if (main_args.filenames.matlab_file != NULL)
    DEBUG_PRINT1 ( "Matlab filename     = %s\n\n",main_args.filenames.matlab_file );
  if (main_args.filenames.measure_file != NULL)
    DEBUG_PRINT1 ( "Measure filename    = %s\n\n",main_args.filenames.measure_file );
  DEBUG_PRINT  ( "Objective function  = ");
  if (main_args.obj_function == xcorr_objective) {
    DEBUG_PRINT2("cross correlation (threshold = %f %f)\n",
		 main_args.threshold[0],main_args.threshold[1]);
  }
  else
    if (main_args.obj_function == zscore_objective) {
      DEBUG_PRINT2( "zscore (threshold = %f %f)\n", main_args.threshold[0],main_args.threshold[1] );
    }
    else
      if (main_args.obj_function == vr_objective) {
	DEBUG_PRINT3( "ratio of variance (threshold = %f %f, groups = %d)\n", 
		     main_args.threshold[0],main_args.threshold[1], main_args.groups);
      }
      else
	if (main_args.obj_function == ssc_objective) {
	  DEBUG_PRINT3( "stochastic sign change (threshold = %f %f, speckle %f %%)\n",
		       main_args.threshold[0],main_args.threshold[1], main_args.speckle);
	}
	else {
	  DEBUG_PRINT("unknown!\n");
	  (void) fprintf(stderr,"Unknown objective function requested.\n");
	  exit(EXIT_FAILURE);
	}


  DEBUG_PRINT1 ( "Transform linear    = %s\n", 
		(get_transform_type(main_args.trans_info.transformation)==LINEAR ? 
		 "TRUE" : "FALSE") );
  DEBUG_PRINT1 ( "Transform inverted? = %s\n", (main_args.trans_info.invert_mapping_flag ? 
						"TRUE" : "FALSE") );
  DEBUG_PRINT1 ( "Transform type      = %d\n", main_args.trans_info.transform_type );

  if (get_transform_type(main_args.trans_info.transformation)==LINEAR) {
    lt = get_linear_transform_ptr(main_args.trans_info.transformation);

    DEBUG_PRINT ( "Transform matrix    = ");
    for_less(i,0,4) DEBUG_PRINT1 ("%9.4f ",Transform_elem(*lt,0,i));
    DEBUG_PRINT ( "\n" );
    DEBUG_PRINT ( "                      ");
    for_less(i,0,4) DEBUG_PRINT1 ("%9.4f ",Transform_elem(*lt,1,i));
    DEBUG_PRINT ( "\n" );
    DEBUG_PRINT ( "                      ");
    for_less(i,0,4) DEBUG_PRINT1 ("%9.4f ",Transform_elem(*lt,2,i));
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
  status = input_volume( main_args.filenames.data, 3, default_dim_names, 
			NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			TRUE, &data, (minc_input_options *)NULL );
  if (status != OK)
    print_error("Cannot input volume '%s'",__FILE__, __LINE__,main_args.filenames.data);
  get_volume_separations(data, step);
  get_volume_sizes(data, sizes);
  get_volume_voxel_range(data, &min_value, &max_value);

  DEBUG_PRINT3 ( "Source volume %3d cols by %3d rows by %d slices\n",
		 sizes[X], sizes[Y], sizes[Z]);
  DEBUG_PRINT3 ( "Source voxel = %8.3f %8.3f %8.3f\n", 
		 step[X], step[Y], step[Z]);
  DEBUG_PRINT2 ( "min/max value= %8.3f %8.3f\n", min_value, max_value);
  get_volume_voxel_range(data, &min_value, &max_value);
  DEBUG_PRINT2 ( "min/max voxel= %8.3f %8.3f\n\n", min_value, max_value);



  status = input_volume( main_args.filenames.model, 3, default_dim_names, 
			NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			TRUE, &model, (minc_input_options *)NULL );
  if (status != OK)
    print_error("Cannot input volume '%s'",
		__FILE__, __LINE__,main_args.filenames.model);
  get_volume_separations(model, step);
  get_volume_sizes(model, sizes);
  get_volume_voxel_range(model, &min_value, &max_value);

  DEBUG_PRINT3 ( "Target volume %3d cols by %3d rows by %d slices\n",
		 sizes[X], sizes[Y], sizes[Z]);
  DEBUG_PRINT3 ( "Target voxel = %8.3f %8.3f %8.3f\n", 
		 step[X], step[Y], step[Z]);
  DEBUG_PRINT2 ( "min/max value= %8.3f %8.3f\n", min_value, max_value);
  get_volume_voxel_range(model, &min_value, &max_value);
  DEBUG_PRINT2 ( "min/max voxel= %8.3f %8.3f\n\n", min_value, max_value);

  
  if (data->n_dimensions!=3) {
    print_error ("File %s has %d dimensions.  Only 3 dims supported.", 
		 __FILE__, __LINE__, main_args.filenames.data, data->n_dimensions);
  }
  if (model->n_dimensions!=3) {
    print_error ("File %s has %d dimensions.  Only 3 dims supported.", 
		 __FILE__, __LINE__, main_args.filenames.model, model->n_dimensions);
  }


  /* ===========================  translate initial transformation matrix into 
                                  transformation parameters */

  init_params( data, model, mask_data, mask_model, &main_args );

  if (main_args.filenames.matlab_file != NULL) {
    init_lattice( data, model, mask_data, mask_model, &main_args );
    make_matlab_data_file( data, model, mask_data, mask_model, comments, &main_args );
    exit( OK );
  }

  DEBUG_PRINT  ("AFTER init_params()\n");
  if (get_transform_type(main_args.trans_info.transformation)==LINEAR) {
    lt = get_linear_transform_ptr(main_args.trans_info.transformation);

    DEBUG_PRINT ( "Transform matrix    = ");
    for_less(i,0,4) DEBUG_PRINT1 ("%9.4f ",Transform_elem(*lt,0,i));
    DEBUG_PRINT ( "\n" );
    DEBUG_PRINT ( "                      ");
    for_less(i,0,4) DEBUG_PRINT1 ("%9.4f ",Transform_elem(*lt,1,i));
    DEBUG_PRINT ( "\n" );
    DEBUG_PRINT ( "                      ");
    for_less(i,0,4) DEBUG_PRINT1 ("%9.4f ",Transform_elem(*lt,2,i));
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
  


  if (main_args.filenames.measure_file != NULL) {


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


    status = open_file(  main_args.filenames.measure_file, WRITE_FILE, BINARY_FORMAT,  &ofd );
    if ( status != OK ) 
      print_error ("filename `%s' cannot be opened.", 
		   __FILE__, __LINE__, main_args.filenames.measure_file);

				/* do zscore and reload */

    main_args.obj_function = zscore_objective;
    obj_func_val = measure_fit( data, model, mask_data, mask_model, &main_args );
    (void)fprintf (ofd, "%f - zscore\n",obj_func_val);
    (void)fflush(ofd);
    DEBUG_PRINT1 ( "%f - zscore\n",obj_func_val);
    
    delete_volume(data);  
    delete_volume(model);  


    status = input_volume( main_args.filenames.data, 3, default_dim_names, 
			  NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			  TRUE, &data, (minc_input_options *)NULL );
    status = input_volume( main_args.filenames.model, 3, default_dim_names, 
			  NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			  TRUE, &model, (minc_input_options *)NULL );

				/* do xcorr   */

    main_args.obj_function = xcorr_objective;
    obj_func_val = measure_fit( data, model, mask_data, mask_model, &main_args );
    (void)fprintf (ofd, "%f - xcorr\n",obj_func_val);
    (void)fflush(ofd);
    DEBUG_PRINT1 ( "%f - xcorr\n",obj_func_val);

				/* do var_ratio */

    main_args.obj_function = vr_objective;
    obj_func_val = measure_fit( data, model, mask_data, mask_model, &main_args );
    (void)fprintf (ofd, "%f - var_ratio\n",obj_func_val);
    (void)fflush(ofd);
    DEBUG_PRINT1 ( "%f - var_ratio\n",obj_func_val);

    delete_volume(data);  
    delete_volume(model);  

				/* do ssc / zero-crossings */

    status = input_volume( main_args.filenames.data, 3, default_dim_names, 
			  NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			  TRUE, &data, (minc_input_options *)NULL );
    status = input_volume( main_args.filenames.model, 3, default_dim_names, 
			  NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			  TRUE, &model, (minc_input_options *)NULL ); 

    main_args.obj_function = ssc_objective;
    obj_func_val = measure_fit( data, model, mask_data, mask_model, &main_args );
    (void)fprintf (ofd, "%f - ssc\n",obj_func_val);
    (void)fflush(ofd);
    DEBUG_PRINT1 ( "%f - ssc\n",obj_func_val);


    status = close_file(ofd);
    if ( status != OK ) 
      print_error ("filename `%s' cannot be closed.", 
		   __FILE__, __LINE__, main_args.filenames.measure_file);


    exit( status ) ;

  }

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

				/* if I have internally inverted the transform,
				   than flip it back forward before the save.   */

  if (main_args.trans_info.invert_mapping_flag) {

    DEBUG_PRINT ("Re-inverting transformation\n");
    create_inverse_general_transform(main_args.trans_info.transformation,
				     &tmp_invert);

    copy_general_transform(&tmp_invert,main_args.trans_info.transformation);

  }
				/* note: - I use the comment string built above! */
				/* open file */

  status = open_file(  main_args.filenames.output_trans, WRITE_FILE, BINARY_FORMAT,  &ofd );
  if ( status != OK ) 
    print_error ("filename `%s' cannot be opened.", 
		 __FILE__, __LINE__, main_args.filenames.output_trans);
  
  				/* save transformation */

  status = output_transform(ofd, comments, main_args.trans_info.transformation );

				/* close file */
  if (status == OK)
    status = close_file(ofd);
  else
    print_error ("Problems writing  `%s'.",
		 __FILE__, __LINE__, main_args.filenames.output_trans);

  
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
   Program_Transformation *transform_info;
   General_transform *transformation;
   General_transform input_transformation;
   FILE *fp;
   int ch, index;

   /* Check for following argument */
   if (nextArg == NULL) {
      (void) fprintf(stderr, 
                     "\"%s\" option requires an additional argument\n",
                     key);
      return FALSE;
   }

   /* Get pointer to transform info structure */
   transform_info = (Program_Transformation *) dst;

   /* Save file name */
   transform_info->file_name = nextArg;

   if (transform_info->transformation == (General_transform *)NULL) 
     ALLOC(transform_info->transformation, 1);

   transformation =  transform_info->transformation;


   /* Open the file */
   if (strcmp(nextArg, "-") == 0) {
      /* Create a temporary for standard input */
      fp=tmpfile();
      if (fp==NULL) {
         (void) fprintf(stderr, "Error opening temporary file.\n");
         exit(EXIT_FAILURE);
      }
      while ((ch=getc(stdin))!=EOF) (void) putc(ch, fp);
      rewind(fp);
   }
   else {
      fp = fopen(nextArg, "r");
      if (fp==NULL) {
         (void) fprintf(stderr, "Error opening transformation file %s.\n",
                        nextArg);
         exit(EXIT_FAILURE);
      }
   }

   /* Read in the file for later use */
   if (transform_info->file_contents == NULL) {
      ALLOC(transform_info->file_contents,TRANSFORM_BUFFER_INCREMENT);
      transform_info->buffer_length = TRANSFORM_BUFFER_INCREMENT;
   }
   for (index = 0; (ch=getc(fp)) != EOF; index++) {
      if (index >= transform_info->buffer_length-1) {
         transform_info->buffer_length += TRANSFORM_BUFFER_INCREMENT;
	 REALLOC(transform_info->file_contents, 
		 transform_info->buffer_length);
      }
      transform_info->file_contents[index] = ch;
   }
   transform_info->file_contents[index] = '\0';
   rewind(fp);

   /* Read the file */
   if (input_transform(fp, &input_transformation)!=OK) {
      (void) fprintf(stderr, "Error reading transformation file.\n");
      exit(EXIT_FAILURE);
   }
   (void) fclose(fp);

   /* set a GLOBAL flag, to show that a transformation has been read in */

   main_args.trans_info.use_default = FALSE;


   copy_general_transform(&input_transformation,transformation);

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
@MODIFIED   : Wed Jun 16 13:21:18 EST 1993 LC
    added: set of main_args.filenames.mask_model and .mask_data

---------------------------------------------------------------------------- */
public int get_mask_file(char *dst, char *key, char *nextArg)
{     /* ARGSUSED */

  Status status;

  if (strncmp ( "-model_mask", key, 2) == 0) {
    ALLOC( mask_model, 1 );
    status = input_volume( nextArg, 3, default_dim_names, 
			  NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			  TRUE, &mask_model, (minc_input_options *)NULL );
    dst = nextArg;
    main_args.filenames.mask_model = nextArg;
  }
  else {
    ALLOC( mask_data, 1);
    status = input_volume( nextArg, 3, default_dim_names, 
			  NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			  TRUE, &mask_data, (minc_input_options *)NULL );
    dst = nextArg;
    main_args.filenames.mask_data = nextArg;
  }

  if (status != OK)
  {
    (void) fprintf(stderr, "Cannot input mask file %s.",nextArg);
    return(FALSE);
  } 

  return TRUE;
  
}




