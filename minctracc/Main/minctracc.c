/* ----------------------------- MNI Header -----------------------------------
   @NAME       : minctracc.c
   @COPYRIGHT  :
              Copyright 1993 Louis Collins, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.

   @CREATED    : February 3, 1992 - louis louis
   @MODIFIED   : $Log: minctracc.c,v $
   @MODIFIED   : Revision 1.14  1995-09-28 11:54:52  louis
   @MODIFIED   : working version, just prior to release 0.9 of mni_autoreg
   @MODIFIED   :
 * Revision 1.13  1995/02/22  08:56:06  louis
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.12  94/05/28  16:18:54  louis
 * working version before modification of non-linear optimiation
 * 
 * Revision 1.11  94/04/26  12:54:30  louis
 * updated with new versions of make_rots, extract2_parameters_from_matrix 
 * that include proper interpretation of skew.
 * 
 * Revision 1.10  94/04/06  11:48:43  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.9  94/02/21  16:35:51  louis
 * version before feb 22 changes
 * 
 * Revision 1.8  93/11/15  16:27:06  louis
 * working version, with new library, with RCS revision stuff,
 * before deformations included
 * 
 * Revision 1.7  93/11/15  13:12:10  louis
 * working version, before deform installation
 * 

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

---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Main/minctracc.c,v 1.14 1995-09-28 11:54:52 louis Exp $";
#endif

#include <limits.h>
#include <volume_io.h>
#include <print_error.h>
#include <minctracc.h>
#include <globals.h>

       Real initial_corr;
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
    sizes[3],i;
  Real
    min_value, max_value, step[3];
  char 
    *comments;
  FILE
    *ofd;
  Real
    obj_func_val;

  prog_name     = argv[0];	

  comments = time_stamp(argc, argv); /* build comment history line for below */


  /* Call ParseArgv to interpret all command line args */
  if (ParseArgv(&argc, argv, argTable, 0) || 
      ((argc!=4) && (strlen(main_args.filenames.measure_file) == 0)) || 
      ((argc!=3) && (strlen(main_args.filenames.measure_file) != 0))
      ) {
    
    print ("Parameters left:\n");
    for_less(i,0,argc)
      print ("%s ",argv[i]);
    print ("\n");
    
    (void)fprintf(stderr, 
		  "\nUsage: %s [<options>] <sourcefile> <targetfile> <output transfile>\n", 
		  prog_name);
    (void)fprintf(stderr,"       %s [-help]\n\n", prog_name);
    exit(EXIT_FAILURE);
  }

  main_args.filenames.data  = argv[1];	/* set up necessary file names */
  main_args.filenames.model = argv[2];
  if (strlen(main_args.filenames.measure_file)==0) 
    main_args.filenames.output_trans = argv[3];

  if ( (strlen(main_args.filenames.measure_file)==0) && !clobber_flag && file_exists(argv[3])) {
    (void)fprintf (stderr,"File %s exists.\n",argv[3]);
    (void)fprintf (stderr,"Use -clobber to overwrite.\n");
    exit(EXIT_FAILURE);
  }

				/* set up linear transformation to identity, if
				   not set up in the argument list */

  if (main_args.trans_info.transformation == (General_transform *)NULL) { 
				/* build identity General_transform */
    make_identity_transform(&ident_trans);
				/* attach it to global w-to-w transformation */
    ALLOC(main_args.trans_info.transformation,1);
    create_linear_transform(main_args.trans_info.transformation,&ident_trans);    
  } 
				/* don't use default PAT if user req -identity */
  if (main_args.trans_info.use_identity)
    main_args.trans_info.use_default = FALSE;

				/* make a copy of the original transformation */
  ALLOC(main_args.trans_info.orig_transformation,1);
  copy_general_transform(main_args.trans_info.transformation,
			 main_args.trans_info.orig_transformation);

  if (main_args.flags.debug) {
    /*
       this is simply to print out debugging info at the beginning of a run
       */

    print ( "===== Debugging information from %s =====\n", prog_name);
    print ( "Data filename       = %s\n", main_args.filenames.data);
    print ( "Model filename      = %s\n", main_args.filenames.model);
    print ( "Data mask filename  = %s\n", main_args.filenames.mask_data);
    print ( "Model mask filename = %s\n", main_args.filenames.mask_model);
    print ( "Input xform name    = %s\n", main_args.trans_info.file_name);
    print ( "Output filename     = %s\n", main_args.filenames.output_trans);
    if (strlen(main_args.filenames.matlab_file) != 0)
      print ( "Matlab filename     = %s\n\n",main_args.filenames.matlab_file );
    if (strlen(main_args.filenames.measure_file) != 0)
      print ( "Measure filename    = %s\n\n",main_args.filenames.measure_file );
    print ( "Step size           = %f %f %f\n",
	          main_args.step[0],
		  main_args.step[1],
		  main_args.step[2]);
    print  ( "Objective function  = ");
    if (main_args.obj_function == xcorr_objective) {
      print("cross correlation (threshold = %f %f)\n",
		   main_args.threshold[0],main_args.threshold[1]);
    }
    else
      if (main_args.obj_function == zscore_objective) {
	print( "zscore (threshold = %f %f)\n", main_args.threshold[0],main_args.threshold[1] );
      }
      else
	if (main_args.obj_function == vr_objective) {
	  print( "ratio of variance (threshold = %f %f, groups = %d)\n", 
		       main_args.threshold[0],main_args.threshold[1], main_args.groups);
	}
	else
	  if (main_args.obj_function == ssc_objective) {
	    print( "stochastic sign change (threshold = %f %f, speckle %f %%)\n",
			 main_args.threshold[0],main_args.threshold[1], main_args.speckle);
	  }
	  else
	    if (main_args.obj_function == mutual_information_objective) {
	      print( "mutual information (groups = %d) \n",
			   main_args.groups);
	    }
	    else {
	      print("unknown!\n");
	      (void)fprintf(stderr,"Unknown objective function requested.\n");
	      exit(EXIT_FAILURE);
	    }


    print ( "Transform linear    = %s\n", 
	   (get_transform_type(main_args.trans_info.transformation)==LINEAR ? 
	    "TRUE" : "FALSE") );
    print ( "Transform inverted? = %s\n", (main_args.trans_info.invert_mapping_flag ? 
					   "TRUE" : "FALSE") );
    print ( "Transform type      = %d\n", main_args.trans_info.transform_type );

    if (get_transform_type(main_args.trans_info.transformation)==LINEAR) {
      lt = get_linear_transform_ptr(main_args.trans_info.transformation);
      print ( "Transform matrix    = ");
      for_less(i,0,4) print ("%9.4f ",Transform_elem(*lt,0,i));
      print ( "\n" );
      print ( "                      ");
      for_less(i,0,4) print ("%9.4f ",Transform_elem(*lt,1,i));
      print ( "\n" );
      print ( "                      ");
      for_less(i,0,4) print ("%9.4f ",Transform_elem(*lt,2,i));
      print ( "\n" );
    }

    if (main_args.trans_info.center[0]!=-DBL_MAX ||
	main_args.trans_info.center[1]!=-DBL_MAX ||
	main_args.trans_info.center[2]!=-DBL_MAX) {
      print ( "Transform center   = %8.3f %8.3f %8.3f\n", 
	      main_args.trans_info.center[0],
	      main_args.trans_info.center[1],
	      main_args.trans_info.center[2] );
    }
    else {
      print ( "Transform center   = %8.3f %8.3f %8.3f\n", 0.0,0.0,0.0);
    }
  
    print ( "Transform trans    = %8.3f %8.3f %8.3f\n", 
		main_args.trans_info.translations[0],
		main_args.trans_info.translations[1],
		main_args.trans_info.translations[2] );
    print ( "Transform rotation = %8.3f %8.3f %8.3f\n", 
		main_args.trans_info.rotations[0],
		main_args.trans_info.rotations[1],
		main_args.trans_info.rotations[2] );
    print ( "Transform scale    = %8.3f %8.3f %8.3f\n\n", 
		main_args.trans_info.scales[0],
		main_args.trans_info.scales[1],
		main_args.trans_info.scales[2] );


  }
  

  ALLOC( data, 1 );		/* read in source data and target model */
  ALLOC( model, 1 );
  ALLOC(data_dx, 1); ALLOC(data_dy,   1);
  ALLOC(data_dz, 1); ALLOC(data_dxyz, 1);
  ALLOC(model_dx,1); ALLOC(model_dy,  1);
  ALLOC(model_dz,1); ALLOC(model_dxyz,1);
  
  if (main_args.trans_info.use_magnitude) {

    /* non-linear optimization is based on the correlation of a local
       sub-lattice between source and target gradient magnitude volumes */

    if (main_args.trans_info.transform_type == TRANS_NONLIN)
      DEBUG_PRINT( "This run will use sub-lattice correlation between the two input vols.\n");


    status = input_volume( main_args.filenames.data, 3, default_dim_names, 
			  NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			  TRUE, &data, (minc_input_options *)NULL );
    if (status != OK)
      print_error_and_line_num("Cannot input volume '%s'",
			       __FILE__, __LINE__,main_args.filenames.data);
    data_dxyz = data;
    
    status = input_volume( main_args.filenames.model, 3, default_dim_names, 
			  NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			  TRUE, &model, (minc_input_options *)NULL );
    if (status != OK)
      print_error_and_line_num("Cannot input volume '%s'",
			       __FILE__, __LINE__,main_args.filenames.model);
    model_dxyz = model;
  }
  else {

    /* non-linear optimization is based on the correlation of the
       Taylor projections stored in the {blurred,dx,dy,dz} data
       volumes, so load all the files based on the basename */

    if (main_args.trans_info.transform_type == TRANS_NONLIN)
      DEBUG_PRINT( "This run will use projection correlation between the 10 input vols.\n");
    status = read_all_data(&data, &data_dx, &data_dy, &data_dz, &data_dxyz,
			   main_args.filenames.data );
    if (status!=OK)
      print_error_and_line_num("Cannot input gradient volumes for  '%s'",
		  __FILE__, __LINE__,main_args.filenames.data);

    status = read_all_data(&model, &model_dx, &model_dy, &model_dz, &model_dxyz,
			   main_args.filenames.model );
    if (status!=OK)
      print_error_and_line_num("Cannot input gradient volumes for  '%s'",
		  __FILE__, __LINE__,main_args.filenames.model);
  }

  get_volume_separations(data, step);
  get_volume_sizes(data, sizes);
  get_volume_real_range(data, &min_value, &max_value);
  DEBUG_PRINT3 ( "Source volume size: %3d  by %3d  by %d \n",
		sizes[X], sizes[Y], sizes[Z]);
  DEBUG_PRINT3 ( "Source voxel = %8.3f %8.3f %8.3f\n", 
		step[X], step[Y], step[Z]);
  DEBUG_PRINT2 ( "min/max value= %8.3f %8.3f\n", min_value, max_value);
  get_volume_voxel_range(data, &min_value, &max_value);
  DEBUG_PRINT2 ( "min/max voxel= %8.3f %8.3f\n\n", min_value, max_value);

  get_volume_separations(model, step);
  get_volume_sizes(model, sizes);
  get_volume_real_range(model, &min_value, &max_value);
  DEBUG_PRINT3 ( "Target volume size: %3d  by %3d  by %d \n",
		sizes[X], sizes[Y], sizes[Z]);
  DEBUG_PRINT3 ( "Target voxel = %8.3f %8.3f %8.3f\n", 
		step[X], step[Y], step[Z]);
  DEBUG_PRINT2 ( "min/max value= %8.3f %8.3f\n", min_value, max_value);
  get_volume_voxel_range(model, &min_value, &max_value);
  DEBUG_PRINT2 ( "min/max voxel= %8.3f %8.3f\n\n\n", min_value, max_value);
  
  if (get_volume_n_dimensions(data)!=3) {
    print_error_and_line_num ("Data file %s has %d dimensions.  Only 3 dims supported.", 
			      __FILE__, __LINE__, main_args.filenames.data, 
			      get_volume_n_dimensions(data));
  }
  if (get_volume_n_dimensions(model)!=3) {
    print_error_and_line_num ("Model file %s has %d dimensions.  Only 3 dims supported.", 
			      __FILE__, __LINE__, main_args.filenames.model, 
			      get_volume_n_dimensions(model));
  }


  /* ===========================  translate initial transformation matrix into 
                                  transformation parameters */

  if (!init_params( data, model, mask_data, mask_model, &main_args )) {
    print_error_and_line_num("%s",__FILE__, __LINE__,
			     "Could not initialize transformation parameters\n");
  }

  if (strlen(main_args.filenames.matlab_file) != 0) {
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
		main_args.trans_info.rotations[0]*180/3.1415927,
		main_args.trans_info.rotations[1]*180/3.1415927,
		main_args.trans_info.rotations[2]*180/3.1415927 );
  DEBUG_PRINT3 ( "Transform scale    = %8.3f %8.3f %8.3f\n", 
		main_args.trans_info.scales[0],
		main_args.trans_info.scales[1],
		main_args.trans_info.scales[2] );
  DEBUG_PRINT3 ( "Transform shear    = %8.3f %8.3f %8.3f\n\n", 
		main_args.trans_info.shears[0],
		main_args.trans_info.shears[1],
		main_args.trans_info.shears[2] );
  

  if (strlen(main_args.filenames.measure_file) != 0) {

#include "measure_code.c"

  }

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


				/* calculate the actual transformation now. */

    if (main_args.trans_info.transform_type == TRANS_NONLIN) {

      build_default_deformation_field(&main_args);
      

      if ( !optimize_non_linear_transformation(data, data_dx, data_dy, data_dz, data_dxyz, 
					       model, model_dx, model_dy, model_dz, model_dxyz, 
					       mask_data, mask_model, &main_args ) ) {
	print_error_and_line_num("Error in optimization of non-linear transformation\n",
				 __FILE__, __LINE__);
	exit(EXIT_FAILURE);
      }
      
    }
    else {
      

      if (!optimize_linear_transformation( data, model, mask_data, mask_model, &main_args )) {
	print_error_and_line_num("Error in optimization of linear transformation\n",
				 __FILE__, __LINE__);
	exit(EXIT_FAILURE);
      }
      
    }

    if (number_dimensions==3) {
      print ("Initial correlation val = %0.5f\n",initial_corr); 
      print ("Final correlation value = %0.5f\n",
           xcorr_objective(data, model, mask_data, mask_model, &main_args ));
    }

  }


			/* if I have internally inverted the transform,
			   than flip it back forward before the save.   */

  if (main_args.trans_info.invert_mapping_flag) {

    DEBUG_PRINT ("Re-inverting transformation\n");
    create_inverse_general_transform(main_args.trans_info.transformation,
				     &tmp_invert);

    copy_general_transform(&tmp_invert,main_args.trans_info.transformation);

  }

  /* ===========================   write out transformation =============== */

  status = output_transform_file(main_args.filenames.output_trans,
				 comments,
				 main_args.trans_info.transformation);
  if (status!=OK) {
    print_error_and_line_num("Error saving transformation file.`\n",
		__FILE__, __LINE__);
    exit(EXIT_FAILURE);
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
/* ARGSUSED */
public int get_transformation(char *dst, char *key, char *nextArg)
{     
   Program_Transformation *transform_info;
   General_transform *transformation;
   General_transform input_transformation;
   FILE *fp;
   int ch, index;
   Status status;

   /* Check for following argument */
   if (nextArg == NULL) {
      (void)fprintf(stderr, 
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
         (void)fprintf(stderr, "Error opening temporary file.\n");
         exit(EXIT_FAILURE);
      }
      while ((ch=getc(stdin))!=EOF) (void) putc(ch, fp);
      rewind(fp);
   }
   else {
     status = open_file_with_default_suffix( nextArg,
					    get_default_transform_file_suffix(),
					    READ_FILE, ASCII_FORMAT, &fp );

      if (status != OK) {
         (void)fprintf(stderr, "Error opening transformation file %s.\n",
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
   if (input_transform(fp, nextArg, &input_transformation)!=OK) {
     (void)fprintf(stderr, "Error reading transformation file.\n");
     exit(EXIT_FAILURE);
   }
   (void) fclose(fp);



   copy_general_transform(&input_transformation, transformation);

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
@CREATED    : Wed May 26 13:05:44 EST 1993 Louis Collins
@MODIFIED   : Wed Jun 16 13:21:18 EST 1993 LC
    added: set of main_args.filenames.mask_model and .mask_data

---------------------------------------------------------------------------- */
/* ARGSUSED */
public int get_mask_file(char *dst, char *key, char *nextArg)
{ 

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
    (void)fprintf(stderr, "Cannot input mask file %s.",nextArg);
    return(FALSE);
  } 

  return TRUE;
  
}


