/* ----------------------------- MNI Header -----------------------------------
   @NAME       : minctracclib.c
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

	This file was created to support the librification of minctracc. At the
	time of creation it contains a library version of minctracc.c, leaving
	that file as a stub.

   @CREATED    : February 15, 2013 - robert brown
---------------------------------------------------------------------------- */

#include <config.h>
#include <float.h>
#include <volume_io.h>
#include <minctracc.h>
#include <objectives.h>
#include "local_macros.h"
#include "globaldefs.h"
#include "extras.h"



/* objective function for nonlinear optimization.
 * Set in get_nonlinear_objective().
 */
static int obj_func0 = -1;

static char *default_dim_names[VIO_N_DIMENSIONS] = 
    { MIzspace, MIyspace, MIxspace };



VIO_General_transform* minctracc( VIO_Volume source, VIO_Volume target, VIO_Volume sourceMask, VIO_Volume targetMask, VIO_General_transform *initialXFM, int iterations, float weight, float simplexSize, float stiffness, float similarity, float sub_lattice, Arg_Data *args) {

	VIO_Transform identityTransform;
	VIO_General_transform *origTransform = NULL;
	VIO_Status status;
	VIO_General_transform tmp_invert;
	
	int parse_flag, measure_matlab_flag, sizes[3],i,num_features;
	VIO_Real min_value, max_value, step[3];
  
 
	// Kill this (and main_args) to test whether the thread safety is working
	data = source; model = target;
	mask_data = sourceMask; mask_model = targetMask;
	main_args = args;
	iteration_limit = iterations;
	iteration_weight = weight;
	simplex_size = simplexSize;
	smoothing_weight = stiffness;
	similarity_cost_ratio = similarity;
	Diameter_of_local_lattice = sub_lattice;
			
	
	// SET UP INPUT TRANSFORMATIONS
	if (initialXFM) {
		args->trans_info.orig_transformation = initialXFM;
		args->trans_info.use_default = FALSE;
	}
	else {
		make_identity_transform(&identityTransform);
		ALLOC(origTransform,1);
	    create_linear_transform(origTransform,&identityTransform);
		args->trans_info.orig_transformation = origTransform;
	} 

	if (args->trans_info.use_identity)
		args->trans_info.use_default = FALSE;

	ALLOC(args->trans_info.transformation,1);
	copy_general_transform(args->trans_info.orig_transformation, args->trans_info.transformation);	
		
	
	switch (args->interpolant_type) {
	case TRICUBIC:
		args->interpolant = tricubic_interpolant;
		break;
	case TRILINEAR:
		args->interpolant = trilinear_interpolant;
		break;
	case N_NEIGHBOUR:
		args->interpolant = nearest_neighbour_interpolant;
		break;
	default:
		(void) fprintf(stderr, "Error determining interpolation type: %d\n",args->interpolant_type);
		return NULL;
	}

	switch (args->obj_function_type) {
	case XCORR:
		args->obj_function = xcorr_objective;
		break;
	case ZSCORE:
		args->obj_function = zscore_objective;
		break;
	case SSC:
		args->obj_function = ssc_objective;
		break;
	case VR:
		args->obj_function = vr_objective;
		break;
	case MUTUAL_INFORMATION:
		args->obj_function = mutual_information_objective;
		break;
	case NORMALIZED_MUTUAL_INFORMATION:
		args->obj_function = normalized_mutual_information_objective;
		break;
	default:
		(void) fprintf(stderr, "Error determining objective function type\n");
		return NULL;
	}
	
	
	get_volume_separations(data, step);
	get_volume_sizes(data, sizes);
	get_volume_minimum_maximum_real_value(data, &min_value, &max_value);
	get_volume_voxel_range(data, &min_value, &max_value);
	get_volume_separations(model, step); 
	get_volume_sizes(model, sizes);
	get_volume_minimum_maximum_real_value(model, &min_value, &max_value);
	get_volume_voxel_range(model, &min_value, &max_value);

	if (!init_params( data, model, mask_data, mask_model, args )) {
		print_error_and_line_num("%s",__FILE__, __LINE__,"Could not initialize transformation parameters\n");
	}
	
	if (args->features.number_of_features == 0) {
		num_features = allocate_a_new_feature(&(args->features));
		args->features.data[0]            = source; 
		args->features.model[0]           = target;
		args->features.data_name[0]       = args->filenames.data;
		args->features.model_name[0]      = args->filenames.model;
		args->features.data_mask[0]       = sourceMask;
		args->features.model_mask[0]      = targetMask;
		args->features.mask_data_name[0]  = args->filenames.mask_data;
		args->features.mask_model_name[0] = args->filenames.mask_model;
		args->features.thresh_data[0]     = args->threshold[0];
		args->features.thresh_model[0]    = args->threshold[1];
		if (args->trans_info.use_magnitude) {
	/*		switch (args->obj_function_type) {
				case NONLIN_XCORR:
					args->features.obj_func[0] = NONLIN_XCORR; break;
				case NONLIN_DIFF:
					args->features.obj_func[0] = NONLIN_DIFF; break;
				case NONLIN_SQDIFF:
					args->features.obj_func[0] = NONLIN_SQDIFF; break;
				case NONLIN_LABEL:
					args->features.obj_func[0] = NONLIN_LABEL; break;
				case NONLIN_CHAMFER:
					args->features.obj_func[0] = NONLIN_CHAMFER; break;
				case NONLIN_OPTICALFLOW:
					args->features.obj_func[0] = NONLIN_OPTICALFLOW; break;
				case NONLIN_CORRCOEFF:
					args->features.obj_func[0] = NONLIN_CORRCOEFF; break;
				default:
					args->features.obj_func[0] = NONLIN_XCORR; 
			} */
			args->features.obj_func[0] = args->obj_function_type;
		} 
		else {
			args->features.obj_func[0]        = NONLIN_OPTICALFLOW;    
		}
		args->features.weight[0]          = 1.0; 
	}
	
	// Go!
	if (args->trans_info.transform_type != TRANS_PAT) {
		init_lattice( data, model, mask_data, mask_model, args );

		if (args->trans_info.transform_type == TRANS_NONLIN) {
			build_default_deformation_field(args);
			if ( !optimize_non_linear_transformation(args) ) {
				print_error_and_line_num("Error in optimization of non-linear transformation\n", __FILE__, __LINE__);
			}
		}
		else {
			if (args->trans_info.rotation_type == TRANS_ROT ) {
				
				if (!optimize_linear_transformation( data, model, mask_data, mask_model, args )) {
					print_error_and_line_num("Error in optimization of linear transformation\n",__FILE__, __LINE__);
				}
			}

			if (args->trans_info.rotation_type == TRANS_QUAT ) {
				if (!optimize_linear_transformation_quater( data, model, mask_data, mask_model, args )) {
					print_error_and_line_num("Error in optimization of linear transformation\n",__FILE__, __LINE__);
				}
			}
		}
	}

	if (args->flags.verbose>0) {
		print ("Initial objective function val = %0.8f\n",initial_corr); 
		print ("Final objective function value = %0.8f\n",final_corr);
	}

	// if I have internally inverted the transform, then flip it back forward before the save.

	if (args->trans_info.invert_mapping_flag) {
		create_inverse_general_transform(args->trans_info.transformation,&tmp_invert);
		copy_general_transform(&tmp_invert,args->trans_info.transformation);
	}
	
	
	if (origTransform) FREE(origTransform);
	return( args->trans_info.transformation );
	
}
  
	
void initializeArgs(Arg_Data *args) {
	// Program filenames
	args->filenames.data = "";
	args->filenames.model = "";
	args->filenames.mask_data = "";
	args->filenames.mask_model = "";
	args->filenames.output_trans = "";
	args->filenames.measure_file = "";
	args->filenames.matlab_file = "";
	
	// Program flags
	args->flags.verbose = 0; args->flags.debug = FALSE;

	// Transformation flags
	args->trans_info.use_identity = FALSE;
	args->trans_info.use_default = TRUE;
	args->trans_info.use_magnitude = TRUE;
	args->trans_info.max_def_magnitude = 50.0;
	args->trans_info.use_simplex = TRUE;
        args->trans_info.use_bfgs = FALSE;
	args->trans_info.use_super = 2;
	args->trans_info.use_local_smoothing = FALSE;
	args->trans_info.use_local_isotropic = TRUE;
	args->trans_info.file_name = "";
	args->trans_info.file_contents = NULL;
	args->trans_info.buffer_length = 0,
	args->trans_info.transformation = (VIO_General_transform *)NULL,
	args->trans_info.orig_transformation = (VIO_General_transform *)NULL,
	args->trans_info.transform_type = TRANS_LSQ7;
	args->trans_info.center[0] = -DBL_MAX; 	args->trans_info.center[1] = -DBL_MAX; 	args->trans_info.center[2] = -DBL_MAX;
	args->trans_info.scales[0] = 1.0; args->trans_info.scales[1] = 1.0; args->trans_info.scales[2] = 1.0;
	args->trans_info.shears[0] = 0.0; args->trans_info.shears[1] = 0.0; args->trans_info.shears[2] = 0.0; 
	args->trans_info.translations[0] = 0.0; args->trans_info.translations[1] = 0.0; args->trans_info.translations[2] = 0.0;
	args->trans_info.quaternions[0] = 0.0; args->trans_info.quaternions[1] = 0.0; args->trans_info.quaternions[2] = 0.0; args->trans_info.quaternions[3] = 1.0;
	args->trans_info.rotations[0] = 0.0; args->trans_info.rotations[1] = 0.0; args->trans_info.rotations[2] = 0.0; 
	args->trans_info.weights[0] = 1.0; args->trans_info.weights[1] = 1.0; args->trans_info.weights[2] = 1.0;
	args->trans_info.weights[3] = 3.1415927/180.0; args->trans_info.weights[4] = 3.1415927/180.0; args->trans_info.weights[5] = 3.1415927/180.0; 
	args->trans_info.weights[6] = 0.02; args->trans_info.weights[7] = 0.02; args->trans_info.weights[8] = 0.02; 
	args->trans_info.weights[9] = 0.02; args->trans_info.weights[10] = 0.02; args->trans_info.weights[11] = 0.02; 
	args->trans_info.invert_mapping_flag = FALSE;
	args->trans_info.rotation_type =	TRANS_ROT;

	// Feature volumes
	args->features.number_of_features = 0;
	args->features.data = NULL; args->features.model = NULL; args->features.data_mask = NULL; args->features.model_mask = NULL;
	
	// Flags
	args->interpolant = trilinear_interpolant;
	args->interpolant_type = TRILINEAR;
	args->obj_function = xcorr_objective;
	args->obj_function_type = XCORR;
	args->optimize_type = OPT_SIMPLEX;
	args->force_lattice = 0;
	args->step[0] = 4.0; args->step[1] = 4.0; args->step[2] = 4.0; 
	args->lattice_width[0] = 24.0; args->lattice_width[1] = 24.0; args->lattice_width[2] = 24.0;
	args->start[0] = 0.0; args->start[1] = 0.0; args->start[2] = 0.0; 
	args->count[0] = 0; args->count[1] = 0; args->count[2] = 0;
	
	args->directions[0].coords[0] = 1.0; args->directions[0].coords[1] = 0.0; args->directions[0].coords[2] = 0.0; 
	args->directions[1].coords[0] = 0.0; args->directions[1].coords[1] = 1.0; args->directions[1].coords[2] = 0.0; 
	args->directions[2].coords[0] = 0.0; args->directions[2].coords[1] = 0.0; args->directions[2].coords[2] = 1.0; 
	args->smallest_vol = 1;
	
	// Estimation flags
	args->trans_flags.estimate_center = FALSE; args->trans_flags.estimate_scale = FALSE; 
	args->trans_flags.estimate_trans = FALSE; args->trans_flags.estimate_rots = FALSE; args->trans_flags.estimate_quats = FALSE; 
	
	// More flags
	args->threshold[0] = 0.0; args->threshold[1] = 0.0;
	args->speckle = 5.0;
	args->groups = 256;
	args->blur_pdf = 3;	
}

/* Command line argument "-nonlinear" may be followed by an optional
 * string e.g. "xcorr" to specify the objective function.  If no
 * objective is specified, the default is xcorr.
 *
 * If nextArg is used to specify the objective function, return 1
 * to inform ParseArgv to skip that argument; else return 0.
 */
int get_nonlinear_objective(char *dst, char *key, char* nextArg)
{
    main_args->trans_info.transform_type = TRANS_NONLIN;

    if (nextArg == NULL) {
        obj_func0 = NONLIN_XCORR;
        return 0;
    }

    if (strcmp( "xcorr", nextArg ) == 0 ) {
        obj_func0 = NONLIN_XCORR;
    } else if (strcmp( "diff", nextArg ) == 0 ) {
        obj_func0 = NONLIN_DIFF;
    } else if (strcmp( "sqdiff", nextArg ) == 0 ) {
        obj_func0 = NONLIN_SQDIFF;
    } else if (strcmp( "label", nextArg ) == 0 ) {
        obj_func0 = NONLIN_LABEL;
    } else if (strcmp( "chamfer", nextArg ) == 0 ) {
        obj_func0 = NONLIN_CHAMFER;
    } else if (strcmp( "opticalflow", nextArg ) == 0 ) {
        obj_func0 = NONLIN_OPTICALFLOW;
    } else if (strcmp( "corrcoeff", nextArg ) == 0 ) {
        obj_func0 = NONLIN_CORRCOEFF;
    } else {
        obj_func0 = NONLIN_XCORR;
        return 0;
    }

    return 1;
}


int free_features(Feature_volumes *features)
{


  if (*(features->data)       != (VIO_Volume)NULL) {delete_volume(*(features->data));      } FREE(features->data); 
  if (*(features->model)      != (VIO_Volume)NULL) {delete_volume(*(features->model));     } FREE(features->model); 

  if (*(features->data_mask)  != (VIO_Volume)NULL) {delete_volume(*(features->data_mask)); } FREE(features->data_mask); 
  if (*(features->model_mask) != (VIO_Volume)NULL) {delete_volume(*(features->model_mask));} FREE(features->model_mask); 

  FREE(features->data_name);
  FREE(features->model_name);
  FREE(features->mask_data_name);
  FREE(features->mask_model_name);
  FREE(features->obj_func);
  FREE(features->weight);
  FREE(features->thresh_data);
  FREE(features->thresh_model);

}


int allocate_a_new_feature(Feature_volumes *features)
{

  int i;

  i = main_args->features.number_of_features;
  main_args->features.number_of_features++;    

  if (i==0) {

    ALLOC(features->data,1);
    ALLOC(features->model,1); 
    ALLOC(features->data_name, 1);
    ALLOC(features->model_name, 1);
    ALLOC(features->data_mask,1);
    ALLOC(features->model_mask,1);
    ALLOC(features->mask_data_name, 1);
    ALLOC(features->mask_model_name, 1);
    ALLOC(features->obj_func, 1);
    ALLOC(features->weight, 1);
    ALLOC(features->thresh_data, 1);
    ALLOC(features->thresh_model, 1);
  }
  else {
    REALLOC(features->data,i+1);
    REALLOC(features->model,i+1);
    REALLOC(features->data_name, i+1);
    REALLOC(features->model_name, i+1);
    REALLOC(features->data_mask,i+1);
    REALLOC(features->model_mask,i+1);
    REALLOC(features->mask_data_name, i+1);
    REALLOC(features->mask_model_name, i+1);
    REALLOC(features->obj_func, i+1);
    REALLOC(features->weight, i+1);
    REALLOC(features->thresh_data, i+1);
    REALLOC(features->thresh_model, i+1);
  }
  return(i);
}


void add_a_feature_for_matching(Feature_volumes *features,
                                VIO_Volume data_vol,
                                VIO_Volume model_vol,
                                VIO_Volume data_mask,
                                VIO_Volume model_mask,
                                char *data_name,
                                char *model_name,
                                char *mask_data_name,
                                char *mask_model_name,
                                char obj_func,
                                VIO_Real weight,
                                VIO_Real thresh_data,
                                VIO_Real thresh_model)
{

  int i;

  i = allocate_a_new_feature(features);

  features->data[i]            = data_vol; 
  features->model[i]           = model_vol;
  features->data_name[i]       = data_name;
  features->model_name[i]      = model_name;
  features->data_mask[i]       = data_mask;
  features->model_mask[i]      = model_mask;
  features->mask_data_name[i]  = mask_data_name;
  features->mask_model_name[i] = mask_model_name;
  features->thresh_data[i]     = thresh_data;
  features->thresh_model[i]    = thresh_model;
  features->obj_func[i]        = obj_func;
  features->weight[i]          = weight;

}




/* MINCTRACCOLDFASHIONED keeps the same interface as the old minctracc main().  The new minctracc
  main() forwards its argc and argv to this function.  Eventually maybe I'll fix this to feed through
  the new minctracc function.
*/

int minctraccOldFashioned ( int argc, char* argv[] )
{
  VIO_Status 
    status;
  VIO_General_transform
    tmp_invert;
  VIO_Transform 
    *lt, ident_trans;
  int
    parse_flag,
    measure_matlab_flag,
    
    sizes[3],i,num_features;
  VIO_Real
    min_value, max_value, step[3];
  char 
    *comments = history_string( argc, argv );
  FILE
    *ofd;
  VIO_Real
    obj_func_val;
  float quat4;
  
  prog_name     = argv[0];        
  
  /* Call ParseArgv to interpret all command line args (returns TRUE if error) */

  parse_flag = ParseArgv(&argc, argv, argTable, 0);

  measure_matlab_flag = 
    (strlen(main_args->filenames.matlab_file)  != 0) ||
    (strlen(main_args->filenames.measure_file) != 0);

  /* assign objective function and interpolant type */
  switch (main_args->interpolant_type) {
  case TRICUBIC:
    main_args->interpolant = tricubic_interpolant;
    break;
  case TRILINEAR:
    main_args->interpolant = trilinear_interpolant;
    break;
  case N_NEIGHBOUR:
    main_args->interpolant = nearest_neighbour_interpolant;
    break;
  default:
    (void) fprintf(stderr, "Error determining interpolation type\n");
    exit(EXIT_FAILURE);
  }

  switch (main_args->obj_function_type) {
  case XCORR:
    main_args->obj_function = xcorr_objective;
    break;
  case ZSCORE:
    main_args->obj_function = zscore_objective;
    break;
  case SSC:
    main_args->obj_function = ssc_objective;
    break;
  case VR:
    main_args->obj_function = vr_objective;
    break;
  case MUTUAL_INFORMATION:
    main_args->obj_function = mutual_information_objective;
    break;
  case NORMALIZED_MUTUAL_INFORMATION:
    main_args->obj_function = normalized_mutual_information_objective;
    break;
  default:
    (void) fprintf(stderr, "Error determining objective function type\n");
    exit(EXIT_FAILURE);
  }
  
  if(main_args->trans_info.use_bfgs)
    main_args->optimize_type=OPT_BFGS;

  if (parse_flag || 
      (measure_matlab_flag && argc!=3) ||
      (!measure_matlab_flag && argc!=4)) {

    print ("Parameters left:\n");
    for(i=0; i<argc; i++)
      print ("%s ",argv[i]);
    print ("\n");
    
    (void)fprintf(stderr, 
                  "\nUsage: %s [<options>] <sourcefile> <targetfile> <output transfile>\n", 
                  prog_name);
    (void)fprintf(stderr,"       %s [-help]\n\n", prog_name);


                                /* helpful hints */
    if ((argc == 4) && 
        ((strlen(main_args->filenames.matlab_file)  != 0) ||
         (strlen(main_args->filenames.measure_file) != 0))) {
      (void)fprintf(stderr, "\nNote: No output transform file needs to be specified for -matlab\n");
      (void)fprintf(stderr, "      or -measure options.\n");
    }

    exit(EXIT_FAILURE);
  }

  main_args->filenames.data  = argv[1];        /* set up necessary file names */
  main_args->filenames.model = argv[2];
  if (strlen(main_args->filenames.measure_file)==0 &&
      strlen(main_args->filenames.matlab_file)==0) 
    main_args->filenames.output_trans = argv[3];


                                /* check to see if they can be overwritten */
  if (!clobber_flag && 
      (strlen(main_args->filenames.measure_file)!=0) && 
      file_exists(main_args->filenames.measure_file)) {
    (void)fprintf (stderr,"Measure file %s exists.\n",main_args->filenames.measure_file);
    (void)fprintf (stderr,"Use -clobber to overwrite.\n");
    exit(EXIT_FAILURE);
  }

  if (!clobber_flag && 
      (strlen(main_args->filenames.matlab_file)!=0) && 
      file_exists(main_args->filenames.matlab_file)) {
    (void)fprintf (stderr,"Matlab file %s exists.\n",main_args->filenames.matlab_file);
    (void)fprintf (stderr,"Use -clobber to overwrite.\n");
    exit(EXIT_FAILURE);
  }

  if (!clobber_flag && 
      (strlen(main_args->filenames.output_trans)!=0) && 
      file_exists(main_args->filenames.output_trans)) {
    (void)fprintf (stderr,"Output file %s exists.\n",main_args->filenames.output_trans);
    (void)fprintf (stderr,"Use -clobber to overwrite.\n");
    exit(EXIT_FAILURE);
  }
  if (strlen(main_args->filenames.matlab_file)  != 0 &&
      strlen(main_args->filenames.measure_file) != 0) {
    (void)fprintf(stderr, "\nWARNING: -matlab and -measure are mutually exclusive.  Only\n");
    (void)fprintf(stderr, "          the matlab option will be executed here.\n");
  }

                                /* set up linear transformation to identity, if
                                   not set up in the argument list */

  if (main_args->trans_info.transformation == (VIO_General_transform *)NULL) { 
                                /* build identity VIO_General_transform */
    make_identity_transform(&ident_trans);
                                /* attach it to global w-to-w transformation */
    ALLOC(main_args->trans_info.transformation,1);
    create_linear_transform(main_args->trans_info.transformation,&ident_trans);    
  } 

                                /* don't use default PAT if user req -identity */
  if (main_args->trans_info.use_identity)
    main_args->trans_info.use_default = FALSE;


                                /* make a copy of the original transformation */

  ALLOC(main_args->trans_info.orig_transformation,1);
  copy_general_transform(main_args->trans_info.transformation,
                         main_args->trans_info.orig_transformation);




  if (main_args->flags.debug) {
    /*
       this is simply to print out debugging info at the beginning of a run
    */

    print ( "===== Debugging information from %s =====\n", prog_name);
    print ( "Data filename       = %s\n", main_args->filenames.data);
    print ( "Model filename      = %s\n", main_args->filenames.model);
    print ( "Data mask filename  = %s\n", main_args->filenames.mask_data);
    print ( "Model mask filename = %s\n", main_args->filenames.mask_model);
    print ( "Input xform name    = %s\n", main_args->trans_info.file_name);
    if (strlen(main_args->filenames.output_trans) != 0)
      print ( "Output filename     = %s\n", main_args->filenames.output_trans);
    if (strlen(main_args->filenames.matlab_file)  != 0)
      print ( "Matlab filename     = %s (num_steps=%d)\n\n",main_args->filenames.matlab_file,Matlab_num_steps );
    if (strlen(main_args->filenames.measure_file) != 0)
      print ( "Measure filename    = %s\n\n",main_args->filenames.measure_file );
    print ( "Step size           = %f %f %f\n",
                  main_args->step[0],
                  main_args->step[1],
                  main_args->step[2]);
    print ( "Sub-lattice dia     = %f %f %f\n",
                  main_args->lattice_width[0],
                  main_args->lattice_width[1],
                  main_args->lattice_width[2]);
    print  ( "Objective function  = ");
    if (main_args->obj_function == xcorr_objective) {
      print("cross correlation (threshold = %f %f)\n",
                   main_args->threshold[0],main_args->threshold[1]);
    }
    else
      if (main_args->obj_function == zscore_objective) {
        print( "zscore (threshold = %f %f)\n", main_args->threshold[0],main_args->threshold[1] );
      }
      else
        if (main_args->obj_function == vr_objective) {
          print( "ratio of variance (threshold = %f %f, groups = %d)\n", 
                       main_args->threshold[0],main_args->threshold[1], main_args->groups);
        }
        else
          if (main_args->obj_function == ssc_objective) {
            print( "stochastic sign change (threshold = %f %f, speckle %f %%)\n",
                         main_args->threshold[0],main_args->threshold[1], main_args->speckle);
          }
          else
            if (main_args->obj_function == mutual_information_objective || main_args->obj_function == normalized_mutual_information_objective ) {
              print( "mutual information (groups = %d) \n",
                           main_args->groups);
            }
            else {
              print("unknown!\n");
              (void)fprintf(stderr,"Unknown objective function requested.\n");
              exit(EXIT_FAILURE);
            }


    print ( "Transform linear    = %s\n", 
           (get_transform_type(main_args->trans_info.transformation)==LINEAR ? 
            "TRUE" : "FALSE") );
    print ( "Transform inverted? = %s\n", (main_args->trans_info.invert_mapping_flag ? 
                                           "TRUE" : "FALSE") );
    print ( "Transform type      = %d\n", main_args->trans_info.transform_type );

    if (get_transform_type(main_args->trans_info.transformation)==LINEAR) {
      lt = get_linear_transform_ptr(main_args->trans_info.transformation);
      print ( "Transform matrix    = ");
      for(i=0; i<4; i++) print ("%9.4f ",Transform_elem(*lt,0,i));
      print ( "\n" );
      print ( "                      ");
      for(i=0; i<4; i++) print ("%9.4f ",Transform_elem(*lt,1,i));
      print ( "\n" );
      print ( "                      ");
      for(i=0; i<4; i++) print ("%9.4f ",Transform_elem(*lt,2,i));
      print ( "\n" );
    }

    if (main_args->trans_info.center[0]!=-DBL_MAX ||
        main_args->trans_info.center[1]!=-DBL_MAX ||
        main_args->trans_info.center[2]!=-DBL_MAX) {
      print ( "Transform center   = %8.3f %8.3f %8.3f\n", 
              main_args->trans_info.center[0],
              main_args->trans_info.center[1],
              main_args->trans_info.center[2] );
    }
    else 
      {
        print ( "Transform center   = %8.3f %8.3f %8.3f\n", 0.0,0.0,0.0);
      }

    /* rotation with quaternion or euler angle */
    if (main_args->trans_info.rotation_type == TRANS_QUAT)
      {
        quat4=sqrt(1-SQR(main_args->trans_info.quaternions[0])-SQR(main_args->trans_info.quaternions[1])-SQR(main_args->trans_info.quaternions[2]));
        print ( "Transform quaternion   = %8.3f %8.3f %8.3f %8.3f \n\n", 
                main_args->trans_info.quaternions[0],
                main_args->trans_info.quaternions[1],
                main_args->trans_info.quaternions[2],
                quat4 );
      }
    if (main_args->trans_info.rotation_type == TRANS_ROT)
      {
        print ( "Transform rotation   = %8.3f %8.3f %8.3f \n\n", 
                main_args->trans_info.rotations[0],
                main_args->trans_info.rotations[1],
                main_args->trans_info.rotations[2] );
      }
    print ( "Transform trans    = %8.3f %8.3f %8.3f\n", 
            main_args->trans_info.translations[0],
            main_args->trans_info.translations[1],
            main_args->trans_info.translations[2] );
    print ( "Transform scale    = %8.3f %8.3f %8.3f\n\n", 
            main_args->trans_info.scales[0],
            main_args->trans_info.scales[1],
            main_args->trans_info.scales[2] );
    

  }


  


  if (main_args->trans_info.use_magnitude) 
    {
      /* non-linear optimization is based on the correlation of a local
         sub-lattice between source and target gradient magnitude volumes 
         for the first two volumes */

      if (main_args->trans_info.transform_type == TRANS_NONLIN)
        DEBUG_PRINT1( "This run will use sub-lattice correlation (type %d) between the two input vols.\n", obj_func0);
    }
  else 
    {
      /* non-linear optimization is based on the direct computation of a
         deformation vector, based on optical flow.  Both volume _MUST_
         be intensity normalized! */

      if (main_args->trans_info.transform_type == TRANS_NONLIN)
        DEBUG_PRINT( "This run will use optical flow.\n");
    }

  ALLOC(data,1);

  status = input_volume( main_args->filenames.data, 3, default_dim_names, 
                         NC_DOUBLE, FALSE, 0.0, 0.0,
                         TRUE, &data, (minc_input_options *)NULL );

  if (status != VIO_OK)
    print_error_and_line_num("Cannot input volume '%s'",
                             __FILE__, __LINE__,main_args->filenames.data);
  data_dxyz = data;
 
  status = input_volume( main_args->filenames.model, 3, default_dim_names, 
                         NC_DOUBLE, FALSE, 0.0, 0.0,
                         TRUE, &model, (minc_input_options *)NULL );
  if (status != VIO_OK)
    print_error_and_line_num("Cannot input volume '%s'",
                             __FILE__, __LINE__,main_args->filenames.model);
  model_dxyz = model;
 

  get_volume_separations(data, step);
  get_volume_sizes(data, sizes);
  DEBUG_PRINT3 ( "Source volume size: %3d  by %3d  by %d \n",
                sizes[VIO_X], sizes[VIO_Y], sizes[VIO_Z]);
  DEBUG_PRINT3 ( "Source voxel size = %8.3f %8.3f %8.3f\n", 
                step[VIO_X], step[VIO_Y], step[VIO_Z]);

  get_volume_minimum_maximum_real_value(data, &min_value, &max_value);
  DEBUG_PRINT2 ( "Source min/max real range = %8.3f %8.3f\n", min_value, max_value);

  get_volume_voxel_range(data, &min_value, &max_value);
  DEBUG_PRINT2 ( "Source min/max voxel= %8.3f %8.3f\n\n", min_value, max_value);

  get_volume_separations(model, step); 
  get_volume_sizes(model, sizes);

  DEBUG_PRINT3 ( "Target volume size: %3d  by %3d  by %d \n",
                sizes[VIO_X], sizes[VIO_Y], sizes[VIO_Z]);
  DEBUG_PRINT3 ( "Target voxel = %8.3f %8.3f %8.3f\n", 
                step[VIO_X], step[VIO_Y], step[VIO_Z]);

  get_volume_minimum_maximum_real_value(model, &min_value, &max_value);
  DEBUG_PRINT2 ( "Target min/max real range= %8.3f %8.3f\n", min_value, max_value);

  get_volume_voxel_range(model, &min_value, &max_value);
  DEBUG_PRINT2 ( "Target min/max voxel = %8.3f %8.3f\n\n\n", min_value, max_value);
  
  if (get_volume_n_dimensions(data)!=3) 
    {
      print_error_and_line_num ("Data file %s has %d dimensions.  Only 3 dims supported.", 
                                __FILE__, __LINE__, main_args->filenames.data, 
                                get_volume_n_dimensions(data));
    }
  if (get_volume_n_dimensions(model)!=3) 
    {
      print_error_and_line_num ("Model file %s has %d dimensions.  Only 3 dims supported.", 
                                __FILE__, __LINE__, main_args->filenames.model, 
                                get_volume_n_dimensions(model));
    }


                                /* shift features to be able to
                                   insert the main source/target
                                   volumes first. */

  num_features = allocate_a_new_feature(&(main_args->features));
  for(i=num_features; i>=1; i--) {
    main_args->features.data[i]            = main_args->features.data[i-1];
    main_args->features.model[i]           = main_args->features.model[i-1];
    main_args->features.data_name[i]       = main_args->features.data_name[i-1];
    main_args->features.model_name[i]      = main_args->features.model_name[i-1];
    main_args->features.data_mask[i]       = main_args->features.data_mask[i-1];
    main_args->features.model_mask[i]      = main_args->features.model_mask[i-1];
    main_args->features.mask_data_name[i]  = main_args->features.mask_data_name[i-1];
    main_args->features.mask_model_name[i] = main_args->features.mask_model_name[i-1];
    main_args->features.thresh_data[i]     = main_args->features.thresh_data[i-1];
    main_args->features.thresh_model[i]    = main_args->features.thresh_model[i-1];
    main_args->features.obj_func[i]        = main_args->features.obj_func[i-1];
    main_args->features.weight[i]          = main_args->features.weight[i-1];
  }

  main_args->features.data[0]            = data; 
  main_args->features.model[0]           = model;
  main_args->features.data_name[0]       = main_args->filenames.data;
  main_args->features.model_name[0]      = main_args->filenames.model;
  main_args->features.data_mask[0]       = mask_data;
  main_args->features.model_mask[0]      = mask_model;
  main_args->features.mask_data_name[0]  = main_args->filenames.mask_data;
  main_args->features.mask_model_name[0] = main_args->filenames.mask_model;
  main_args->features.thresh_data[0]     = main_args->threshold[0];
  main_args->features.thresh_model[0]    = main_args->threshold[1];
  if (main_args->trans_info.use_magnitude) {
    main_args->features.obj_func[0]        = obj_func0;
  } 
  else {
    main_args->features.obj_func[0]        = NONLIN_OPTICALFLOW;    
  }
  main_args->features.weight[0]          = 1.0;

  /* ===========================  translate initial transformation matrix into 
                                  transformation parameters */

    if (!init_params( data, model, mask_data, mask_model, main_args )) {
      print_error_and_line_num("%s",__FILE__, __LINE__,
                             "Could not initialize transformation parameters\n");
    }



  if (strlen(main_args->filenames.matlab_file) != 0) 
    {
      init_lattice( data, model, mask_data, mask_model, main_args );
      make_matlab_data_file( data, model, mask_data, mask_model, comments, main_args );
      exit( VIO_OK );
    }

  if (strlen(main_args->filenames.measure_file) != 0) 
    {

#include "measure_code.c"

    /* measure code finishes with 
       exit(status); */
    }

  DEBUG_PRINT  ("AFTER init_params()\n");
  if (get_transform_type(main_args->trans_info.transformation)==LINEAR) {
    lt = get_linear_transform_ptr(main_args->trans_info.transformation);
    
    DEBUG_PRINT ( "Transform matrix    = ");
    for(i=0; i<4; i++) DEBUG_PRINT1 ("%9.4f ",Transform_elem(*lt,0,i));
    DEBUG_PRINT ( "\n" );
    DEBUG_PRINT ( "                      ");
    for(i=0; i<4; i++) DEBUG_PRINT1 ("%9.4f ",Transform_elem(*lt,1,i));
    DEBUG_PRINT ( "\n" );
    DEBUG_PRINT ( "                      ");
    for(i=0; i<4; i++) DEBUG_PRINT1 ("%9.4f ",Transform_elem(*lt,2,i));
    DEBUG_PRINT ( "\n" );
  }
  DEBUG_PRINT ( "\n" );
  
  DEBUG_PRINT3 ( "Transform center   = %8.3f %8.3f %8.3f\n", 
                main_args->trans_info.center[0],
                main_args->trans_info.center[1],
                main_args->trans_info.center[2] );
  if (main_args->trans_info.rotation_type == TRANS_ROT){
    DEBUG_PRINT3 ( "Transform rotations  = %8.3f %8.3f %8.3f\n", 
                   main_args->trans_info.rotations[0],
                   main_args->trans_info.rotations[1],
                   main_args->trans_info.rotations[2] );
  }
  if (main_args->trans_info.rotation_type == TRANS_QUAT)
    {
      quat4=sqrt(1-SQR(main_args->trans_info.quaternions[0])-SQR(main_args->trans_info.quaternions[1])-SQR(main_args->trans_info.quaternions[2]));
      DEBUG_PRINT4 ( "Transform quaternions  = %8.3f %8.3f %8.3f %8.3f\n", 
                     main_args->trans_info.quaternions[0],
                     main_args->trans_info.quaternions[1],
                     main_args->trans_info.quaternions[2],
                     quat4 );
    }
  

  DEBUG_PRINT3 ( "Transform trans    = %8.3f %8.3f %8.3f\n", 
                main_args->trans_info.translations[0],
                main_args->trans_info.translations[1],
                main_args->trans_info.translations[2] );
  DEBUG_PRINT3 ( "Transform scale    = %8.3f %8.3f %8.3f\n", 
                main_args->trans_info.scales[0],
                main_args->trans_info.scales[1],
                main_args->trans_info.scales[2] );
  DEBUG_PRINT3 ( "Transform shear    = %8.3f %8.3f %8.3f\n\n", 
                main_args->trans_info.shears[0],
                main_args->trans_info.shears[1],
                main_args->trans_info.shears[2] );
  


                                /* do not do any optimization if the transformation
                                   requested is the Principal Axes Transformation 
                                   then:
                                   =======   do linear fitting =============== */
  
  if (main_args->trans_info.transform_type != TRANS_PAT) {
    
                                /* initialize the sampling lattice and figure out
                                   which of the two volumes is smaller.           */
    
    init_lattice( data, model, mask_data, mask_model, main_args );

    if (main_args->smallest_vol == 1) {
      DEBUG_PRINT("Source volume is smallest\n");
    }
    else {
      DEBUG_PRINT("Target volume is smallest\n");
    }
    DEBUG_PRINT3 ( "Lattice step size  = %8.3f %8.3f %8.3f\n",
                  main_args->step[0],main_args->step[1],main_args->step[2]);
    DEBUG_PRINT3 ( "Lattice start      = %8.3f %8.3f %8.3f\n",
                  main_args->start[0],main_args->start[1],main_args->start[2]);
    DEBUG_PRINT3 ( "Lattice count      = %8d %8d %8d\n\n",
                  main_args->count[0],main_args->count[1],main_args->count[2]);


                                /* calculate the actual transformation now. */

    if (main_args->trans_info.transform_type == TRANS_NONLIN) {

      build_default_deformation_field(main_args);
      

      if ( !optimize_non_linear_transformation( main_args ) ) {
        print_error_and_line_num("Error in optimization of non-linear transformation\n",
                                 __FILE__, __LINE__);
        exit(EXIT_FAILURE);
      }
      
    }
    else {
      
      if (main_args->trans_info.rotation_type == TRANS_ROT )
        {
          if (!optimize_linear_transformation( data, model, mask_data, mask_model, main_args )) {
            print_error_and_line_num("Error in optimization of linear transformation\n",
                                 __FILE__, __LINE__);
        exit(EXIT_FAILURE);
          }
        }
      
      
      if (main_args->trans_info.rotation_type == TRANS_QUAT )
        {
          if (!optimize_linear_transformation_quater( data, model, mask_data, mask_model, main_args )) {
            print_error_and_line_num("Error in optimization of linear transformation\n",
                                 __FILE__, __LINE__);
        exit(EXIT_FAILURE);
          }
        }
      
      
    }

    if (number_dimensions==3 && main_args->flags.verbose>0) {
      print ("Initial objective function val = %0.8f\n",initial_corr); 
      print ("Final objective function value = %0.8f\n",final_corr);
    }

  }


                        /* if I have internally inverted the transform,
                           than flip it back forward before the save.   */

  if (main_args->trans_info.invert_mapping_flag) {

    DEBUG_PRINT ("Re-inverting transformation\n");
    create_inverse_general_transform(main_args->trans_info.transformation,
                                     &tmp_invert);

    copy_general_transform(&tmp_invert,main_args->trans_info.transformation);

  }




  /* ===========================   write out transformation =============== */


  status = output_transform_file(main_args->filenames.output_trans,
                                 comments,
                                 main_args->trans_info.transformation);
     
  if (status!=VIO_OK) {
    print_error_and_line_num("Error saving transformation file.`\n",
                __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  if( comments ) {
    FREE( comments );
    comments = NULL;
  }
  
  delete_general_transform(main_args->trans_info.transformation);
  FREE(main_args->trans_info.transformation);
  delete_general_transform(main_args->trans_info.orig_transformation);
  FREE(main_args->trans_info.orig_transformation);
  free_features (&(main_args->features));
  FREE(main_args->trans_info.file_contents);
  // Note: don't know how to free ALLOC(data) above. (Claude).

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
int get_transformation(char *dst, char *key, char *nextArg)
{     
   Program_Transformation *transform_info;
   VIO_General_transform *transformation;
   VIO_General_transform input_transformation;
   FILE *fp;
   int ch, index;
   VIO_Status status;

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

   if (transform_info->transformation == (VIO_General_transform *)NULL) 
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

      if (status != VIO_OK) {
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
   if (input_transform(fp, nextArg, &input_transformation)!=VIO_OK) {
     (void)fprintf(stderr, "Error reading transformation file.\n");
     exit(EXIT_FAILURE);
   }
   (void) fclose(fp);



   copy_general_transform(&input_transformation, transformation);

   delete_general_transform(&input_transformation);

   /* set a GLOBAL flag, to show that a transformation has been read in */

   main_args->trans_info.use_default = FALSE;

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
    added: set of main_args->filenames.mask_model and .mask_data

---------------------------------------------------------------------------- */
/* ARGSUSED */
int get_mask_file(char *dst, char *key, char *nextArg)
{ 

  VIO_Status status;

  if (strncmp ( "-model_mask", key, 2) == 0) {
    /*    ALLOC( mask_model, 1 );*/
    status = input_volume( nextArg, 3, default_dim_names, 
                          NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                          TRUE, &mask_model, (minc_input_options *)NULL );
    dst = nextArg;
    main_args->filenames.mask_model = nextArg;
  }
  else {
    /*    ALLOC( mask_data, 1);*/
    status = input_volume( nextArg, 3, default_dim_names, 
                          NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                          TRUE, &mask_data, (minc_input_options *)NULL );
    dst = nextArg;
    main_args->filenames.mask_data = nextArg;
  }

  if (status != VIO_OK)
  {
    (void)fprintf(stderr, "Cannot input mask file %s.",nextArg);
    return(FALSE);
  } 

  return TRUE;
  
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_feature volumes
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
    added: set of main_args->filenames.mask_model and .mask_data

---------------------------------------------------------------------------- */
/* ARGSUSED */
int get_feature_volumes(char *dst, char *key, int argc, char **argv)
{ 
  int 
    i,
    args_used,
    weight_index,
    obj_func_index;
  char 
    *end_ptr;
  VIO_Real
    tmp;
  VIO_Status status;

  VIO_Volume 
    data_vol,
    model_vol,
    data_mask,
    model_mask;
  char
    *data_name,
    *model_name,
    *mask_data_name,
    *mask_model_name,
    obj_func;
  VIO_Real 
    weight,
    thresh_data,
    thresh_model;


  if ( argc >=2 && argv[0] != NULL && argv[1] != NULL ) {

    data_name = argv[0];

    model_name = argv[1];

    obj_func        = NONLIN_XCORR;
    weight          = 1.0;
    thresh_data     = -DBL_MAX;
    thresh_model    = -DBL_MAX;
    data_mask       = (VIO_Volume)NULL;
    model_mask      = (VIO_Volume)NULL;
    mask_data_name  = (char *)NULL;
    mask_model_name = (char *)NULL;
                      /* look for optional objective function 
                   specification and optional weighting value */
    args_used   =2;
    if ( argc>2 ) {
                                /* objective functions */
      obj_func_index = weight_index = 2;
      if ( strncmp(argv[obj_func_index], "xcorr", 2)==0 ) {
        obj_func =  NONLIN_XCORR;
        weight_index++; args_used++;
      }
      if ( strncmp(argv[obj_func_index], "sqdiff", 2)==0 ) {
        obj_func =  NONLIN_SQDIFF;
        weight_index++; args_used++;
      }
      if ( strncmp(argv[obj_func_index], "diff", 2)==0 ) {
        obj_func =  NONLIN_DIFF;
        weight_index++; args_used++;
      }
      if ( strncmp(argv[obj_func_index], "label", 2)==0 ) {
        obj_func =  NONLIN_LABEL;
        weight_index++; args_used++;
      }
      if ( strncmp(argv[obj_func_index], "chamfer", 2)==0 ) {
        obj_func =  NONLIN_CHAMFER;
        weight_index++; args_used++;
      }
      if ( strncmp(argv[obj_func_index], "corrcoeff", 2)==0 ) {
        obj_func = NONLIN_CORRCOEFF;
        weight_index++; args_used++;
      }
      if ( strncmp(argv[obj_func_index], "opticalflow", 2)==0 ) {
        obj_func =  NONLIN_OPTICALFLOW;
        weight_index++; args_used++;
      }

                                /* weighting value */
      tmp = strtod(argv[weight_index],&end_ptr) ;
      if (end_ptr != argv[weight_index]) {
        weight = tmp;
        args_used++;
      }
        
    }

    if (obj_func == NONLIN_LABEL) { /* if the feature is a label, then load data as is */

      status = input_volume(data_name, 3, default_dim_names, 
			    NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			    TRUE, &data_vol, 
			    (minc_input_options *)NULL );
      if (status != VIO_OK) {
	(void)fprintf(stderr, "Cannot input feature %s.\n",data_name);
	return(-1);
      } 
      status = input_volume(model_name, 3, default_dim_names, 
			    NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			    TRUE, &model_vol, 
			    (minc_input_options *)NULL );
      if (status != VIO_OK) {
	(void)fprintf(stderr, "Cannot input feature %s.\n",model_name);
	return(-1);
      } 

    }
    else {			/* if feature is not a label, then force load as DOUBLEs */

      status = input_volume(data_name, 3, default_dim_names, 
			    NC_DOUBLE, FALSE, 0.0, 0.0,
			    TRUE, &data_vol, 
			    (minc_input_options *)NULL );
      if (status != VIO_OK) {
	(void)fprintf(stderr, "Cannot input feature %s.\n",data_name);
	return(-1);
      } 
      status = input_volume(model_name, 3, default_dim_names, 
			    NC_DOUBLE, FALSE, 0.0, 0.0,
			    TRUE, &model_vol, 
			    (minc_input_options *)NULL );
      if (status != VIO_OK) {
	(void)fprintf(stderr, "Cannot input feature %s.\n",model_name);
	return(-1);
      } 

    }



    add_a_feature_for_matching(&(main_args->features),
                               data_vol, model_vol, data_mask, model_mask,
                               data_name, model_name, 
                               mask_data_name, mask_model_name,
                               obj_func, weight,
                               thresh_data, thresh_model);


    i = main_args->features.number_of_features-1;
    print ("Features %d: %s %s %d %f\n", i,
           main_args->features.data_name[i],
           main_args->features.model_name[i],
           (int)main_args->features.obj_func[i],
           main_args->features.weight[i]);
    
  }
  else {
    fprintf (stderr,"the -feature option requires at least two arguments.\n");
    return (-1);
  }

  for(i=0; i<argc-args_used; i++) {
    argv[i] = argv[i+args_used];
  }
  argc -= args_used;

  return (argc);                        /* VIO_OK */
  
}

/* kate: indent-mode cstyle; indent-width 2; replace-tabs on; */
