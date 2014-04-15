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

   @CREATED    : February 3, 1992 - louis collins
   @MODIFIED   : $Log: minctracc.c,v $
   @MODIFIED   : Revision 96.22  2009/05/26 18:03:07  claude
   @MODIFIED   : free more memory after usage
   @MODIFIED   :
   @MODIFIED   : Revision 96.21  2009/05/22 15:49:19  claude
   @MODIFIED   : fixed memory bug freeing initial transform
   @MODIFIED   :
   @MODIFIED   : Revision 96.20  2009/04/03 18:36:59  louis
   @MODIFIED   : made changes to use only DOUBLES for input source and model volumes, and for all estimation of deformation fields
   @MODIFIED   :
   @MODIFIED   : Revision 96.19  2009/03/13 19:51:31  claude
   @MODIFIED   : fixed bug in offsets for minctracc and free memory upon exit
   @MODIFIED   :
   @MODIFIED   : Revision 96.18  2008/10/08 15:17:49  louis
   @MODIFIED   : added -nmi option for linear normalized mutual information
   @MODIFIED   :
   @MODIFIED   : Revision 96.17  2006/11/30 09:07:31  rotor
   @MODIFIED   :  * many more changes for clean minc 2.0 build
   @MODIFIED   :
   @MODIFIED   : Revision 96.15  2006/06/04 07:02:35  rotor
   @MODIFIED   :  * Fixed 64 bit function pointer and ParseArgv problem with an enum for
   @MODIFIED   :       objective function type and interpolation type (thanks jason)
   @MODIFIED   :
   @MODIFIED   : Revision 96.14  2005/07/20 20:45:48  rotor
   @MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
   @MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
   @MODIFIED   :     * Removed old VOLUME_IO cruft #defines
   @MODIFIED   :     * Fixed up all Makefile.am's in subdirs
   @MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
   @MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
   @MODIFIED   :
   @MODIFIED   : Revision 96.13  2005/06/28 18:56:14  rotor
   @MODIFIED   :  * added masking for feature volumes (irina and patricia)
   @MODIFIED   :
   @MODIFIED   : Revision 96.12  2004/03/18 06:51:03  rotor
   @MODIFIED   :  * changed make_model from csh to sh
   @MODIFIED   :  * removed an extraneous printf from minctracc
   @MODIFIED   :
   @MODIFIED   : Revision 96.11  2004/02/12 05:54:22  rotor
   @MODIFIED   :  * removed public/private defs
   @MODIFIED   :
   @MODIFIED   : Revision 96.10  2004/02/04 20:42:33  lenezet
   @MODIFIED   : *** empty log message ***
   @MODIFIED   :
   @MODIFIED   : Revision 96.9  2003/02/04 06:08:44  stever
   @MODIFIED   : Add support for correlation coefficient and sum-of-squared difference.
   @MODIFIED   :
   @MODIFIED   : Revision 96.8  2002/12/13 21:16:11  lenezet
   @MODIFIED   : nonlinear in 2D has changed. The option -2D-non-lin is no more necessary. The grid transform has been adapted to feet on the target volume whatever is size. The Optimization is done on the dimensions for which "count" is greater than 1.
   @MODIFIED   :
   @MODIFIED   : Revision 96.7  2002/11/20 21:38:31  lenezet
   @MODIFIED   :
   @MODIFIED   : Fix the code to take in consideration the direction cosines especially in the grid transform.
   @MODIFIED   : Add an option to choose the maximum expected deformation magnitude.
   @MODIFIED   :
   @MODIFIED   : Revision 96.6  2002/08/14 19:54:42  lenezet
   @MODIFIED   :  quaternion option added for the rotation
   @MODIFIED   :
   @MODIFIED   : Revision 96.5  2002/03/26 14:15:38  stever
   @MODIFIED   : Update includes to <volume_io/foo.h> style.
   @MODIFIED   :
   @MODIFIED   : Revision 96.4  2002/03/07 19:08:34  louis
   @MODIFIED   : Added -lattice_diameter as an optionto minctracc to account for a
   @MODIFIED   : problem with the automated calculation of the sub-lattice diameter.
   @MODIFIED   : It used to be step*3*2 - which was pretty big, when step = 8mm.
   @MODIFIED   :
   @MODIFIED   : Now, the sub lattice diameter can be input on the command line, and I
   @MODIFIED   : suggest a lattice size 3 times greater than the step size.
   @MODIFIED   :
   @MODIFIED   : If not on the command line, the default is = 24mm.
   @MODIFIED   :
   @MODIFIED   : Revision 96.3  2000/03/15 08:42:41  stever
   @MODIFIED   : Code cleanup: all functions prototyped (except ParseArgs.c), no useless declarations, etc
   @MODIFIED   :
   @MODIFIED   : Revision 96.2  2000/02/20 04:01:03  stever
   @MODIFIED   : * use new history_string() function to generate history strings
   @MODIFIED   :   when outputting MNI files (.mnc, .xfm)
   @MODIFIED   : * removed unused vax routines from Proglib
   @MODIFIED   : * tuned configure script; CPPFLAGS and LDFLAGS are now left alone,
   @MODIFIED   :   for the installer to use
   @MODIFIED   :
   @MODIFIED   : Revision 96.1  1999/10/25 19:52:19  louis
   @MODIFIED   : final checkin before switch to CVS
   @MODIFIED   :
 * Revision 96.0  1996/08/21  18:21:51  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:21:49  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:15:42  louis
 * Never released version 0.95
 *
 * Revision 1.16  1996/08/12  14:15:40  louis
 * Pre-release
 *
 * Revision 1.15  1995/09/28  13:24:26  collins
 * added multiple feature volume loading with get_feature_volumes()
 * to be able to correlate multiple features at the same time.
 *
 * Revision 1.14  1995/09/28  11:54:52  collins
 * working version, just prior to release 0.9 of mni_autoreg
 *
 * Revision 1.13  1995/02/22  08:56:06  collins
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
static char minctracc_rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Main/minctracc.c,v 96.22 2009/05/26 18:03:07 claude Exp $";
#endif

#include <config.h>
#include <float.h>
#include <ParseArgv.h>
#include <volume_io.h>
#include <minctracc.h>
#include <globals.h>
#include <objectives.h>
#include "local_macros.h"

/* objective function for nonlinear optimization.
 * Set in get_nonlinear_objective().
 */
static int obj_func0 = -1;

VIO_Real initial_corr, final_corr;
static char *default_dim_names[VIO_N_DIMENSIONS] = 
    { MIzspace, MIyspace, MIxspace };

/* Why are these declared here?  They aren't apparently used.
 * -smr
 *
 * static   const VIO_STR      TRANSFORM_FILE_HEADER = "MNI VIO_Transform File";
 * static   const VIO_STR      LINEAR_TYPE = "Linear";
 * static   const VIO_STR      TYPE_STRING = "Transform_Type";
 * static   const VIO_STR      LINEAR_TRANSFORM_STRING = "Linear_Transform";
 * static   const VIO_STR      GRID_TRANSFORM_STRING = "Grid_Transform";
 * static   const VIO_STR      DISPLACEMENT_VOLUME = "Displacement_Volume";
 */

/*************************************************************************/
int main ( int argc, char* argv[] )
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
    (strlen(main_args.filenames.matlab_file)  != 0) ||
    (strlen(main_args.filenames.measure_file) != 0);

  /* assign objective function and interpolant type */
  switch (main_args.interpolant_type) {
  case TRICUBIC:
    main_args.interpolant = tricubic_interpolant;
    break;
  case TRILINEAR:
    main_args.interpolant = trilinear_interpolant;
    break;
  case N_NEIGHBOUR:
    main_args.interpolant = nearest_neighbour_interpolant;
    break;
  default:
    (void) fprintf(stderr, "Error determining interpolation type\n");
    exit(EXIT_FAILURE);
  }

  switch (main_args.obj_function_type) {
  case XCORR:
    main_args.obj_function = xcorr_objective;
    break;
  case ZSCORE:
    main_args.obj_function = zscore_objective;
    break;
  case SSC:
    main_args.obj_function = ssc_objective;
    break;
  case VR:
    main_args.obj_function = vr_objective;
    break;
  case MUTUAL_INFORMATION:
    main_args.obj_function = mutual_information_objective;
    break;
  case NORMALIZED_MUTUAL_INFORMATION:
    main_args.obj_function = normalized_mutual_information_objective;
    break;
  default:
    (void) fprintf(stderr, "Error determining objective function type\n");
    exit(EXIT_FAILURE);
  }

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
        ((strlen(main_args.filenames.matlab_file)  != 0) ||
         (strlen(main_args.filenames.measure_file) != 0))) {
      (void)fprintf(stderr, "\nNote: No output transform file needs to be specified for -matlab\n");
      (void)fprintf(stderr, "      or -measure options.\n");
    }

    exit(EXIT_FAILURE);
  }

  main_args.filenames.data  = argv[1];        /* set up necessary file names */
  main_args.filenames.model = argv[2];
  if (strlen(main_args.filenames.measure_file)==0 &&
      strlen(main_args.filenames.matlab_file)==0) 
    main_args.filenames.output_trans = argv[3];


                                /* check to see if they can be overwritten */
  if (!clobber_flag && 
      (strlen(main_args.filenames.measure_file)!=0) && 
      file_exists(main_args.filenames.measure_file)) {
    (void)fprintf (stderr,"Measure file %s exists.\n",main_args.filenames.measure_file);
    (void)fprintf (stderr,"Use -clobber to overwrite.\n");
    exit(EXIT_FAILURE);
  }

  if (!clobber_flag && 
      (strlen(main_args.filenames.matlab_file)!=0) && 
      file_exists(main_args.filenames.matlab_file)) {
    (void)fprintf (stderr,"Matlab file %s exists.\n",main_args.filenames.matlab_file);
    (void)fprintf (stderr,"Use -clobber to overwrite.\n");
    exit(EXIT_FAILURE);
  }

  if (!clobber_flag && 
      (strlen(main_args.filenames.output_trans)!=0) && 
      file_exists(main_args.filenames.output_trans)) {
    (void)fprintf (stderr,"Output file %s exists.\n",main_args.filenames.output_trans);
    (void)fprintf (stderr,"Use -clobber to overwrite.\n");
    exit(EXIT_FAILURE);
  }
  if (strlen(main_args.filenames.matlab_file)  != 0 &&
      strlen(main_args.filenames.measure_file) != 0) {
    (void)fprintf(stderr, "\nWARNING: -matlab and -measure are mutually exclusive.  Only\n");
    (void)fprintf(stderr, "          the matlab option will be executed here.\n");
  }

                                /* set up linear transformation to identity, if
                                   not set up in the argument list */

  if (main_args.trans_info.transformation == (VIO_General_transform *)NULL) { 
                                /* build identity VIO_General_transform */
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
    if (strlen(main_args.filenames.output_trans) != 0)
      print ( "Output filename     = %s\n", main_args.filenames.output_trans);
    if (strlen(main_args.filenames.matlab_file)  != 0)
      print ( "Matlab filename     = %s (num_steps=%d)\n\n",main_args.filenames.matlab_file,Matlab_num_steps );
    if (strlen(main_args.filenames.measure_file) != 0)
      print ( "Measure filename    = %s\n\n",main_args.filenames.measure_file );
    print ( "Step size           = %f %f %f\n",
                  main_args.step[0],
                  main_args.step[1],
                  main_args.step[2]);
    print ( "Sub-lattice dia     = %f %f %f\n",
                  main_args.lattice_width[0],
                  main_args.lattice_width[1],
                  main_args.lattice_width[2]);
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
            if (main_args.obj_function == mutual_information_objective || main_args.obj_function == normalized_mutual_information_objective ) {
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
      for(i=0; i<4; i++) print ("%9.4f ",Transform_elem(*lt,0,i));
      print ( "\n" );
      print ( "                      ");
      for(i=0; i<4; i++) print ("%9.4f ",Transform_elem(*lt,1,i));
      print ( "\n" );
      print ( "                      ");
      for(i=0; i<4; i++) print ("%9.4f ",Transform_elem(*lt,2,i));
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
    else 
      {
        print ( "Transform center   = %8.3f %8.3f %8.3f\n", 0.0,0.0,0.0);
      }

    /* rotation with quaternion or euler angle */
    if (main_args.trans_info.rotation_type == TRANS_QUAT)
      {
        quat4=sqrt(1-SQR(main_args.trans_info.quaternions[0])-SQR(main_args.trans_info.quaternions[1])-SQR(main_args.trans_info.quaternions[2]));
        print ( "Transform quaternion   = %8.3f %8.3f %8.3f %8.3f \n\n", 
                main_args.trans_info.quaternions[0],
                main_args.trans_info.quaternions[1],
                main_args.trans_info.quaternions[2],
                quat4 );
      }
    if (main_args.trans_info.rotation_type == TRANS_ROT)
      {
        print ( "Transform rotation   = %8.3f %8.3f %8.3f \n\n", 
                main_args.trans_info.rotations[0],
                main_args.trans_info.rotations[1],
                main_args.trans_info.rotations[2] );
      }
    print ( "Transform trans    = %8.3f %8.3f %8.3f\n", 
            main_args.trans_info.translations[0],
            main_args.trans_info.translations[1],
            main_args.trans_info.translations[2] );
    print ( "Transform scale    = %8.3f %8.3f %8.3f\n\n", 
            main_args.trans_info.scales[0],
            main_args.trans_info.scales[1],
            main_args.trans_info.scales[2] );
    

  }


  


  if (main_args.trans_info.use_magnitude) 
    {
      /* non-linear optimization is based on the correlation of a local
         sub-lattice between source and target gradient magnitude volumes 
         for the first two volumes */

      if (main_args.trans_info.transform_type == TRANS_NONLIN)
        DEBUG_PRINT1( "This run will use sub-lattice correlation (type %d) between the two input vols.\n", obj_func0);
    }
  else 
    {
      /* non-linear optimization is based on the direct computation of a
         deformation vector, based on optical flow.  Both volume _MUST_
         be intensity normalized! */

      if (main_args.trans_info.transform_type == TRANS_NONLIN)
        DEBUG_PRINT( "This run will use optical flow.\n");
    }

  ALLOC(data,1);

  status = input_volume( main_args.filenames.data, 3, default_dim_names, 
                         NC_DOUBLE, FALSE, 0.0, 0.0,
                         TRUE, &data, (minc_input_options *)NULL );

  if (status != OK)
    print_error_and_line_num("Cannot input volume '%s'",
                             __FILE__, __LINE__,main_args.filenames.data);
  data_dxyz = data;
 
  status = input_volume( main_args.filenames.model, 3, default_dim_names, 
                         NC_DOUBLE, FALSE, 0.0, 0.0,
                         TRUE, &model, (minc_input_options *)NULL );
  if (status != OK)
    print_error_and_line_num("Cannot input volume '%s'",
                             __FILE__, __LINE__,main_args.filenames.model);
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
                                __FILE__, __LINE__, main_args.filenames.data, 
                                get_volume_n_dimensions(data));
    }
  if (get_volume_n_dimensions(model)!=3) 
    {
      print_error_and_line_num ("Model file %s has %d dimensions.  Only 3 dims supported.", 
                                __FILE__, __LINE__, main_args.filenames.model, 
                                get_volume_n_dimensions(model));
    }


                                /* shift features to be able to
                                   insert the main source/target
                                   volumes first. */

  num_features = allocate_a_new_feature(&(main_args.features));
  for(i=num_features; i>=1; i--) {
    main_args.features.data[i]            = main_args.features.data[i-1];
    main_args.features.model[i]           = main_args.features.model[i-1];
    main_args.features.data_name[i]       = main_args.features.data_name[i-1];
    main_args.features.model_name[i]      = main_args.features.model_name[i-1];
    main_args.features.data_mask[i]       = main_args.features.data_mask[i-1];
    main_args.features.model_mask[i]      = main_args.features.model_mask[i-1];
    main_args.features.mask_data_name[i]  = main_args.features.mask_data_name[i-1];
    main_args.features.mask_model_name[i] = main_args.features.mask_model_name[i-1];
    main_args.features.thresh_data[i]     = main_args.features.thresh_data[i-1];
    main_args.features.thresh_model[i]    = main_args.features.thresh_model[i-1];
    main_args.features.obj_func[i]        = main_args.features.obj_func[i-1];
    main_args.features.weight[i]          = main_args.features.weight[i-1];
  }

  main_args.features.data[0]            = data; 
  main_args.features.model[0]           = model;
  main_args.features.data_name[0]       = main_args.filenames.data;
  main_args.features.model_name[0]      = main_args.filenames.model;
  main_args.features.data_mask[0]       = mask_data;
  main_args.features.model_mask[0]      = mask_model;
  main_args.features.mask_data_name[0]  = main_args.filenames.mask_data;
  main_args.features.mask_model_name[0] = main_args.filenames.mask_model;
  main_args.features.thresh_data[0]     = main_args.threshold[0];
  main_args.features.thresh_model[0]    = main_args.threshold[1];
  if (main_args.trans_info.use_magnitude) {
    main_args.features.obj_func[0]        = obj_func0;
  } 
  else {
    main_args.features.obj_func[0]        = NONLIN_OPTICALFLOW;    
  }
  main_args.features.weight[0]          = 1.0;

  /* ===========================  translate initial transformation matrix into 
                                  transformation parameters */

    if (!init_params( data, model, mask_data, mask_model, &main_args )) {
      print_error_and_line_num("%s",__FILE__, __LINE__,
                             "Could not initialize transformation parameters\n");
    }



  if (strlen(main_args.filenames.matlab_file) != 0) 
    {
      init_lattice( data, model, mask_data, mask_model, &main_args );
      make_matlab_data_file( data, model, mask_data, mask_model, comments, &main_args );
      exit( OK );
    }

  if (strlen(main_args.filenames.measure_file) != 0) 
    {

#include "measure_code.c"

    /* measure code finishes with 
       exit(status); */
    }

  DEBUG_PRINT  ("AFTER init_params()\n");
  if (get_transform_type(main_args.trans_info.transformation)==LINEAR) {
    lt = get_linear_transform_ptr(main_args.trans_info.transformation);
    
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
                main_args.trans_info.center[0],
                main_args.trans_info.center[1],
                main_args.trans_info.center[2] );
  if (main_args.trans_info.rotation_type == TRANS_ROT){
    DEBUG_PRINT3 ( "Transform rotations  = %8.3f %8.3f %8.3f\n", 
                   main_args.trans_info.rotations[0],
                   main_args.trans_info.rotations[1],
                   main_args.trans_info.rotations[2] );
  }
  if (main_args.trans_info.rotation_type == TRANS_QUAT)
    {
      quat4=sqrt(1-SQR(main_args.trans_info.quaternions[0])-SQR(main_args.trans_info.quaternions[1])-SQR(main_args.trans_info.quaternions[2]));
      DEBUG_PRINT4 ( "Transform quaternions  = %8.3f %8.3f %8.3f %8.3f\n", 
                     main_args.trans_info.quaternions[0],
                     main_args.trans_info.quaternions[1],
                     main_args.trans_info.quaternions[2],
                     quat4 );
    }
  

  DEBUG_PRINT3 ( "Transform trans    = %8.3f %8.3f %8.3f\n", 
                main_args.trans_info.translations[0],
                main_args.trans_info.translations[1],
                main_args.trans_info.translations[2] );
  DEBUG_PRINT3 ( "Transform scale    = %8.3f %8.3f %8.3f\n", 
                main_args.trans_info.scales[0],
                main_args.trans_info.scales[1],
                main_args.trans_info.scales[2] );
  DEBUG_PRINT3 ( "Transform shear    = %8.3f %8.3f %8.3f\n\n", 
                main_args.trans_info.shears[0],
                main_args.trans_info.shears[1],
                main_args.trans_info.shears[2] );
  


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
      

      if ( !optimize_non_linear_transformation( &main_args ) ) {
        print_error_and_line_num("Error in optimization of non-linear transformation\n",
                                 __FILE__, __LINE__);
        exit(EXIT_FAILURE);
      }
      
    }
    else {
      
      if (main_args.trans_info.rotation_type == TRANS_ROT )
        {
          if (!optimize_linear_transformation( data, model, mask_data, mask_model, &main_args )) {
            print_error_and_line_num("Error in optimization of linear transformation\n",
                                 __FILE__, __LINE__);
        exit(EXIT_FAILURE);
          }
        }
      
      
      if (main_args.trans_info.rotation_type == TRANS_QUAT )
        {
          if (!optimize_linear_transformation_quater( data, model, mask_data, mask_model, &main_args )) {
            print_error_and_line_num("Error in optimization of linear transformation\n",
                                 __FILE__, __LINE__);
        exit(EXIT_FAILURE);
          }
        }
      
      
    }

    if (number_dimensions==3) {
      print ("Initial objective function val = %0.8f\n",initial_corr); 
      print ("Final objective function value = %0.8f\n",final_corr);
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
  if( comments ) {
    FREE( comments );
    comments = NULL;
  }
  
  delete_general_transform(main_args.trans_info.transformation);
  FREE(main_args.trans_info.transformation);
  delete_general_transform(main_args.trans_info.orig_transformation);
  FREE(main_args.trans_info.orig_transformation);
  free_features (&(main_args.features));
  FREE(main_args.trans_info.file_contents);
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

   delete_general_transform(&input_transformation);

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
int get_mask_file(char *dst, char *key, char *nextArg)
{ 

  VIO_Status status;

  if (strncmp ( "-model_mask", key, 2) == 0) {
    /*    ALLOC( mask_model, 1 );*/
    status = input_volume( nextArg, 3, default_dim_names, 
                          NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                          TRUE, &mask_model, (minc_input_options *)NULL );
    dst = nextArg;
    main_args.filenames.mask_model = nextArg;
  }
  else {
    /*    ALLOC( mask_data, 1);*/
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
    added: set of main_args.filenames.mask_model and .mask_data

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
      if (status != OK) {
	(void)fprintf(stderr, "Cannot input feature %s.\n",data_name);
	return(-1);
      } 
      status = input_volume(model_name, 3, default_dim_names, 
			    NC_UNSPECIFIED, FALSE, 0.0, 0.0,
			    TRUE, &model_vol, 
			    (minc_input_options *)NULL );
      if (status != OK) {
	(void)fprintf(stderr, "Cannot input feature %s.\n",model_name);
	return(-1);
      } 

    }
    else {			/* if feature is not a label, then force load as DOUBLEs */

      status = input_volume(data_name, 3, default_dim_names, 
			    NC_DOUBLE, FALSE, 0.0, 0.0,
			    TRUE, &data_vol, 
			    (minc_input_options *)NULL );
      if (status != OK) {
	(void)fprintf(stderr, "Cannot input feature %s.\n",data_name);
	return(-1);
      } 
      status = input_volume(model_name, 3, default_dim_names, 
			    NC_DOUBLE, FALSE, 0.0, 0.0,
			    TRUE, &model_vol, 
			    (minc_input_options *)NULL );
      if (status != OK) {
	(void)fprintf(stderr, "Cannot input feature %s.\n",model_name);
	return(-1);
      } 

    }



    add_a_feature_for_matching(&(main_args.features),
                               data_vol, model_vol, data_mask, model_mask,
                               data_name, model_name, 
                               mask_data_name, mask_model_name,
                               obj_func, weight,
                               thresh_data, thresh_model);


    i = main_args.features.number_of_features-1;
    print ("Features %d: %s %s %d %f\n", i,
           main_args.features.data_name[i],
           main_args.features.model_name[i],
           (int)main_args.features.obj_func[i],
           main_args.features.weight[i]);
    
  }
  else {
    fprintf (stderr,"the -feature option requires at least two arguments.\n");
    return (-1);
  }

  for(i=0; i<argc-args_used; i++) {
    argv[i] = argv[i+args_used];
  }
  argc -= args_used;

  return (argc);                        /* OK */
  
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
    main_args.trans_info.transform_type = TRANS_NONLIN;

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

  i = main_args.features.number_of_features;
  main_args.features.number_of_features++;    

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
