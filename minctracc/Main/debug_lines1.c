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

