/* ----------------------------- MNI Header -----------------------------------
@NAME       : measure_code.c
@DESCRIPTION: this file is included into minctracc.c, and contains code
              to init the sampling lattice and calculate the objective function
              value for each of the objective functions possible,
              writing the resulting values out to a file.
              
@MODIFIED   : $Log: measure_code.c,v $
@MODIFIED   : Revision 1.3  2006-11-30 09:07:31  rotor
@MODIFIED   :  * many more changes for clean minc 2.0 build
@MODIFIED   :
@MODIFIED   : Revision 1.1  1999/10/25 19:52:19  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :
 * Revision 1.4  1995/02/22  08:56:06  collins
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.3  94/04/06  11:48:42  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.2  93/11/15  16:27:05  louis
 * working version, with new library, with RCS revision stuff,
 * before deformations included
 * 
 * Revision 1.1  93/11/15  13:12:23  louis
 * Initial revision
 * 
---------------------------------------------------------------------------- */


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


    status = open_file(  main_args->filenames.measure_file, WRITE_FILE, BINARY_FORMAT,  &ofd );
    if ( status != VIO_OK ) 
      print_error_and_line_num ("filename `%s' cannot be opened.", 
                   __FILE__, __LINE__, main_args->filenames.measure_file);

                                /* do zscore and reload */

    main_args->obj_function = zscore_objective;
    obj_func_val = measure_fit( data, model, mask_data, mask_model, main_args );
    (void)fprintf (ofd, "%f - zscore\n",obj_func_val);
    (void)fflush(ofd);
    DEBUG_PRINT1 ( "%f - zscore\n",obj_func_val);
    
    delete_volume(data);  
    delete_volume(model);  


    status = input_volume( main_args->filenames.data, 3, default_dim_names, 
                          NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                          TRUE, &data, (minc_input_options *)NULL );
    status = input_volume( main_args->filenames.model, 3, default_dim_names, 
                          NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                          TRUE, &model, (minc_input_options *)NULL );

                                /* do xcorr   */

    main_args->obj_function = xcorr_objective;
    obj_func_val = measure_fit( data, model, mask_data, mask_model, main_args );
    (void)fprintf (ofd, "%f - xcorr\n",obj_func_val);
    (void)fflush(ofd);
    DEBUG_PRINT1 ( "%f - xcorr\n",obj_func_val);

                                /* do var_ratio */

    main_args->obj_function = vr_objective;
    obj_func_val = measure_fit( data, model, mask_data, mask_model, main_args );
    (void)fprintf (ofd, "%f - var_ratio\n",obj_func_val);
    (void)fflush(ofd);
    DEBUG_PRINT1 ( "%f - var_ratio\n",obj_func_val);

    delete_volume(data);  
    delete_volume(model);  

                                /* do ssc / zero-crossings */

    status = input_volume( main_args->filenames.data, 3, default_dim_names, 
                          NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                          TRUE, &data, (minc_input_options *)NULL );
    status = input_volume( main_args->filenames.model, 3, default_dim_names, 
                          NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                          TRUE, &model, (minc_input_options *)NULL ); 

    main_args->obj_function = ssc_objective;
    obj_func_val = measure_fit( data, model, mask_data, mask_model, main_args );
    (void)fprintf (ofd, "%f - ssc\n",obj_func_val);
    (void)fflush(ofd);
    DEBUG_PRINT1 ( "%f - ssc\n",obj_func_val);


    status = close_file(ofd);
    if ( status != VIO_OK ) 
      print_error_and_line_num ("filename `%s' cannot be closed.", 
                   __FILE__, __LINE__, main_args->filenames.measure_file);


    exit( status ) ;

