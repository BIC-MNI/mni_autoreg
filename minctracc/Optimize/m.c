
/**************************************************************************/

public Status do_non_linear_optimization(Volume d1,
					 Volume d1_dx, 
					 Volume d1_dy, 
					 Volume d1_dz, 
					 Volume d1_dxyz,
					 Volume d2,
					 Volume d2_dx, 
					 Volume d2_dy, 
					 Volume d2_dz, 
					 Volume d2_dxyz,
					 Volume m1,
					 Volume m2, 
					 Arg_Data *globals)
{
  General_transform
    *all_until_last, 
    *additional_warp,
    *current_warp;
  
  Volume
    additional_vol,
    current_vol,
    additional_mag;

  Volume
    additional_dx,additional_dy,additional_dz;

  Deform_field 
    *current;
  long
    timer1,timer2,
    nfunk_total;
  int 
    loop_start[MAX_DIMENSIONS], 
    loop_end[MAX_DIMENSIONS],
    iters,
    i,j,k,m,n,
    nodes_done, nodes_tried, nodes_seen, over,
    nfunks,
    nfunk1, nodes1,
    sizes[MAX_DIMENSIONS];
  Real 
    steps[MAX_DIMENSIONS],
    steps_data[MAX_DIMENSIONS],
    min, max, sum, sum2, mag, mean_disp_mag, std, var, 
    voxel,val_x, val_y, val_z,
    def_x, def_y, def_z, 
    wx,wy,wz,
    mx,my,mz,
    tx,ty,tz, 
    displace,
    zero,
    threshold1,
    threshold2,
    result;
  progress_struct
    progress;

	       	/* set up globals for communication with other routines */
  Gd1     = d1;
  Gd1_dx  = d1_dx; 
  Gd1_dy  = d1_dy; 
  Gd1_dz  = d1_dz; 
  Gd1_dxyz= d1_dxyz;
  Gd2     = d2;
  Gd2_dx  = d2_dx; 
  Gd2_dy  = d2_dy; 
  Gd2_dz  = d2_dz; 
  Gd2_dxyz= d2_dxyz;
  Gm1     = m1;
  Gm2     = m2; 
  Gglobals= globals;


  Ga1xyz = vector(1,512);	/* allocate space for the global data for */
  Ga2xyz = vector(1,512);	/* the local neighborhood lattice values */
  
  SX = vector(1,512);		/* and coordinates in source volume  */
  SY = vector(1,512);
  SZ = vector(1,512);
  TX = vector(1,512);		/* and coordinates in target volume  */
  TY = vector(1,512);
  TZ = vector(1,512);

				/* extract the deformation field to be optimized
				   from the current transformation     */
  split_up_the_transformation(globals->trans_info.transformation,
			      &all_until_last,
			      &current_warp);
  
  if (current_warp == (General_transform *)NULL) { /* exit if no deformation */
    print_error_and_line_num("Cannot find the deformation field transform to optimize",
		__FILE__, __LINE__);
  }

  if (globals->flags.debug) {	/* print some debugging info    */
    print("orig transform is %d long\n",
	  get_n_concated_transforms(globals->trans_info.transformation));
    print("all_until_last is %d long\n",
	  get_n_concated_transforms(all_until_last));
  }

  copy_general_transform(current_warp,
			 additional_warp);

  current_vol =    current_warp->displacement_volume;
  additional_vol = additional_warp->displacement_volume;

  get_volume_sizes(additional_vol, sizes);
  get_volume_separations(additional_vol, steps);

  get_volume_separations(d2_dxyz, steps_data);

  if (steps_data[0]!=0.0) {
    Gsimplex_size= ABS(steps[0]/steps_data[0]);
    if (ABS(Gsimplex_size) < ABS(steps_data[0])) {
      print ("*** WARNING ***\n");
      print ("Simplex size will be smaller than data voxel size (%f < %f)\n",
	     Gsimplex_size,steps_data[0]);
    }
  }
  else
    print_error_and_line_num("Zero step size for gradient data2: %f %f %f\n", 
		__FILE__, __LINE__,steps_data[0],steps_data[1],steps_data[2]);


  Gcost_radius = 8*Gsimplex_size*Gsimplex_size*Gsimplex_size;
  Glinear_transform = all_until_last; /* set up linear part of transformation   */
  additional_mag = copy_volume_definition(current->dz, NC_UNSPECIFIED, FALSE, 0.0,0.0);




  /*******************************************************************************/
  /*   initialize this iterations warp                                           */

  zero = CONVERT_VALUE_TO_VOXEL(additional_vol, 0.0);	
  for_less(i,0,sizes[0])
    for_less(j,0,sizes[1])
      for_less(k,0,sizes[2])
	for_less(m,0,sizes[3])
	  for_less(n,0,sizes[4]){
	    SET_VOXEL(additional_vol, i,j,k,m,n, zero);
	  }


  zero = CONVERT_VALUE_TO_VOXEL(additional_mag, 0.0);	
  for_less(i,0,sizes[0])
    for_less(j,0,sizes[1])
      for_less(k,0,sizes[2]) {
	SET_VOXEL(additional_mag, i,j,k,0,0, zero);
      }
  mean_disp_mag = 0.0;


  /*******************************************************************************/
  /*    set the threshold to be 10% of the maximum gradient magnitude            */
  /*    for each source and target volumes                                       */

  if (globals->threshold[0]==0.0)
    threshold1 = 0.10 * get_maximum_magnitude(d1_dxyz);
  else
    threshold1 = globals->threshold[0] * get_maximum_magnitude(d1_dxyz);

  if (globals->threshold[1]==0.0)
    threshold2 = 0.10 * get_maximum_magnitude(d2_dxyz);
  else
    threshold2 = globals->threshold[1] * get_maximum_magnitude(d2_dxyz);

  if (threshold1<0.0 || threshold2<0.0) {
    print_error_and_line_num("Gradient magnitude threshold error: %f %f\n", __FILE__, __LINE__, 
		threshold1, threshold2);
  }
  if (globals->flags.debug) {	
    print("Source vol threshold = %f\n", threshold1);
    print("Target vol threshold = %f\n", threshold2);
    print("Iteration limit      = %d\n", iteration_limit);
    print("Iteration weight     = %f\n", iteration_weight);
  }

  /*******************************************************************************/
  /*    start the iterations to estimate the warp                     

	for each iteration {
	   for each voxel in the deformation field {
	      get the original coordinate node
	      find the best deformation for the node, taking into consideration
	        the previous deformation
	      store this best deformation
	   }

	   for each voxel in the deformation field {
 	      add WEIGHTING*best_deformation to the current deformation
	   }

	}
  */

  for_less(i,0,MAX_DIMENSIONS) {
    if (sizes[i]>3) {
      loop_start[i] = 1;
      loop_end[i] = sizes[i]-1;
    }
    else {
      loop_start[i]=0;
      loop_end[i] = sizes[i];
    }
  }

  if (globals->flags.debug) {
    print("loop: (%d %d) (%d %d) (%d %d)\n",
	  loop_start[0],loop_end[0],loop_start[1],loop_end[1],loop_start[2],loop_end[2]);
  }

  if (globals->trans_info.use_super) {
    make_super_sampled_data(current->dx, &Gsuper_dx,  &Gsuper_dy,  &Gsuper_dz);
  }
  else {
    Gsuper_dx = current->dx;
    Gsuper_dy = current->dy;
    Gsuper_dz = current->dz;
  }


  for_less(iters,0,iteration_limit) {

    print("Iteration %2d of %2d\n",iters+1, iteration_limit);

    if (globals->trans_info.use_super)
      interpolate_super_sampled_data(current->dx, current->dy, current->dz,
				     Gsuper_dx,  Gsuper_dy,  Gsuper_dz,
				     number_dimensions);
  
    nodes_done = 0; nodes_tried = 0; nodes_seen=0; displace = 0.0; over = 0;
    nfunk_total = 0;

    sum = sum2 = 0.0;
    min = 1000.0;
    max = -1000.0;

    initialize_progress_report( &progress, FALSE, 
			       (loop_end[0]-loop_start[0])*(loop_end[1]-loop_start[1]) + 1,
			       "Estimating deformations" );

    for_less(i,loop_start[0],loop_end[0]) {

      timer1 = time(NULL);
      nfunk1 = 0; nodes1 = 0;

      for_less(j,loop_start[1],loop_end[1]) {
	for_less(k,loop_start[2],loop_end[2]){
	  
	  nodes_seen++;
	  def_x = def_y = def_z = 0.0;
	  
				/* get the lattice coordinate */
	  convert_3D_voxel_to_world(current->dx, 
				    (Real)i, (Real)j, (Real)k,
				    &wx, &wy, &wz);

				/* get the warp to be added to the target point */
	  GET_VOXEL_3D(voxel, current->dx, i,j,k); 
	  val_x = CONVERT_VOXEL_TO_VALUE(current->dx , voxel ); 
	  GET_VOXEL_3D(voxel, current->dy, i,j,k); 
	  val_y = CONVERT_VOXEL_TO_VALUE(current->dy , voxel ); 
	  GET_VOXEL_3D(voxel, current->dz, i,j,k); 
	  val_z = CONVERT_VOXEL_TO_VALUE(current->dz , voxel ); 

				/* add the warp to the target lattice point */
	  wx += val_x; wy += val_y; wz += val_z;

	  if (point_not_masked(m2, wx, wy, wz) &&
	      get_value_of_point_in_volume(wx,wy,wz, Gd2_dxyz)>threshold2) {


			    /* now get the mean warped position of the target's neighbours */
/*
!!!
	    if (!get_average_warp_of_neighbours(current_warp,
					   i,j,k,
					   &mx, &my, &mz))
	      {
	      }
*/
    
	    general_inverse_transform_point(globals->trans_info.transformation,
					    wx,wy,wz,
					    &tx,&ty,&tz);


	    result = optimize_3D_deformation_for_single_node(steps[0], 
							     threshold1,
							     tx,ty,tz,
							     mx, my, mz,
							     &def_x, &def_y, &def_z, 
							     iters, iteration_limit,
							     &nfunks,
							     number_dimensions);
	    

	    if (result == -40.0) {
	      nodes_tried++;
	      result = 0.0;
	    } else {
	      if (ABS(result) > 0.95*steps[0])
		over++;

	      nodes_done++;

	      displace += ABS(result);
	      sum += ABS(result);
	      sum2 += ABS(result) * ABS(result);

	      if (ABS(result)>max) max = ABS(result);
	      if (ABS(result)<min) min = ABS(result);

	      nfunk_total += nfunks;

	      nfunk1 += nfunks; nodes1++;

	    }
	    
	  }
	  
	  mag = sqrt(def_x*def_x + def_y*def_y + def_z*def_z);
	  mag = CONVERT_VALUE_TO_VOXEL(additional_dx, mag); 
	  SET_VOXEL_3D(additional_mag, i,j,k, mag);
	  
	  def_x = CONVERT_VALUE_TO_VOXEL(additional_dx, def_x); 
	  SET_VOXEL_3D(additional_dx, i,j,k, def_x);
	  def_y = CONVERT_VALUE_TO_VOXEL(additional_dy, def_y); 
	  SET_VOXEL_3D(additional_dy, i,j,k, def_y);
	  def_z = CONVERT_VALUE_TO_VOXEL(additional_dz, def_z); 
	  SET_VOXEL_3D(additional_dz, i,j,k, def_z);
	  
	}

	update_progress_report( &progress, 
			       (loop_end[1]-loop_start[1])*(i-loop_start[0])+(j-loop_start[1])+1 );
      }
      timer2 = time(NULL);

      if (globals->flags.debug) 
	print ("slice: (%d : %d) = %d sec -- nodes=%d av funks %f\n",
	       i+1-loop_start[0], loop_end[0]-loop_start[0]+1, timer2-timer1, 
	       nodes1,
	       nodes1==0? 0.0:(float)nfunk1/(float)nodes1);
      
    }
    terminate_progress_report( &progress );

    if (globals->flags.debug) {
      if (nodes_done>0) {
	mean_disp_mag = displace/nodes_done;
	var = ((sum2 * nodes_done) - sum*sum) / ((float)nodes_done*(float)(nodes_done-1));
	std = sqrt(var);
	nfunks = nfunk_total / nodes_done;
      }
      else {
	mean_disp_mag=0.0; std = 0.0;
      }
      print ("Nodes seen = %d, tried = %d, done = %d, avg disp = %f +/- %f\n",
	     nodes_seen, nodes_tried, nodes_done, mean_disp_mag, std);
      print ("av nfunks = %d , over = %d, max disp = %f, min disp = %f\n", nfunks, over, max, min);

      nodes_tried = 0;
      for_less(i,0,sizes[0])
	for_less(j,0,sizes[1])
	  for_less(k,0,sizes[2]){
	    GET_VOXEL_3D(voxel, additional_mag, i,j,k); 
	    mag = CONVERT_VOXEL_TO_VALUE(additional_mag , voxel ); 
	    if (mag >= (mean_disp_mag+std))
	      nodes_tried++;
	  }
      print ("there are %d of %d over (mean+1std) disp.\n", nodes_tried, nodes_done);
  

    }


				/* update the current warp, so that the
				   next iteration will use all the data
				   calculated thus far.
				   (the result goes into additional_d* )   */

    add_additional_warp_to_current(additional_warp,
				   current_warp,
				   iteration_weight);

				/* smooth the warp (result into current->d*)  */

    smooth_the_warp(current_warp,
		    additional_warp,
		    additional_mag, 0.0);

    
/*                                 clamp the data so that the 1st derivative of
				   the deformation field does not exceed 1.0*step
				   in magnitude 

     clamp_warp_deriv(current->dx, current->dy, current->dz); 
*/

    
    if (iters<iteration_limit-1 && Gglobals->trans_info.use_simplex==TRUE) {
      for_less(i,loop_start[0],loop_end[0]) {
	for_less(j,loop_start[1],loop_end[1]) {
	  for_less(k,loop_start[2],loop_end[2]){
	    
	    def_x = def_y = def_z = 0.0;
	    
	    GET_VOXEL_3D(voxel, additional_mag, i,j,k); 
	    mag = CONVERT_VOXEL_TO_VALUE(additional_mag , voxel ); 
	    if (mag >= (mean_disp_mag+std)) {
	      
	      /* get the lattice coordinate */
	      convert_3D_voxel_to_world(current->dx, 
					(Real)i, (Real)j, (Real)k,
					&wx, &wy, &wz);
	      
	      /* get the warp to be added to the target point */
	      GET_VOXEL_3D(voxel, current->dx, i,j,k); 
	      val_x = CONVERT_VOXEL_TO_VALUE(current->dx , voxel ); 
	      GET_VOXEL_3D(voxel, current->dy, i,j,k); 
	      val_y = CONVERT_VOXEL_TO_VALUE(current->dy , voxel ); 
	      GET_VOXEL_3D(voxel, current->dz, i,j,k); 
	      val_z = CONVERT_VOXEL_TO_VALUE(current->dz , voxel ); 
	      
	      /* add the warp to the target lattice point */
	      wx += val_x; wy += val_y; wz += val_z;
	      
	      if ( point_not_masked(m2, wx, wy, wz)) {
		
		/* now get the mean warped position of the target's neighbours */
		
/*
!!!
		get_average_warp_of_neighbours(current->dx, current->dy, current->dz, 
					       i,j,k,
					       &mx, &my, &mz);
*/
		
		general_inverse_transform_point(globals->trans_info.transformation,
						wx,wy,wz,
						&tx,&ty,&tz);
		
		result = optimize_3D_deformation_for_single_node(steps[0], 
								 threshold1,
								 tx,ty,tz,
								 mx, my, mz,
								 &def_x, &def_y, &def_z, 
								 iters, iteration_limit,
								 &nfunks,
								 number_dimensions);
		
	      }
	    }
	      
	    def_x = CONVERT_VALUE_TO_VOXEL(additional_dx, def_x); 
	    SET_VOXEL_3D(additional_dx, i,j,k, def_x);
	    def_y = CONVERT_VALUE_TO_VOXEL(additional_dy, def_y); 
	    SET_VOXEL_3D(additional_dy, i,j,k, def_y);
	    def_z = CONVERT_VALUE_TO_VOXEL(additional_dz, def_z); 
	    SET_VOXEL_3D(additional_dz, i,j,k, def_z);
	    
	  }
	}
      }
      
      
      add_additional_warp_to_current(additional_warp,
				     current_warp,
				     iteration_weight);
      
      /* smooth the warp (result into current->d*)  */
      
      smooth_the_warp(current_warp,
		      additional_warp,
		      additional_mag, (Real)(mean_disp_mag+std));


    }


				/* reset the next iteration's warp. */


    zero = CONVERT_VALUE_TO_VOXEL(additional_dx, 0.0);
    for_less(i,0,sizes[0])
      for_less(j,0,sizes[1])
	for_less(k,0,sizes[2]){
	  SET_VOXEL_3D(additional_dx, i,j,k, zero);
	  SET_VOXEL_3D(additional_dy, i,j,k, zero);
	  SET_VOXEL_3D(additional_dz, i,j,k, zero);

	}
  

    if (globals->flags.debug && 
	globals->flags.verbose == 3)
      save_data(globals->filenames.output_trans, 
		iters+1, iteration_limit, 
		current_warp);
      
    }
    
    /* free up allocated temporary deformation volumes */
  if (globals->trans_info.use_super) {
    (void)free_volume_data(Gsuper_dx); (void)delete_volume(Gsuper_dx);
    (void)free_volume_data(Gsuper_dy); (void)delete_volume(Gsuper_dy);
    (void)free_volume_data(Gsuper_dz); (void)delete_volume(Gsuper_dz);
  }

  (void)free_volume_data(additional_dx); (void)delete_volume(additional_dx);
  (void)free_volume_data(additional_dy); (void)delete_volume(additional_dy);
  (void)free_volume_data(additional_dz); (void)delete_volume(additional_dz);
  (void)free_volume_data(additional_mag); (void)delete_volume(additional_mag);

  free_vector(Ga1xyz ,1,512);
  free_vector(Ga2xyz ,1,512);
  free_vector(TX ,1,512);
  free_vector(TY ,1,512);
  free_vector(TZ ,1,512);
  free_vector(SX ,1,512);
  free_vector(SY ,1,512);
  free_vector(SZ ,1,512);
    

  return (OK);



}




