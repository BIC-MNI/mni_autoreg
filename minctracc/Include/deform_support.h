/*  this file was created automatically with:
 
extract_functions -public Optimize/deform_support.c > Include/deform_support.h
 
   and these lines added manually!
*/

public void  get_voxel_spatial_loop_limits(Volume volume,
					   int start[],	
					   int end[]);

public BOOLEAN get_average_warp_vector_from_neighbours(General_transform *trans,
						       int voxel[],
						       int avg_type,
						       Real *mx, Real *my, Real *mz);
public BOOLEAN get_average_warp_of_neighbours(General_transform *trans,
					      int voxel[],
					      Real mean_pos[]);
public void add_additional_warp_to_current(General_transform *additional,
					   General_transform *current,
					   Real weight);
public void smooth_the_warp(General_transform *smoothed,
			    General_transform *current,
			    Volume warp_mag, Real thres) ;
public void extrapolate_to_unestimated_nodes(General_transform *current,
					     General_transform *additional,
					     Volume estimated_flag_vol) ;

public Real get_value_of_point_in_volume(Real xw, Real yw, Real zw, 
					  Volume data);

public void split_up_the_transformation(General_transform *total,
					General_transform **all_until_last,
					General_transform **last_warp) ;

public void  set_feature_value_threshold(Volume d1, 
					 Volume d2,
					 Real *global_thres1, 
					 Real *global_thres2, 
					 Real *threshold1, 
					 Real *threshold2);

public void build_two_perpendicular_vectors(Real orig[], 
					     Real p1[], 
					     Real p2[]);
public float xcorr_objective_with_def(Volume d1,
                                      Volume d2,
                                      Volume m1,
                                      Volume m2, 
                                      Arg_Data *globals);
