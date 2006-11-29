void  get_voxel_spatial_loop_limits(VIO_Volume volume,
                                           int start[],        
                                           int end[]);

VIO_BOOL get_average_warp_vector_from_neighbours(VIO_General_transform *trans,
                                                       int voxel[],
                                                       int avg_type,
                                                       VIO_Real *mx, VIO_Real *my, VIO_Real *mz);
VIO_BOOL get_average_warp_of_neighbours(VIO_General_transform *trans,
                                              int voxel[],
                                              VIO_Real mean_pos[]);
void add_additional_warp_to_current(VIO_General_transform *additional,
                                           VIO_General_transform *current,
                                           VIO_Real weight);
void smooth_the_warp(VIO_General_transform *smoothed,
                            VIO_General_transform *current,
                            VIO_Volume warp_mag, VIO_Real thres) ;
void extrapolate_to_unestimated_nodes(VIO_General_transform *current,
                                             VIO_General_transform *additional,
                                             VIO_Volume estimated_flag_vol) ;

VIO_Real get_value_of_point_in_volume(VIO_Real xw, VIO_Real yw, VIO_Real zw, 
                                          VIO_Volume data);

void split_up_the_transformation(VIO_General_transform *total,
                                        VIO_General_transform **all_until_last,
                                        VIO_General_transform **last_warp) ;

void  set_feature_value_threshold(VIO_Volume d1, 
                                         VIO_Volume d2,
                                         VIO_Real *global_thres1, 
                                         VIO_Real *global_thres2, 
                                         VIO_Real *threshold1, 
                                         VIO_Real *threshold2);

void build_two_perpendicular_vectors(VIO_Real orig[], 
                                             VIO_Real p1[], 
                                             VIO_Real p2[]);
float xcorr_objective_with_def(VIO_Volume d1,
                                      VIO_Volume d2,
                                      VIO_Volume m1,
                                      VIO_Volume m2, 
                                      Arg_Data *globals);
