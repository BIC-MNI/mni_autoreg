
#include <bicpl.h>


int point_not_masked(VIO_Volume volume, 
                            VIO_Real wx, VIO_Real wy, VIO_Real wz);

void get_volume_XYZV_indices(VIO_Volume data, int xyzv[]);


int main(int argc, char *argv[])
{
  VIO_General_transform 
    *non_lin_def,
    def;
  VIO_Volume 
    mask,
    vol;

  object_struct 
    *obj;
  lines_struct  
    *lines;
  Point 
    p;
  FILE 
    *fp;
  VIO_Status 
    stat;
  VIO_progress_struct
    progress;

  VIO_Real
     mult,
    voxel[VIO_MAX_DIMENSIONS],
    vx,vy,vz,
    wx,wy,wz,
    dx,dy,dz;
  int 
    i, count,
    index[VIO_MAX_DIMENSIONS],
    xyzv[VIO_MAX_DIMENSIONS],
    sizes[VIO_MAX_DIMENSIONS];


  if (argc<3 || argc>5) {
    print("usage: %s in_def.xfm output.obj [mask.mnc] [mult]\n", argv[0]);
    exit(EXIT_FAILURE);
  }


  stat = input_transform_file(argv[1], &def);
  if (stat != OK) {
    print("error: cannot input %s.\n", argv[1]);
    exit(EXIT_FAILURE);
  }
  
  stat = open_file( argv[2] , WRITE_FILE, BINARY_FORMAT, &fp);
  if (stat != OK) {
    print("error: cannot open %s for output.\n", argv[2]);
    exit(EXIT_FAILURE);
  }
  
  mask = (VIO_Volume)NULL;
  if (argc>=4) {
    stat = input_volume(argv[3], 3, NULL, 
                        NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                        TRUE, &mask, (minc_input_options *)NULL );
    if (stat != OK) {
      print("error: cannot read mask %s.\n", argv[3]);
      exit(EXIT_FAILURE);
    }
    
  }
  mult = 1.0;
  if (argc==5) {
    mult = atof(argv[4]);    
  }

  non_lin_def = (VIO_General_transform *)NULL;
  for(i=0; i<get_n_concated_transforms(&def; i++)) {
    if (get_transform_type( get_nth_general_transform(&def,i) ) 
        == GRID_TRANSFORM)
      non_lin_def = get_nth_general_transform(&def,i);
  }
  
  if (non_lin_def == (VIO_General_transform *)NULL) {
    print("error: cannot find deformation field in %s.\n", argv[1]);
    exit(EXIT_FAILURE);
  }
  
  vol = non_lin_def->displacement_volume;
  
  get_volume_sizes(vol, sizes);
  get_volume_XYZV_indices(vol, xyzv);
  
  obj   = create_object(LINES);
  lines = get_lines_ptr(obj);
  initialize_lines(lines, YELLOW);


  initialize_progress_report(&progress, FALSE, sizes[0]*sizes[1]*sizes[2]+1,
                             "Building vectors");
  count = 0;
  
  for(i=0; i<MAX_DIMENSIONS; i++) index[i] = 0;
  
  for(index[xyzv[VIO_X]]=0; index[xyzv[VIO_X]]<sizes[xyzv[VIO_X]]; index[xyzv[VIO_X]]++)
    for(index[xyzv[VIO_Y]]=0; index[xyzv[VIO_Y]]<sizes[xyzv[VIO_Y]]; index[xyzv[VIO_Y]]++)
      for(index[xyzv[VIO_Z]]=0; index[xyzv[VIO_Z]]<sizes[xyzv[VIO_Z]]; index[xyzv[VIO_Z]]++) {
        

        for(i=0; i<MAX_DIMENSIONS; i++) voxel[i]=index[i];
        convert_voxel_to_world(vol, voxel, &wx, &wy, &wz);

        if (point_not_masked(mask, wx, wy,wz)) {
          index[ xyzv[Z+1] ] = X;
          dx = get_volume_real_value(vol,
                       index[0],index[1],index[2],index[3],index[4]);
          index[ xyzv[Z+1] ] = Y;
          dy = get_volume_real_value(vol,
                       index[0],index[1],index[2],index[3],index[4]);
          index[ xyzv[Z+1] ] = Z;
          dz = get_volume_real_value(vol,
                       index[0],index[1],index[2],index[3],index[4]);

          if ( (fabs(dx) + fabs(dy) + fabs(dz) ) > 0.01 ) {
            start_new_line(lines);
            fill_Point(p, wx, wy, wz);
            add_point_to_line(lines, &p);
            fill_Point(p, wx+dx*mult, wy+dy*mult, wz+dz*mult);
            add_point_to_line(lines, &p);
          }
        }
            
        count ++;
        update_progress_report( &progress, count );
        
      }
  
  terminate_progress_report(&progress);
  
  print ("Saving data...\n");
  
  
  output_object(fp, BINARY_FORMAT, obj);

  close_file(fp);
  
  delete_object(obj);
  
  exit(EXIT_SUCCESS);
}

int point_not_masked(VIO_Volume volume, 
                            VIO_Real wx, VIO_Real wy, VIO_Real wz)
{

  int
    i,
    inside,
    index[VIO_MAX_DIMENSIONS],
    sizes[VIO_MAX_DIMENSIONS];
  VIO_Real
    world[VIO_MAX_DIMENSIONS],
    voxel[VIO_MAX_DIMENSIONS],
    result;

  
  if (volume!=(VIO_Volume)NULL) {
    get_volume_sizes(volume,sizes);
    convert_world_to_voxel(volume, wx, wy, wz, voxel);
    
    inside = TRUE;
    for(i=0; i<get_volume_n_dimensions(volume; i++))
      inside = ( inside && (voxel[i] > -0.5) && (voxel[i] < sizes[i]-0.5));

    if ( inside ) {
      
      for(i=0; i<get_volume_n_dimensions(volume; i++)) {
        index[i] = ROUND( voxel[i] );
      }
      result =  get_volume_real_value(volume,
                       index[0],index[1],index[2],index[3],index[4]);

      if (result>0.0)
        return(TRUE);
      else
        return(FALSE);
    }
    else
      return(FALSE) ;
  }
  else
    return(TRUE) ;
}


void get_volume_XYZV_indices(VIO_Volume data, int xyzv[])
{
  
  int 
    axis, i, vol_dims;
  char 
    **data_dim_names;
  
  vol_dims       = get_volume_n_dimensions(data);
  data_dim_names = get_volume_dimension_names(data);
  
  for(i=0; i<N_DIMENSIONS+1; i++) xyzv[i] = -1;
  for(i=0; i<vol_dims; i++) {
    if (convert_dim_name_to_spatial_axis(data_dim_names[i], &axis )) {
      xyzv[axis] = i; 
    } 
    else {     /* not a spatial axis */
      xyzv[Z+1] = i;
    }
  }
  delete_dimension_names(data,data_dim_names);
  
}



