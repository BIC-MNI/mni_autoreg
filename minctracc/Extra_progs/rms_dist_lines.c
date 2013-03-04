#include <bicpl.h>
#include <config.h>

double minimum_distance_point_to_object (Point           *p, 
                                         object_struct   *targ,
                                         VIO_BOOL         *end_point)
{
  double 
    dx,dy,dz, 
    min_d, dist;
  int 
    j, num_pts_targ;
  Point 
    *the_points_targ;

  num_pts_targ = get_object_points( targ,&the_points_targ);
  *end_point = FALSE;
  

  if (num_pts_targ > 0) {
    min_d = DBL_MAX;

    for(j=0; j<num_pts_targ; j++) { 

      dx =  Point_x(the_points_targ[j]) -  Point_x(*p);
      dy =  Point_y(the_points_targ[j]) -  Point_y(*p);
      dz =  Point_z(the_points_targ[j]) -  Point_z(*p);
      
      dist = dx*dx + dy*dy + dz*dz;
      
      if (dist < min_d) {
          min_d = dist;
          if ((j == 0) || (j == num_pts_targ-1)) 
              *end_point = TRUE;
          else 
          {
              *end_point = FALSE;
          }
          
      }
      
    }
    min_d = sqrt(min_d);
  }
  else
    min_d = DBL_MAX;

  return(min_d);

}

double calc_rms_distance_between_objects(object_struct   *src, 
                                         object_struct   *targ)
{
  int 
    i,j, num_pts,num_pts_src, num_pts_targ;
  double 
    var, std,rms,  min_d_sum, min_d_sum2, min_d;
  VIO_BOOL end_point;
  
  Point *the_points_src, *the_points_targ;

  num_pts_src  = get_object_points( src, &the_points_src);
  num_pts_targ = get_object_points( targ, &the_points_targ);

  rms = 0.0; min_d_sum2 = min_d_sum = 0.0; num_pts = 0;
  

  
  if (num_pts_targ > num_pts_src) {
  if (num_pts_src > 0) {
      for(i=0; i<num_pts_src; i++) {
        
        min_d = minimum_distance_point_to_object(&the_points_src[i], targ, &end_point);
        
        if (!end_point) 
        {

            min_d_sum  += min_d;
            min_d_sum2 += min_d * min_d;
            num_pts++;
        }
        
      }
      
/*      rms = sqrt( min_d_sum2 / num_pts_src );*/
    }
  }
  else {
    if (num_pts_targ > 0) {
      for(i=0; i<num_pts_targ; i++) {
        
        min_d = minimum_distance_point_to_object(&the_points_targ[i], src, &end_point);
        
        if (!end_point) 
        {
            min_d_sum  += min_d;
            min_d_sum2 += min_d * min_d;
            num_pts++;
        }
        
      }
      
    }
  }

  if (num_pts>0) 
  {
      rms = min_d_sum / num_pts ; /* this should use min_d_sum2 for RMS */
  }
  
  return(rms);


}

void apply_transform_to_object(VIO_General_transform *xform, 
                               object_struct   *object )
{
  int i,num_points;
  Point *the_points;
  VIO_Real tx,ty,tz;

  num_points = get_object_points( object, &the_points);
  for(i=0; i<num_points; i++) {
    general_transform_point(xform, 
                            Point_x(the_points[i]),
                            Point_y(the_points[i]),
                            Point_z(the_points[i]),
                            &tx, &ty, &tz);
    Point_x(the_points[i]) = tx; 
    Point_y(the_points[i]) = ty; 
    Point_z(the_points[i]) = tz; 
  }
}


int main(int argc, char *argv[])
{
  VIO_General_transform 
    xform;

  object_struct 
    **list_of_objs,**list_of_objs2;
  
  Object_types   type;

  File_formats   
    format;
  VIO_progress_struct
    progress;

  int 
    i,
    count,
    num_objects,
    num_objects2;

  double dist;

  if (argc != 4) {
    print("usage: %s xform.xfm source.obj target.obj \n", argv[0]);
    exit(EXIT_FAILURE);
  }

  if (input_transform_file(argv[1], &xform) != OK) {
    print("error: cannot input %s.\n", argv[1]);
    exit(EXIT_FAILURE);
  }
  
  format = BINARY_FORMAT;
  if (input_graphics_file(argv[2], 
                          &format,
                          &num_objects,
                          &list_of_objs) != OK) {
    print("error: problems reading %s.\n", argv[2]);
    exit(EXIT_FAILURE);
  }
  if (input_graphics_file(argv[3], 
                          &format,
                          &num_objects2,
                          &list_of_objs2) != OK) {
    print("error: problems reading %s.\n", argv[3]);
    exit(EXIT_FAILURE);
  }

  for(i=0; i<num_objects; i++) {
    type = get_object_type( list_of_objs[i] );
    switch (type) {
    case MARKER:print ("markers unsupported, sorry\n"); break;
    case MODEL:print ("models unsupported, sorry\n"); break;
    case PIXELS:print ("pixels unsupported, sorry\n"); break;
    case LINES: 
      apply_transform_to_object(&xform,  list_of_objs[i] );
      break;
    case POLYGONS:
      print ("polygons found, will transform coords\n"); 
      apply_transform_to_object(&xform,  list_of_objs[i] );
      break;
    case QUADMESH:
      print ("quadmesh found, will transform coords\n"); 
      apply_transform_to_object(&xform,  list_of_objs[i] );
      break;
    case TEXT:print ("text unsupported, sorry\n"); break;
    default:
      print_error( "Unrecognized object type %d\n", type );
    }
  }

  dist = calc_rms_distance_between_objects(list_of_objs[0], list_of_objs2[0]);

  if (dist!=0.0)
      print ("%f\n",dist);

  
  
  for(i=0; i<num_objects; i++)
    delete_object(list_of_objs[i]);

  for(i=0; i<num_objects2; i++)
    delete_object(list_of_objs2[i]);
  
  exit(EXIT_SUCCESS);
}




















