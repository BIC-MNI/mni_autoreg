

#include <bicpl.h>


void apply_transform_to_polygons(VIO_General_transform *xform, 
                                 object_struct   *polygon_object )
{
  int i,num_points;
  Point *the_points;
  VIO_Real tx,ty,tz;

  num_points = get_object_points( polygon_object, &the_points);
  print ("There are %d points to transform.\n", num_points);
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
    **list_of_objs;                /* should become an array of
                                   object_struct pointers    */
  Object_types   type;

  File_formats   
    format;
  VIO_progress_struct
    progress;

  int 
    i,
    count,
    num_objects;

  if (argc != 4) {
    print("usage: %s xform.xfm input.obj output.obj \n", argv[0]);
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

/*
  initialize_progress_report(&progress, FALSE, 1000,
                             "Building vectors");
  count = 0;
  count ++;
  update_progress_report( &progress, count );
  terminate_progress_report(&progress);
*/
  
  for(i=0; i<num_objects; i++) {
    print ("object %d of %d:\n",i+1,num_objects);
    type = get_object_type( list_of_objs[i] );
    switch (type) {
    case LINES: print ("lines unsupported, sorry\n"); break;
    case MARKER:print ("markers unsupported, sorry\n"); break;
    case MODEL:print ("models unsupported, sorry\n"); break;
    case PIXELS:print ("pixels unsupported, sorry\n"); break;
    case POLYGONS:
      print ("polygons found, will transform coords\n"); 
      apply_transform_to_polygons(&xform,  list_of_objs[i] );
      break;
    case QUADMESH:print ("quadmesh unsupported, sorry\n"); break;
    case TEXT:print ("text unsupported, sorry\n"); break;
    default:
      print_error( "Unrecognized object type %d\n", type );
    }
  }

  print ("Saving data...\n");
  if (output_graphics_file(argv[3], format, num_objects, list_of_objs) != OK) {
    print("error: problems saving to %s.\n", argv[3]);
    exit(EXIT_FAILURE);
  }
  
  
  for(i=0; i<num_objects; i++)
    delete_object(list_of_objs[i]);
  
  exit(EXIT_SUCCESS);
}

