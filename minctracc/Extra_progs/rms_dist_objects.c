/*
# ------------------------------ MNI Header ----------------------------------
#@NAME       : rms_dist_objects
#@INPUT      : one transform, two objects
#@OUTPUT     : 
               the rms minimum distance for all points on the first object
               (after transformation) and the second object.
#@RETURNS    : 0 if no error, VIO_Status otherwise
#@DESCRIPTION: 
#@CREATED    : Tuesday, Feb 10, 1998, Louis Collins
#@MODIFIED   : from rms_dist_lines
#@VERSION    : $Id: rms_dist_objects.c,v 1.3 2006-11-29 09:09:31 rotor Exp $
#-----------------------------------------------------------------------------
*/

#include <bicpl.h>
#include <config.h>


VIO_BOOL end_point(object_struct   *object, 
                  int             obj_index)
{
   VIO_BOOL is_an_end;
   int     i,j,m,n;
   
   is_an_end = FALSE;


    switch (get_object_type( object )) {
    case MARKER: print ("end_point(): markers unsupported, sorry\n"); break;
    case MODEL:  print ("end_point(): models unsupported, sorry\n"); break;
    case PIXELS: print ("end_point(): pixels unsupported, sorry\n"); break;
    case TEXT:   print ("end_point(): text unsupported, sorry\n"); break;
    case LINES: 
       get_lines_ptr(object);
      break;
    case POLYGONS:
      is_an_end = FALSE;
      break;
    case QUADMESH:
      get_quadmesh_n_objects( get_quadmesh_ptr(object), &m, &n );
      i = obj_index / n;
      j = obj_index % n;

/*      print ("%3d %3d %3d %3d ", i,j,m,n);*/

      is_an_end = (i==0 || j==0 || i ==(m-1) || j ==(n-1));
      break;
    default:
      print_error( "Unrecognized object type %d\n", get_object_type( object ) );
    }

   return(is_an_end);
}




double calc_rms_distance_between_objects(object_struct   **src, int num_src ,
                                         object_struct   **targ, int num_targ)
{
   VIO_Real 
      dist;
   double 
      var, std,rms,  min_d_sum, min_d_sum2, min_d, d;

   Point 
      *the_points_src,  point, closest_point;

   object_struct 
      *closest_targ;

   int 
      i,j,num_pts, num_pts_src, s,t, object_index,
      closest_t, closest_index;

   VIO_progress_struct
      progress;

   VIO_BOOL flag;

   rms = 0.0; min_d_sum2 = min_d_sum = 0.0; num_pts = 0;



   for(s=0; s<num_src; s++) {      /* for all objects in the source list */

      num_pts_src  = get_object_points( src[s], &the_points_src);

      initialize_progress_report(&progress, TRUE, 
                                 num_pts_src+1,
                                 "RMS" );

      print ("num_pts_src = %d\n",num_pts_src);

      for(i=0; i<num_pts_src; i++) { /* for all points in the object */

         min_d = 1e20;

         for(t=0; t<num_targ; t++) { /* search through all targ objects */

            d = find_closest_point_on_object(&the_points_src[i], targ[t],
                                             &object_index, &point);
            if (d<min_d) {
               min_d = d;
               closest_t = t;
               closest_targ = targ[t];
               closest_point = point;
               closest_index = object_index;
            }
      
         }

/*print ("[%2d] %7.2f %7.2f %7.2f -> (%7.2f)  [%2d/%2d] %7.2f %7.2f %7.2f ",
       i,
       Point_x(the_points_src[i]),
       Point_y(the_points_src[i]),
       Point_z(the_points_src[i]),
       min_d,
       closest_t,
       closest_index,
       Point_x(closest_point),
       Point_y(closest_point),
       Point_z(closest_point));
*/
flag = end_point(closest_targ,closest_index);
/*
if (flag) {
   print ("end\n");
}
else {
   print ("\n");
}
*/


         if (min_d < 1e20 && ! flag) {
            min_d_sum  += min_d;
            min_d_sum2 += min_d * min_d;
            num_pts++;
         }
         update_progress_report( &progress, i);
      }
      terminate_progress_report( &progress );


   }

   print ("num_pts = %d\n",num_pts);

   if (num_pts>0) {
      rms = sqrt(min_d_sum2 / num_pts) ; 
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
     switch (get_object_type(list_of_objs[i])) {
        case MARKER:   print ("markers \n"); break;
        case MODEL:    print ("models \n"); break;
        case PIXELS:   print ("pixels\n"); break;
        case LINES:    print ("lines\n"); break;
        case POLYGONS: print ("polygons\n"); break;
        case QUADMESH: print ("quadmesh\n"); break;
        case TEXT:     print ("text\n"); break;
        default:
        print_error( "Unrecognized object type %d\n", 
                    get_object_type(list_of_objs[i]) );
     }
     apply_transform_to_object(&xform,  list_of_objs[i] );
  }

  for(i=0; i<num_objects2; i++) {
     if (get_object_type(list_of_objs2[i]) == QUADMESH) {
        create_quadmesh_bintree(get_quadmesh_ptr(list_of_objs2[i]) ,200);
     }
  }


  dist = calc_rms_distance_between_objects(list_of_objs, num_objects ,
                                           list_of_objs2,num_objects2);

  print ("%f\n",dist);

  
  
  for(i=0; i<num_objects; i++)
    delete_object(list_of_objs[i]);

  for(i=0; i<num_objects2; i++)
    delete_object(list_of_objs2[i]);
  
  exit(EXIT_SUCCESS);
}




















