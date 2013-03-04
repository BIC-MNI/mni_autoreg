
#include <bicpl.h>

int main(int argc, char *argv[])
{
  object_struct *obj;
  lines_struct  *lines;
  Point p;
  FILE *fp;

  if (argc!=2) {
    print("usage: %s output.obj\n", argv[0]);
    exit(EXIT_FAILURE);
  }



  obj   = create_object(LINES);
  lines = get_lines_ptr(obj);


  initialize_lines(lines, YELLOW);

  start_new_line(lines);

  fill_Point(p, 0.0, 0.0, 0.0);
  add_point_to_line(lines, &p);
  fill_Point(p, 0.0, 2.0, 1.0);
  add_point_to_line(lines, &p);

  start_new_line(lines);

  fill_Point(p, 0.0, 0.0, 0.0);
  add_point_to_line(lines, &p);
  fill_Point(p, 2.0, 2.0, 0.0);
  add_point_to_line(lines, &p);


  open_file( argv[1] , WRITE_FILE, BINARY_FORMAT, &fp);
  
  output_object(fp, BINARY_FORMAT, obj);

  close_file(fp);

  delete_object(obj);

 exit(EXIT_SUCCESS);
}

