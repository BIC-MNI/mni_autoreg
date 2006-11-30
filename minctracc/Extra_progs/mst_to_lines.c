
#include <bicpl.h>

typedef struct {          /* element of minimum spanning tree */
   int index;
   int parent_index;
   float xyz[3];
   float err;
} ML;

#define NPTS 50000

int main(int argc, char *argv[])
{

  static ML mst[NPTS];              /* mst of at most NPTS points */

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

  long int
    i, count, total;

  if (argc!=3) {
    print("usage: %s in_mst output.obj\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  stat = open_file( argv[1] , READ_FILE, ASCII_FORMAT, &fp);
  if (stat != OK) {
    print("error: cannot open %s for input.\n", argv[1]);
    exit(EXIT_FAILURE);
  }

  count = 0L;
  total = 0L;
  while( fscanf( fp, "%d%d%f%f%f%f", 
                &mst[total].index,
                &mst[total].parent_index,
                &mst[total].xyz[0],
                &mst[total].xyz[1],
                &mst[total].xyz[2],
                &mst[total].err) != EOF ) {
    ++total;
    if ( total >= NPTS ) {
      printf("\ntoo much data!");
      fclose( fp );
      exit( 0 );
    }
  }
  if ( close_file( fp ) != OK ){
    print("error: cannot close %s.\n", argv[1]);
    exit(EXIT_FAILURE);
  }

  
  stat = open_file( argv[2] , WRITE_FILE, BINARY_FORMAT, &fp);
  if (stat != OK) {
    print("error: cannot open %s for output.\n", argv[2]);
    exit(EXIT_FAILURE);
  }
  
  obj   = create_object(LINES);
  lines = get_lines_ptr(obj);
  initialize_lines(lines, YELLOW);


  initialize_progress_report(&progress, FALSE, total+1,
                             "Building vectors");
  


  for(i=1; i<total; i++) {

  start_new_line(lines);
  
  count = mst[i].index;
  fill_Point(p, mst[count].xyz[0],mst[count].xyz[1],mst[count].xyz[2]);
  add_point_to_line(lines, &p);

  count = mst[i].parent_index;
  fill_Point(p, mst[count].xyz[0],mst[count].xyz[1],mst[count].xyz[2]);
  add_point_to_line(lines, &p);
  
  update_progress_report( &progress, i );
        
  }
  
  terminate_progress_report(&progress);
  
  print ("Saving data...\n");
  
  
  output_object(fp, ASCII_FORMAT, obj);

  close_file(fp);
  
  delete_object(obj);
  
  exit(EXIT_SUCCESS);
}

