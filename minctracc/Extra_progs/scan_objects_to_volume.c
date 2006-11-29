/* ----------------------------- MNI Header -----------------------------------
@NAME       : scan_objects_to_volume
@INPUT      : a list of graphical .obj objects and a volume template
@OUTPUT     : a minc volume, with all objects voxelated.
@RETURNS    : 
@DESCRIPTION: makes a label volume from a bunch of objects.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb 13, 1998
@MODIFIED   : 
---------------------------------------------------------------------------- */
#include  <volume_io.h>
#include  <ParseArgv.h>
#include  <time_stamp.h>
#include  <bicpl.h>

void parse_arguments(int argc, char* argv[]);

/* global variables */

VIO_Status     status  = OK;
int        verbose = FALSE;
int        clobber = FALSE;
int        debug   = 0;

int        *labels;
char       **filenames;

char       *template;
char       *output_filename;

int        num_objects;

char       *history;          /* command line added to volume's history */

ArgvInfo argTable[] = {
{"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
     "Show progress"},
{"-clobber", ARGV_CONSTANT, (char *) TRUE, (char *) &clobber,
    "Overwrite output file."},
{"-debug", ARGV_INT, (char *) NULL, (char *) &debug,
    "Show debugging information"},
{NULL, ARGV_END, NULL, NULL, NULL}
};


int  main(
    int   argc,
    char  *argv[] )
{
   VIO_Volume               volume, label_volume;
   File_formats         format;
   int                  object, obj, n_objects, scan_value;
   object_struct        **objects;
   VIO_Real                 max_distance = 1.0;


                                /* parse the flags on the command line and
                                   set
                                      num_objects  

                                   and build up the two arrays: 

                                      filenames[0..num_objects-1]
                                      labels   [0..num_objects-1]
                                */

   parse_arguments(argc, argv);

                                /* read in the template and create a
                                   null label volume.                 */
   
   if( input_volume_header_only(template, 3,
                                File_order_dimension_names, 
                                &volume, NULL) != OK )
      return( 1 );
   
   label_volume = create_label_volume( volume, NC_BYTE );

   set_all_volume_label_data( label_volume, 0 );

   
   for(object=0; object<num_objects; object++) {
      
      if( input_graphics_file(filenames[object],
                              &format, &n_objects, &objects ) != OK )
         return( 1 );
      
      print( "Scanning\n" );
      
      for(obj=0; obj<n_objects; obj++)  {
         scan_object_to_volume(objects[obj],
                               volume, label_volume, labels[object], max_distance );
      }
      
      print( "Done scanning\n" );

   }


   (void) output_volume(output_filename, NC_UNSPECIFIED, FALSE,
                        0.0, 0.0, label_volume, history, NULL );
   
    return( 0 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : parse_arguments
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb. 28, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void parse_arguments(int argc, char* argv[])
{

  int k;      /* volume filename index */
  int m = 2;  /* used to offset class argument */

  /* form the history string - before parseargv is called */
  history = time_stamp(argc, argv);

  /* Call ParseArgv */
  if ( ParseArgv(&argc, argv, argTable, 0) || (argc < 3 )) {
     (void) fprintf(stderr, 
                    "\nUsage: %s <options> template.mnc output.mnc <object> <label> ... \n", 
                    argv[0]);
     (void) fprintf(stderr,   
                    "       %s [-help]\n\n", argv[0]);
     exit(EXIT_FAILURE);
  }

  template        = argv[1];
  output_filename = argv[2];

  if ( debug > 1) {

     fprintf(stdout, "the number of arguments left = %d\n", argc);

     for ( k = 0; k < argc; k++)
        fprintf(stdout, "argv[%d] = %s\n", k, argv[k]);
  }

  
  /* warn of clobbering if output volume file exists */
  if (!clobber && file_exists(output_filename)) {
     printf("File %s exists.\n", output_filename);
     printf("Use -clobber to overwrite.\n");
     exit(EXIT_FAILURE);
  }

  /* make sure each volume has a class label */
  if ( ( argc - 1) % 2 != 0 ) {
     fprintf(stderr, "Please specify a label for each object.\n");   
     exit(EXIT_FAILURE); 
  }
     
  num_objects = ( argc - 2 - 1 ) / 2;  /* count out progname, template and output filename */
 
  if ( debug > 1 ) 
    fprintf(stdout, "the number of objects = %d\n", num_objects);
 
  
  ALLOC( filenames, num_objects );

  ALLOC( labels, num_objects );

  for(k=0; k<num_objects; k++) {
     
     if (!file_exists(argv[(k*m)+3])) {
        
        (void) fprintf(stderr, "filename `%s' not found. \n", argv[(k*m)+3]);
        exit(EXIT_FAILURE);
     }
     
     /* get the filename and class label */
     filenames[k] = argv[(k*m)+3];
     
     labels[k] = atoi(argv[(k*m)+4]);
     
     if ( debug > 1 ) {
        
        fprintf( stdout, "filenames[%d] = %s\n", k, filenames[k]);
        fprintf( stdout, "labels[%d] = %d\n", k, labels[k]);
     }
     
for(=; <; ++)
  

} /* parse_arguments() */


