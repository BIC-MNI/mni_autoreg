/* ----------------------------- MNI Header -----------------------------------
@NAME       : crispify
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: makes a crisp labelled volume out of fuzzy ones
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Aug 25, 1996 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
#include <volume_io.h>
#include <ParseArgv.h>
#include <time_stamp.h>

/* function prototypes */

void parse_arguments(int argc, char* argv[]);
void load_volume(char *, VIO_Volume * );
void create_empty_crisp_volume(VIO_Volume volume_example);
void scan_fuzzy_volumes(VIO_Volume in_vol, char* label, VIO_Volume max_vol, VIO_Volume crisp_vol);
void write_crisp_volume(void);

/* global variables */

VIO_Status     status;
int        verbose = FALSE;
int        clobber = FALSE;
int        debug = 0;

char       **fuzzy_label;
char       **input_filename;
int        fuzzy_vol_index = 0;           /* index into fuzzy vol array */

VIO_Volume     in_volume, max_volume, crisp_volume;
int        num_volumes;
int        max_class_index = 0;           /* highest label value */
char       *crisp_vol_filename;           /* name of the crisp volume */
int        fuzz_vol;                      /* index to  number of fuzzy volumes */
int        block_sizes[3] = { 1, -1, -1}; /* set the block size for slice access */

int        first_volume_sizes[VIO_MAX_DIMENSIONS];       /* 1D array to hold sizes for 1st vol */
int        first_volume_num_dims;     /* to hold num of dimensions */
char       **first_volume_dim_names;  /* 1D array of char* to hold the dim names */
char       *history;          /* command line added to volume's history */

ArgvInfo argTable[] = {


  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
     "Show progress"},
  
  {"-clobber", ARGV_CONSTANT, (char *) TRUE, (char *) &clobber,
     "Overwrite output file."},
  
  {"-debug", ARGV_INT, (char *) NULL, (char *) &debug,
     "Show debugging information"},

  {"-volume", ARGV_STRING, (char *) NULL, (char *) &crisp_vol_filename,
     "Specify a volume to write crisp labels."},

  
  {NULL, ARGV_END, NULL, NULL, NULL}
};



int main(int argc, char *argv[])
{
   
   int i,j,k;
   
                                /* parse the flags on the command line and
                                   set
                                      num_volumes  

                                   and build up the two arrays: 

                                      input_filename[0..num_volumes-1]
                                      fuzzy_label   [0..num_volumes-1]
                                */
   parse_arguments(argc, argv);

                                /* get first volume into memory */
   fuzzy_vol_index = 0;
   load_volume( input_filename[0], &in_volume );

                                /* init the discrete crisp volume */
   create_empty_crisp_volume(in_volume);
   for(i=0; i<first_volume_sizes[0]; i++)
      for(j=0; j<first_volume_sizes[1]; j++)
         for(k=0; k<first_volume_sizes[2]; k++)
            set_volume_real_value(crisp_volume,
                                  i,j,k,0,0, 0.0);

                                /* init the max volume */
   max_volume = copy_volume_definition(in_volume,
                                       NC_FLOAT,
                                       FALSE,
                                       0.0, 0.0 );
   for(i=0; i<first_volume_sizes[0]; i++)
      for(j=0; j<first_volume_sizes[1]; j++)
         for(k=0; k<first_volume_sizes[2]; k++)
            set_volume_real_value(max_volume,
                                  i,j,k,0,0, 0.0);


                                /* write in the first volume labels */
   scan_fuzzy_volumes(in_volume, fuzzy_label[0], max_volume, crisp_volume);

  /* allocate memory for first volume sizes */
   
   for(fuzz_vol=1; fuzz_vol<num_volumes; fuzz_vol++) {

      fuzzy_vol_index++;
      load_volume( input_filename[fuzz_vol], &in_volume );
      scan_fuzzy_volumes(in_volume, fuzzy_label[fuzz_vol], max_volume, crisp_volume);
   }
   

   write_crisp_volume();
   
   delete_volume(in_volume);
   delete_volume(max_volume);
   delete_volume(crisp_volume);
   
   exit(EXIT_SUCCESS);

} /* main */





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
  if ( ParseArgv(&argc, argv, argTable, 0) || (argc < 2 )) {
    (void) fprintf(stderr, 
                   "\nUsage: %s <options> <infile> <class> ... \n", 
                   argv[0]);
    (void) fprintf(stderr,   
                   "       %s [-help]\n\n", argv[0]);
    exit(EXIT_FAILURE);
  }


  if ( debug > 1) {

    fprintf(stdout, "the number of arguments left = %d\n", argc);

    for ( k = 0; k < argc; k++)
      fprintf(stdout, "argv[%d] = %s\n", k, argv[k]);
  }

  /* make sure you specify a new tag filename */
  if (!crisp_vol_filename) {
    
    printf("Please specify crisp volume filename with -crisp_volume  switch\n");
    exit(EXIT_FAILURE);
  }
  
  /* warn of clobbering if crisp volume file exists */
  if (!clobber && file_exists(crisp_vol_filename)) {
    
    printf("File %s exists.\n", crisp_vol_filename);
    printf("Use -clobber to overwrite.\n");
    exit(EXIT_FAILURE);
  }

  /* make sure each volume has a class label */
  if ( ( argc - 1) % 2 != 0 ) {

    fprintf(stderr, "Please specify a class label.\n");   
    exit(EXIT_FAILURE); 
  }
     
  num_volumes = ( argc - 1 ) / 2;  /* count out progname */
 
  if ( debug > 1 ) 
    fprintf(stdout, "the number of volumes = %d\n", num_volumes);
 
  
  ALLOC( input_filename, num_volumes );

  ALLOC( fuzzy_label, num_volumes );

  for(k=0; k<num_volumes; k++) {
    
    if (!file_exists(argv[(k*m)+1])) {

      (void) fprintf(stderr, "filename `%s' not found. \n", argv[(k*m)+1]);
      exit(EXIT_FAILURE);
    }

    /* get the filename and class label */
    input_filename[k] = argv[(k*m)+1];

    fuzzy_label[k] = argv[(k*m)+2];
 
    if ( debug > 1 ) {

      fprintf( stdout, "input_filename[%d] = %s\n", k, input_filename[k]);
      fprintf( stdout, "fuzzy_label[%d] = %s\n", k, fuzzy_label[k]);
    }
    
  }

  /* determine the maximum label value input by user */
  for(k=0; k<num_volumes; k++) {

    if ( max_class_index < atoi(fuzzy_label[k]) ) 
        max_class_index = atoi(fuzzy_label[k]);
  }

} /* parse_arguments() */



/* ----------------------------- MNI Header -----------------------------------
@NAME       : load_volume
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
void load_volume(char *in_filename, VIO_Volume *volume)
{


  if (verbose)
    printf ("Processing volume %s\n", in_filename);



  /* load the volume */
  status = input_volume( in_filename, 3, NULL, NC_UNSPECIFIED, 
                        FALSE, 0.0, 0.0,
                        (fuzzy_vol_index >= 1 ? FALSE: TRUE), volume, (minc_input_options *) NULL ) ;

  if( status != OK )
    exit(EXIT_FAILURE);

  /* if loading the very first volume, get its sizes, number of dims and dim 
     names, dim starts and dim steps */
  if ( fuzzy_vol_index == 0 ) {
    
    int k; /* local counter */
    
    get_volume_sizes(*volume, first_volume_sizes);
    first_volume_num_dims = get_volume_n_dimensions(*volume);
    first_volume_dim_names = get_volume_dimension_names(*volume);
    
    if ( debug > 2 ) {
      fprintf(stdout, "Vol number of dims. = %d\n", first_volume_num_dims);
      
      fprintf(stdout, "Vol dimension names = ");
      for(k=0; k<first_volume_num_dims; k++) 
        fprintf(stdout, "%s ", first_volume_dim_names[k]);
      fprintf(stdout, "\n");
      
    }
    
  }

  /* if you have more than one volume, check to see if volumes are
     of same size in each dim and make sure the dimension orders are
     also the same. this is done by volume_size_is_ok() function */
  
  if ( fuzzy_vol_index >= 1 &&  !volume_size_is_ok( *volume )){
    
    (void) fprintf(stderr,"in volume %s\n", in_filename);
    exit(EXIT_FAILURE);
  }
  


} /* load_volume(...) */




/* ----------------------------- MNI Header -----------------------------------
@NAME       : scan_fuzzy_volumes
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: scans fuzzy volumes to find highest fuzzy label
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Aug 25, 1996 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void scan_fuzzy_volumes(VIO_Volume in_vol, char* label, VIO_Volume max_vol, VIO_Volume crisp_vol)
{

  int v1, v2, v3, j;

  VIO_Real max_fuzzy_val;      /* var to hold the highest fuzzy value */
  int  class_val;          /* var to hold the corresponding class label */ 
  VIO_Real fuzzy_value;        /* var to hold the current fuzzy value */


  class_val = atoi(label);

  if (verbose)
    (void) fprintf(stdout, "Scanning fuzzy  volume for label %d \n", class_val);


  /* take the size of volume 0, since they should all be the same */
  for(v1=0; v1<first_volume_sizes[0]; v1++) {

    if ( verbose ) 
      write(2,"*",1);

    for(v2=0; v2<first_volume_sizes[1]; v2++) {

      if ( debug > 10 )
        write(2,"-",1);
      
      for(v3=0; v3<first_volume_sizes[2]; v3++) {

        if ( debug >= 20 )
          write(2,"+",1);

        fuzzy_value   = get_volume_real_value(in_vol,  v1, v2, v3, 0, 0);
        
        max_fuzzy_val = get_volume_real_value(max_vol, v1, v2, v3, 0, 0);

        
        if ( max_fuzzy_val <  fuzzy_value ) {
           
           set_volume_real_value(crisp_vol, v1, v2, v3, 0, 0, (VIO_Real)class_val);
           set_volume_real_value(max_vol,   v1, v2, v3, 0, 0, fuzzy_value);
                                
        }
        
      } /* v3 */

    } /* v2 */

  } /* v1 */

} /* scan_fuzzy_volumes(...)*/



/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_empty_crisp_volume
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRION   : create and initialized an empty crisp volume
@METHOD     : 
@GLOBAL     : 
@CALLS      : 
@CREATED    : February 6, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
-------------------------------------------------------------------------- */
void create_empty_crisp_volume(VIO_Volume volume_example)
{

  /* create the classification volume here */   

  if (verbose) 
    fprintf( stdout ,"Creating crisp output volume\n");


 
  crisp_volume = copy_volume_definition(volume_example,
                                        NC_BYTE,
                                        FALSE,
                                        0.0, (VIO_Real) max_class_index );

  set_volume_voxel_range(crisp_volume, 0.0, (VIO_Real) max_class_index );
  set_volume_real_range(crisp_volume,  0.0, (VIO_Real) max_class_index );



} /* create_empty_crisp_volume */


/* ----------------------------- MNI Header -----------------------------------
@NAME       : volume_size_is_ok
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: verifies that volume sizes are OK
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb 10, 1996 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
int volume_size_is_ok( VIO_Volume loaded_volume) 
{


  int    *loaded_volume_sizes;
  int    loaded_volume_num_dims;
  STRING *loaded_volume_dim_names;
  
  /* allocate memory for first volume sizes */
  ALLOC(loaded_volume_sizes, VIO_MAX_DIMENSIONS); 

  /* get dim size, nums, order */
  get_volume_sizes(loaded_volume, loaded_volume_sizes);
  loaded_volume_num_dims = get_volume_n_dimensions(loaded_volume);
  loaded_volume_dim_names = get_volume_dimension_names(loaded_volume);
  
  if ( debug > 2 ) {

    int k; /* local counter */      
    
    fprintf(stdout, "Vol number of dims. = %d\n", loaded_volume_num_dims);
    
    fprintf(stdout, "Vol dimension names = ");
    for(k=0; k<loaded_volume_num_dims; k++) 
      fprintf(stdout, "%s ", loaded_volume_dim_names[k]);
    fprintf(stdout, "\n");
    
  }

  /* all the volume dimensions should be the same as the first volume */

  /* check for number of dimensions mismatch */
  if (loaded_volume_num_dims != first_volume_num_dims ) {

    (void) fprintf(stderr,"Error - Number of dimensions mismatch ");
    return FALSE;    
  }

     
  /* check for volume size mismatches */
  if (loaded_volume_sizes[VIO_X] != first_volume_sizes[VIO_X]) {

    (void) fprintf(stderr,"Error - VIO_Volume size mismatch in X dimension ");
    return FALSE;
  }

  if (loaded_volume_sizes[VIO_Y] != first_volume_sizes[VIO_Y]) {

    (void) fprintf(stderr,"Error - VIO_Volume size mismatch in Y dimension ");
    return FALSE;
  }

  if (loaded_volume_sizes[VIO_Z] != first_volume_sizes[VIO_Z]) {

    (void) fprintf(stderr,"Error - VIO_Volume size mismatch in Z dimension ");
    return FALSE;
  }


  /* check for dimensions order mismatch */
  if ( strcmp(loaded_volume_dim_names[0], first_volume_dim_names[0]) ) {
      
    (void) fprintf(stderr,"Error - First dimension order mismatch ");
    return FALSE;
  }

  /* if there are more than 1 dimension - check dim order of second dim */
  if ( loaded_volume_num_dims > 1) 
    if ( strcmp(loaded_volume_dim_names[1], first_volume_dim_names[1]) ) {
      
      (void) fprintf(stderr,"Error - Second dimension order mismatch ");
      return FALSE;
    }

  /* if there are more than 2 dimensions - check dim order of third dim*/
  if ( loaded_volume_num_dims > 2) 
    if ( strcmp(loaded_volume_dim_names[2], first_volume_dim_names[2]) ) {
      
      (void) fprintf(stderr,"Error - Third dimension order mismatch ");
      return FALSE;
    }

  /* if there are more then 3 dimensions - warn and die */
  if ( loaded_volume_num_dims > 3) {

    (void) fprintf(stderr,"Support is limited to 3 spatial dimensions ");
    exit(EXIT_FAILURE);
  }

  /* free the reserved memory locations */
  delete_dimension_names( loaded_volume, loaded_volume_dim_names );
  FREE(loaded_volume_sizes);

  return TRUE;
            
} /* volume_size_is_ok */
        

/* ----------------------------- MNI Header -----------------------------------
@NAME       : write_crisp_volume()
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: write out the crisp volume into a minc file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : February 6, 1995 ( Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void write_crisp_volume(void)
{

  /* write the crisp volume to a file */
  if (verbose) 
    fprintf(stdout,"\nWriting crisp volume %s to file ...\n", crisp_vol_filename);
   


  status = output_volume( crisp_vol_filename, 
                         NC_BYTE, 
                         FALSE, 
                         0.0, 0.0,
                         crisp_volume, 
                         history, 
                         (minc_output_options *) NULL ) ;


  if ( status != OK )
     exit(EXIT_FAILURE);
                              

} /* write_crisp_volume */

