/* ----------------------------- MNI Header -----------------------------------
@NAME       :  mincchamfer.c 
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: read in a mnc volume, and compute the chamfer distance
              
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Wed Jul  7 20:50:20 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */


#define XSIZE 256
#define YSIZE 256

#include "mincchamfer.h"
#include "chamfer.h"

int main ( int argc, char* argv[] )
{   
   static char *default_dim_names[VIO_N_DIMENSIONS] = { MIzspace, MIyspace, MIxspace };

   char 
      *infilename,
      *outfilename;
   
   VIO_Volume 
      data;
   
   FILE 
      *ifd,*ofd;
   
   VIO_Status 
      status;

   char* history = history_string( argc, argv );

   /* set globals */
   max_dist = 50;
   debug = FALSE;
   verbose = TRUE;
   prog_name = argv[0];
   
   
   /* Call ParseArgv to interpret all command line args */
   
   if (ParseArgv(&argc, argv, argTable, 0) || (argc!=3)) {
      (void) fprintf(stderr,
                     "\nUsage: %s [<options>] <inputfile> <outputfile>\n",
                     prog_name);
      (void) fprintf(stderr,"       %s [-help]\n\n", prog_name);
      exit(EXIT_FAILURE);
   }
   
   infilename  = argv[1];        /* set up necessary file names */
   outfilename = argv[2];
   
   if (debug) {
      printf ("input file :%s\n",infilename);
      printf ("output file:%s\n",outfilename);
   }
   
   
                                /* check files to be used.  */
   
   status = open_file( infilename , READ_FILE, BINARY_FORMAT, &ifd );
   if ( status != OK ) 
      print_error_and_line_num ("filename `%s' not found.", __FILE__, __LINE__, infilename);
   status = close_file(ifd);
   
   status = open_file( outfilename , WRITE_FILE, BINARY_FORMAT, &ofd );
   if ( status != OK ) 
      print_error_and_line_num ("filename `%s' cannot be opened.", __FILE__, __LINE__, outfilename);
   status = close_file(ofd);
   remove(outfilename);
   
   
   /* read input data volume */
   
   if (verbose)
      printf ("Reading in input data.\n");
   status = input_volume(infilename, 3, default_dim_names, 
                         NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                         TRUE, &data, (minc_input_options *)NULL );
   
   if ( status != OK ) 
      print_error_and_line_num ("Problems reading `%s'.", __FILE__, __LINE__, infilename);
   
   
                                /* call the chamfer estimation */
   
   status = compute_chamfer(data,max_dist);
   if (status != OK)
      print_error_and_line_num("problems computing chamfer...",__FILE__, __LINE__);
   
                                /* output volume data!  */

   if (verbose)
      printf ("Writing volumetric data.\n");
   if (debug)
      printf ("output file:%s\n",outfilename);
   
   status = output_modified_volume(outfilename, NC_UNSPECIFIED, FALSE, 0.0, 0.0, data, 
                                   infilename, history, (minc_output_options *)NULL);
   
   if (status != OK)
      print_error_and_line_num("problems writing chamfer...",__FILE__, __LINE__);
   
   
   delete_volume(data);
   

   return(OK);
}

