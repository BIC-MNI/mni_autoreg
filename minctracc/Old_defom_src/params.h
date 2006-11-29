/* flags */
         if ( strncmp(argv[argn],"-fixed_scale",MAX(strlen(argv[argn]),3)) == 0 ) {
            fixed_scale_flag = TRUE;
            argn--;
         } else

         if ( strncmp(argv[argn],"-invert",MAX(strlen(argv[argn]),4)) == 0 ) {
            invert_mapping_flag = TRUE;
            argn--;
         } else

/* filenames */
         if ( strncmp(argv[argn],"-bmask",MAX(strlen(argv[argn]),4)) == 0 ) {
            if ( (mask2filename = copy_string( argv[argn+1] )) == NULL )
               { printf ("(-bmask)"); print_usage( usage ); }
            status = open_file( mask2filename , READ_FILE, BINARY_FORMAT,  &bfd );
            if ( status != OK ) 
               print_error ("mask2 filename `%s' not found.", __FILE__, __LINE__, mask2filename, 0,0,0,0);
         } else
         if ( strncmp(argv[argn],"-infile",MAX(strlen(argv[argn]),4)) == 0 ) {
            if ( (infilename = copy_string( argv[argn+1] )) == NULL )
               { printf ("(-infile)"); print_usage( usage ); }
            status = open_file( infilename , READ_FILE, BINARY_FORMAT,  &ifd );
            if ( status != OK ) 
               print_error ("input filename `%s' not found.", __FILE__, __LINE__, infilename, 0,0,0,0);
         } else
         if ( strncmp(argv[argn],"-modelfile",MAX(strlen(argv[argn]),3)) == 0 ) {
           if ( (modelfilename = copy_string( argv[argn+1] )) == NULL )
             { printf ("(-model)"); print_usage( usage ); }
           status = open_file( modelfilename ,READ_FILE, BINARY_FORMAT, &mfd );
           if ( status != OK ) 
             print_error ("model filename `%s' not found.", __FILE__, __LINE__, modelfilename, 0,0,0,0);
         } else
         if ( strncmp(argv[argn],"-outfile",MAX(strlen(argv[argn]),2)) == 0 ) {
           if ( (outfilename = copy_string( argv[argn+1] )) == NULL )
             { printf ("(-outf)"); print_usage( usage ); }
           status = open_file( outfilename ,WRITE_FILE, BINARY_FORMAT, &ofd );
           if ( status != OK ) 
             print_error ("output filename `%s' cannot be opened.", __FILE__, __LINE__, outfilename, 0,0,0,0);
         } else

 /* data rotations: */
         if ( strncmp(argv[argn],"-rotat",MAX(strlen(argv[argn]),3)) == 0 ) {
            if ( sscanf( argv[argn+1], "%f", &data_rot_x) != 1 )
               { printf ("(-rotat rx)"); print_usage( usage ); }
            if ( sscanf( argv[argn+2], "%f", &data_rot_y) != 1 )
               { printf ("(-rotat ry)"); print_usage( usage ); }
            if ( sscanf( argv[argn+3], "%f", &data_rot_z) != 1 )
               { printf ("(-rotat rz)"); print_usage( usage ); }
            argn += 2;
            rot_flag = TRUE;
         } else
 /* data scaling factors: */
         if ( strncmp(argv[argn],"-scale",MAX(strlen(argv[argn]),3)) == 0 ) {
            if ( sscanf( argv[argn+1], "%f", &data_scale_x) != 1 )
               { printf ("(-scale sx)"); print_usage( usage ); }
            if ( sscanf( argv[argn+2], "%f", &data_scale_y) != 1 )
               { printf ("(-scale sy)"); print_usage( usage ); }
            if ( sscanf( argv[argn+3], "%f", &data_scale_z) != 1 )
               { printf ("(-scale sz)"); print_usage( usage ); }
            argn += 2;
            scale_flag = TRUE;
         } else

 /* data rotation center: */
         if ( strncmp(argv[argn],"-center",MAX(strlen(argv[argn]),3)) == 0 ) {
            if ( sscanf( argv[argn+1], "%f", &data_center_x) != 1 )
               { printf ("(-center cx)"); print_usage( usage ); }
            if ( sscanf( argv[argn+2], "%f", &data_center_y) != 1 )
               { printf ("(-center cy)"); print_usage( usage ); }
            if ( sscanf( argv[argn+3], "%f", &data_center_z) != 1 )
               { printf ("(-center cz)"); print_usage( usage ); }
            argn += 2;
            center_flag = TRUE;
         } else

 /* data translations: */
         if ( strncmp(argv[argn],"-trans",MAX(strlen(argv[argn]),3)) == 0 ) {
            if ( sscanf( argv[argn+1], "%f", &data_trans_x) != 1 )
               { printf ("(-trans tx)"); print_usage( usage ); }
            if ( sscanf( argv[argn+2], "%f", &data_trans_y) != 1 )
               { printf ("(-trans ty)"); print_usage( usage ); }
            if ( sscanf( argv[argn+3], "%f", &data_trans_z) != 1 )
               { printf ("(-trans tz)"); print_usage( usage ); }
            argn += 2;
            trans_flag = TRUE;
         } else


/* tolerence for fitting */
         if ( strncmp(argv[argn],"-tol",MAX(strlen(argv[argn]),3)) == 0 ) {
            if ( sscanf( argv[argn+1], "%f", &ftol ) != 1 )
               { printf ("(-tol)"); print_usage( usage ); }
         } else
/* simplex size for fitting */
         if ( strncmp(argv[argn],"-simplex",MAX(strlen(argv[argn]),3)) == 0 ) {
            if ( sscanf( argv[argn+1], "%f", &simplex_size ) != 1 )
               { printf ("(-simplex)"); print_usage( usage ); }
         } else

/* kernel size */
         if ( strncmp(argv[argn],"-k<ernelsize",MAX(strlen(argv[argn]),2)) == 0 ) {
            if ( sscanf( argv[argn+1], "%f", &kernel_size ) != 1 )
               { printf ("(-kernel)"); print_usage( usage ); }
         } else


            { 
              printf ("(last else argn=%d,arc=%d)\n",argn,argc); 
              printf ("(%s\n)\n",argv[argn]); 
              print_usage( usage ); 
            }

         argn += 2;

