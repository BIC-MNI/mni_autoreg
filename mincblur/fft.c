/* ----------------------------- MNI Header -----------------------------------
@NAME       : mincdxyz2.c
@INPUT      : the base name for files containing dx,dy,dz,dxyz and blur volumes.
@OUTPUT     : a volume whose values represent the second derivative of the
              blurred intensity volume along the gradient direction (i.e. 
	      perpendicular to the gradient magnitude iso-surface.
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Wed Sep 22 13:31:00 EST 1993
@MODIFIED   : $Log: fft.c,v $
@MODIFIED   : Revision 1.1  1993-10-15 13:52:44  louis
@MODIFIED   : Initial revision
@MODIFIED   :
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/mincblur/fft.c,v 1.1 1993-10-15 13:52:44 louis Exp $";
#endif

#include <mincdxyz2.h>
#include <limits.h>

main(int argc, char *argv[])
{

  Volume
    data_dx, data_dy, data_dz, data_blur, data_dxyz,
    grad2;

  STRING
    file_blur,
    file_dx,
    file_dy,
    file_dz,
    file_dxyz;

  char
    *infile_base,
    *outfilename,
    *history;
  int
    i,j,k,
    size[3];
  Real
    value,
    min, max,
    thick[3];
  Status
    status;
  
  General_transform *transform, copied_transform;

				/* set up global variables. */
  debug = FALSE;
  verbose = 1;
  clobber_flag = FALSE;
  fwhm = standard = 0.0;

  history = time_stamp(argc, argv);
  prog_name = argv[0];


				/* Call ParseArgv to interpret all command line args */

  if (ParseArgv(&argc, argv, argTable, 0) || (argc!=4)) {
    (void) fprintf(stderr, 
		   "\nUsage: %s [<options>] <input_basename> <fwhm> <outputfile> \n", 
		   prog_name);
    (void) fprintf(stderr,"       %s [-help]\n\n", prog_name);
    exit(EXIT_FAILURE);
  }

  sscanf(argv[2],"%lf",&fwhm);


  if ( fwhm<=0.0) {
    print_error ("Must specify positive fwhm command line.", 
		 __FILE__, __LINE__);
  }


  infile_base = argv[1];	/* set up necessary file names */
  outfilename = argv[3]; 

  file_blur[0] = '\0';
  strcpy(file_blur,infile_base);
  strcat(file_blur,"_blur.mnc\0");

  file_dx[0] = '\0';
  strcpy(file_dx,infile_base);
  strcat(file_dx,"_dx.mnc\0");

  file_dy[0] = '\0';
  strcpy(file_dy,infile_base);
  strcat(file_dy,"_dy.mnc\0");

  file_dz[0] = '\0';
  strcpy(file_dz,infile_base);
  strcat(file_dz,"_dz.mnc\0");

  file_dxyz[0] = '\0';
  strcpy(file_dxyz,infile_base);
  strcat(file_dxyz,"_dxyz.mnc\0");

  print ("file_blur : %s\n", file_blur);
  print ("file_dx : %s\n", file_dx);
  print ("file_dy : %s\n", file_dy);
  print ("file_dz : %s\n", file_dz);
  print ("file_dxyz : %s\n", file_dxyz);

  if ( !file_exists(file_blur) ) print_error("File %s not found\n", __FILE__, __LINE__,file_blur);
  if ( !file_exists(file_dx) )   print_error("File %s not found\n", __FILE__, __LINE__,file_dx);
  if ( !file_exists(file_dy) )   print_error("File %s not found\n", __FILE__, __LINE__,file_dy);
  if ( !file_exists(file_dz) )   print_error("File %s not found\n", __FILE__, __LINE__,file_dz);
  if ( !file_exists(file_dxyz) ) print_error("File %s not found\n", __FILE__, __LINE__,file_dxyz);

  if (!clobber_flag && file_exists(outfilename)) {
    print ("File %s exists.\n",outfilename);
    print ("Use -clobber to overwrite.\n");
    return ERROR;
  }

  status = input_volume(file_blur, 3, default_dim_names, NC_UNSPECIFIED,
			FALSE, 0.0, 0.0, TRUE, &data_blur,NULL );

  print ("file_blur : %s\n", file_blur);
  print ("file_dx : %s\n", file_dx);
  print ("file_dy : %s\n", file_dy);
  print ("file_dz : %s\n", file_dz);
  print ("file_dxyz : %s\n", file_dxyz);

  if (status==OK) status = input_volume(file_dx, 3,  default_dim_names, NC_UNSPECIFIED, 
					FALSE, 0.0, 0.0, TRUE, &data_dx, NULL );
  if (status==OK) status = input_volume(file_dy, 3,  default_dim_names, NC_UNSPECIFIED, 
					FALSE, 0.0, 0.0, TRUE, &data_dy, NULL );
  if (status==OK) status = input_volume(file_dz, 3,  default_dim_names, NC_UNSPECIFIED, 
					FALSE, 0.0, 0.0, TRUE, &data_dz, NULL );
  if (status==OK) status = input_volume(file_dxyz,3, default_dim_names, NC_UNSPECIFIED, 
					FALSE, 0.0, 0.0, TRUE, &data_dxyz , NULL);

  if (status != OK)
    print_error ("Error loading data volumes.\n", __FILE__, __LINE__);


  grad2 = copy_volume_definition( data_blur, NC_FLOAT, FALSE, -FLT_MAX, FLT_MAX);

  status = build_2nd_derivative_volume(data_blur, 
				       data_dx, data_dy, data_dz, data_dxyz,
				       grad2, fwhm);

  /* find min and max in volume before save */
  
  get_volume_sizes(grad2, size);
  min = FLT_MAX;
  max = -FLT_MAX;
  for_less(i,0,size[0])
    for_less(j,0,size[1])
      for_less(k,0,size[2]) {
	GET_VOXEL_3D(value, grad2, i,j,k);
	value = CONVERT_VOXEL_TO_VALUE(grad2, value);
	if (value<min)
	  min=value;
	if (value>max)
	  max=value;
      }

print ("min and max: %f %f\n",min, max);

  set_volume_voxel_range(grad2, min, max);
  set_volume_real_range(grad2,  min, max);


  if (status==OK)
    status = output_modified_volume(outfilename, NC_UNSPECIFIED, FALSE, 0.0, 0.0,
				    grad2, file_blur, history, NULL);
  else
    print_error ("Error building 2nd derivative volume.\n", __FILE__, __LINE__);

  if (status!=OK)
    print_error ("Error saving 2nd derivative volume.\n", __FILE__, __LINE__);
    

  delete_volume(data_blur);
  delete_volume(data_dx);
  delete_volume(data_dy);
  delete_volume(data_dz);
  delete_volume(data_dxyz);
  delete_volume(grad2);

  return(status);

}
