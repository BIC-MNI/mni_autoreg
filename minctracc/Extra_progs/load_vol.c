
#include <internal_volume_io.h>

static char *default_dim_names[N_DIMENSIONS] =
   { MIzspace, MIyspace, MIxspace };


char *prog_name;


main(int argc, char *argv[])
{
  int 
    i,j,k,
    sizes[3];
  Real
     start[3],pt[3];
  Status status;

  Volume 
    data1;


  status = input_volume(argv[1], 3, default_dim_names, 
                        NC_UNSPECIFIED, FALSE, 0.0,0.0,
			TRUE, &data1, (minc_input_options *)NULL); 


  if (status!=OK) {
    print ("Error reading %s.\n",argv[1]);
    exit(EXIT_FAILURE);
  }
  
  start[0] = start[1] = start[2] = 0.0;
  convert_voxel_to_world(data1, start,
                         &pt[0], &pt[1], &pt[2]);

  print ("start: %9.3f %9.3f %9.3f\n", pt[0], pt[1], pt[2]);
  
}
