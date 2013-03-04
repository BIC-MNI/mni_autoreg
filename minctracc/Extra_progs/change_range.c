
#include <volume_io.h>

static char *default_dim_names[VIO_N_DIMENSIONS] =
   { MIzspace, MIyspace, MIxspace };


char *prog_name;


main(int argc, char *argv[])
{
  int 
    i,j,k,
    sizes[3];
  VIO_Real
    min,max, v1, v2;
  VIO_Status status;

  VIO_Volume 
    data1, data2;
  char *f1, *f2;

  if (argc<4) {
    print ("usage:  change_range in.mnc out.mnc min_val max_val\n");
    exit(EXIT_FAILURE);
  }

  prog_name = argv[0];
  f1        = argv[1];
  f2        = argv[2];  
  min = atof(argv[3]);
  max = atof(argv[4]);

  status = input_volume(f1, 3, default_dim_names, NC_UNSPECIFIED, FALSE, 0.0,0.0,
                        TRUE, &data1, (minc_input_options *)NULL); 
  if (status!=OK) {
    print ("Error reading %s.\n",f1);
    exit(EXIT_FAILURE);
  }
  

  data2 = copy_volume_definition(data1, NC_UNSPECIFIED, FALSE, 0.0, 0.0);
  data2->data = data1->data;

  set_volume_real_range(data2, min, max);

  get_volume_sizes(data1,sizes);

  for(i=0; i<sizes[0]; i++) {
    for(j=0; j<sizes[1]; j++) {
      for(k=0; k<sizes[2]; k++) {

        GET_VALUE_3D( v1 ,  data1, i, j, k);
        if (v1>max) v1 = max;
        if (v1<min) v1 = min;
        v2 = CONVERT_VALUE_TO_VOXEL(data2, v1);
        SET_VOXEL_3D( data2, i, j, k, v2);

      }
    }
  }

  status = output_modified_volume(f2, NC_UNSPECIFIED, FALSE,
                                  min,max, data2, f1, NULL, NULL);

}
