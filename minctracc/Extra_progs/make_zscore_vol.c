
#include <volume_io.h>
#include <config.h>

static char *default_dim_names[VIO_N_DIMENSIONS] =
   { MIzspace, MIyspace, MIxspace };


char *prog_name;

void make_zscore_volume(VIO_Volume d1, VIO_Volume m1, 
                               VIO_Real *threshold);


main(int argc, char *argv[])
{
  int 
    count,
    i,j,k, p,r,s, flag,
    sizes1[3],sizes2[3], sizes3[3];
  VIO_Real
    thresh, mean, std;
  VIO_Real
    min1,max1,min2,max2, mean1, std1, mean2, std2,
    corr, v1, v2, v3, s1, s2, s12,
    u,v,w,  x,y,z;
  VIO_Status status;

  VIO_Volume 
    data1, data2, mask;
  char *f1, *f2, *mf;

  if (argc<4) {
    print ("usage:  xcorr_vol.c invol.mnc outvol.mnc threshold [mask.mnc]\n");
    exit(EXIT_FAILURE);
  }

  prog_name = argv[0];
  f1        = argv[1];
  f2        = argv[2];
  thresh    = atof(argv[3]);
  if (argc==5)
    mf      = argv[4];

  status = input_volume(f1, 3, default_dim_names, NC_UNSPECIFIED, FALSE, 0.0,0.0,
                        TRUE, &data1, (minc_input_options *)NULL); 
  if (status!=OK) {
    print ("Error reading %s.\n",f1);
    exit(EXIT_FAILURE);
  }
  
  if (argc==5) {
    status = input_volume(mf, 3, default_dim_names, NC_UNSPECIFIED, FALSE, 0.0,0.0,
                          TRUE, &mask, (minc_input_options *)NULL); 
    
    if (status!=OK) {
      print ("Error reading %s.\n",mf);
      exit(EXIT_FAILURE);
    }
    get_volume_sizes(mask,sizes3); 
  }
  else 
    mask = (VIO_Volume)NULL;
  
  get_volume_sizes(data1,sizes1);
  
  make_zscore_volume(data1, mask, &thresh);

  status = output_modified_volume(f2, NC_UNSPECIFIED, TRUE, 0.0, 0.0,
                                  data1, f1, (char *)NULL,
                                  (minc_output_options *)NULL);

  if (status==OK)
    exit(EXIT_SUCCESS);
  else {
      print ("Error saving %s.\n",f2);
      exit(EXIT_FAILURE);
  }

}


