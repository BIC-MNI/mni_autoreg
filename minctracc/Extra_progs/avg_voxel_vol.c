
#include <volume_io.h>

static char *default_dim_names[N_DIMENSIONS] =
   { MIzspace, MIyspace, MIxspace };


char *prog_name;


main(int argc, char *argv[])
{
  int 
    i,j,k, p,r,s, flag,
    sizes1[3],sizes2[3], sizes3[3];
  Real
    std,corr, v1, v2, v3, s1, s2, s12,
    u,v,w,  x,y,z;
  Status status;

  Volume 
    data1, data2, mask;
  char *f1, *f2, *mf;
  long int count;

  if (argc<3) {
    print ("usage:  xcorr_vol.c vol1.mnc  [mask.mnc]\n");
    exit(EXIT_FAILURE);
  }

  prog_name = argv[0];
  f1        = argv[1];
  if (argc==3)
    mf      = argv[2];

  status = input_volume(f1, 3, default_dim_names, NC_UNSPECIFIED, FALSE, 0.0,0.0,
			TRUE, &data1, (minc_input_options *)NULL); 
  if (status!=OK) {
    print ("Error reading %s.\n",f1);
    exit(EXIT_FAILURE);
  }
  

  if (argc==3) {
    status = input_volume(mf, 3, default_dim_names, NC_UNSPECIFIED, FALSE, 0.0,0.0,
			  TRUE, &mask, (minc_input_options *)NULL); 
    
    if (status!=OK) {
      print ("Error reading %s.\n",mf);
      exit(EXIT_FAILURE);
    }
    get_volume_sizes(mask,sizes3); 
  }
  else 
    mask = (Volume)NULL;
  
  get_volume_sizes(data1,sizes1);

  
  s1 = s2 = s12 = 0.0; count = 0;
  for_less(i,0,sizes1[0]) {
    for_less(j,0,sizes1[1]) {
      for_less(k,0,sizes1[2]) {

	GET_VALUE_3D( v1 ,  data1, i, j, k);

	flag = TRUE;

	if (mask != NULL) {
	  convert_3D_voxel_to_world(data1, i,j,k, &x, &y, &z);
	  convert_3D_world_to_voxel(mask,  x,y,z, &u, &v, &w);
	  p = ROUND(u);
	  r = ROUND(v);
	  s = ROUND(w);
	  if (p>=0 &&  p<sizes3[0] &&
	      r>=0 &&  r<sizes3[1] &&
	      s>=0 &&  s<sizes3[2]) {
	    GET_VALUE_3D( v3 ,  mask, p, r, s);
	    if (v3 < 0.5)
	      flag = FALSE;
	  }
	  else 
	    flag = FALSE;
	}

	if (flag) {
	  s1  += v1;
	  s2  += v1*v1;
	  count++;
	}
      }
    }
  }


  
  if (count==0) {
    print ("No voxels!");
  }
  else {
    corr = s1 / count;
    std = sqrt( (count*s2 - s1*s1)/(count * (count-1)) );
    print ("%f %f\n",corr,std);
  }

}
