
#include <volume_io.h>

static char *default_dim_names[VIO_N_DIMENSIONS] =
   { MIzspace, MIyspace, MIxspace };


char *prog_name;


int main(int argc, char *argv[])
{
  int 
    i,j,k, p,r,s, flag,
    sizes1[3],sizes2[3], sizes3[3];
  VIO_Real
    corr, v1, v2, v3, s1, s2, s12,
    u,v,w,  x,y,z;
  VIO_Status status;

  VIO_Volume 
    data1, data2, mask;
  char *f1, *f2, *mf;

  if (argc<3) {
    print ("usage:  xcorr_vol vol1.mnc vol2.mnc [mask.mnc]\n");
    exit(EXIT_FAILURE);
  }

  prog_name = argv[0];
  f1        = argv[1];
  f2        = argv[2];  
  if (argc==4)
    mf      = argv[3];

  status = input_volume(f1, 3, default_dim_names, NC_UNSPECIFIED, FALSE, 0.0,0.0,
                        TRUE, &data1, (minc_input_options *)NULL); 
  if (status!=OK) {
    print ("Error reading %s.\n",f1);
    exit(EXIT_FAILURE);
  }
  
  status = input_volume(f2, 3, default_dim_names, NC_UNSPECIFIED, FALSE, 0.0,0.0,
                        TRUE, &data2, (minc_input_options *)NULL); 
 
  if (status!=OK) {
    print ("Error reading %s.\n",f2);
    exit(EXIT_FAILURE);
  }

  if (argc==4) {
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
  get_volume_sizes(data2,sizes2);
  
  if (sizes1[0] != sizes2[0] || sizes1[1] != sizes2[1] || sizes1[2] != sizes2[2]) {
    print ("Size mismatch between %s and %s.\n",f1,f2);
    print ("%d,%d,%d != %d,%d,%d\n", 
           sizes1[0],sizes1[1],sizes1[2],
           sizes2[0],sizes2[1],sizes2[2]);
    exit(EXIT_FAILURE);
  }
  
  s1 = s2 = s12 = 0.0;
  for(i=0; i<sizes1[0]; i++) {
    for(j=0; j<sizes1[1]; j++) {
      for(k=0; k<sizes1[2]; k++) {

        GET_VALUE_3D( v1 ,  data1, i, j, k);
        GET_VALUE_3D( v2 ,  data2, i, j, k);

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
          s1  += v1*v1;
          s2  += v2*v2;
          s12 += v1*v2;
        }
      }
    }
  }

  if (s1==0.0 || s12==0.0 || s2==0.0) {
    print ("Null value somewhere: s1=%f, s2=%f, s12=%f.\n",s1,s2,s12);
    exit(EXIT_FAILURE);
  }
  else {
    corr = s12 / (sqrt(s1)*sqrt(s2));
    print ("%15.12f\n",corr);
    exit(EXIT_SUCCESS);
  }

}
