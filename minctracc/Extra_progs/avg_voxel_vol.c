#include <volume_io.h>

static char *default_dim_names[VIO_N_DIMENSIONS] =
   { MIzspace, MIyspace, MIxspace };


char *prog_name;


main(int argc, char *argv[])
{
  int 
    histo[256],
    hist_flag,
    index[VIO_MAX_DIMENSIONS],
    i,j,k, p,r,s, flag,
    sizes1[3],sizes2[3], sizes3[3];
  VIO_Real
    r_max, r_min,
    voxel[VIO_MAX_DIMENSIONS],
    std,mean_value, v1, v2, v3, s1, s2, smag, s12,
    u,v,w,  x,y,z;
  VIO_Status status;

  VIO_Volume 
    data1, data2, mask;
  char *f1, *f2, *mf;
  long int count;

  if (argc<2) {
    print ("usage:  avg_voxel_val vol1.mnc  [mask.mnc] [-hist]\n");
    exit(EXIT_FAILURE);
  }

  prog_name = argv[0];
  f1        = argv[1];
  hist_flag = FALSE;


  mf = NULL;

  if (argc>2) {
    
    if ( strcmp( argv[argc-1], "-hist")==0 ) {
      hist_flag = TRUE;
      argc--;
    }
    if (argc>2)
      mf      = argv[2];
  }

                                /* read in data file */
  status = input_volume(f1, 3, default_dim_names, NC_UNSPECIFIED, FALSE, 0.0,0.0,
                        TRUE, &data1, (minc_input_options *)NULL); 
  if (status!=OK) {
    print ("Error reading %s.\n",f1);
    exit(EXIT_FAILURE);
  }
  
                                /* read in optional mask file */
  if (mf != NULL) {
    print("Reading mask %s\n",mf);
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
  get_volume_minimum_maximum_real_value(data1, &r_min, &r_max);
  
  smag = s1 = s2 = s12 = 0.0; count = 0;
  for(i=0; i<MAX_DIMENSIONS; i++) index[i] = 0;

  for(i=0; i<256; i++) histo[i] = 0;

  for(index[X]=0; index[X]<sizes1[X]; index[X]++)
    for(index[Y]=0; index[Y]<sizes1[Y]; index[Y]++)
      for(index[Z]=0; index[Z]<sizes1[Z]; index[Z]++) {
        
        v1 = get_volume_real_value(data1,
                                   index[0],index[1],index[2],index[3],index[4]);
        
        flag = TRUE;
        
        if (mask != NULL) {

          for(i=0; i<MAX_DIMENSIONS; i++) voxel[i] = (VIO_Real)index[i];

          convert_voxel_to_world(data1, voxel, &x, &y, &z);
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
          smag += fabs(v1);
          count++;
        }
      }

  
  if (count==0) {
    print ("No voxels!");
  }
  else {

    if (!hist_flag) {
      mean_value = s1 / count;
      print ("mean value = %f\n", mean_value);
      
      mean_value = smag / count;
      print ("mean mag   = %f\n", mean_value);
      
      v1 = (double)( (count*s2 - s1*s1)/(1.0*count * (count-1)));
      std = sqrt( v1 );
      
      print ("std = %f\n",std);
    }
      
      
    mean_value = s1 / count;
    v1 = (double)( (count*s2 - s1*s1)/(1.0*count * (count-1)));
    std = sqrt( v1 );
    
/*
    r_min = mean_value - 5*std;
    r_max = mean_value + 5*std;
*/
      
    print ("%f %f\n",r_min, r_max);
    if (hist_flag) {
      for(index[X]=0; index[X]<sizes1[X]; index[X]++)
        for(index[Y]=0; index[Y]<sizes1[Y]; index[Y]++)
          for(index[Z]=0; index[Z]<sizes1[Z]; index[Z]++) {
            
            v1 = get_volume_real_value(data1,
                                       index[0],index[1],index[2],index[3],index[4]);
            
            flag = TRUE;
            
            if (mask != NULL) {
              
              for(i=0; i<MAX_DIMENSIONS; i++) voxel[i] = (VIO_Real)index[i];
              
              convert_voxel_to_world(data1, voxel, &x, &y, &z);
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
              smag += fabs(v1);
              count++;
              
              if (hist_flag) {
                i = 255 * (v1 - r_min) / (r_max - r_min);
                if (i<0) i = 0;
                else
                  if (i>255) i = 255;
                histo[i]++;
              }
            }
          }
      
      for(i=0; i<256; i++) {
        print ("%d\n",histo[i]);
      }
      
    }
  }
  

}
