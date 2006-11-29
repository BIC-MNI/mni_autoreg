#include <volume_io.h>
#include <config.h>

static char *default_dim_names[VIO_N_DIMENSIONS] =
   { MIzspace, MIyspace, MIxspace };


VIO_BOOL vol_cog(VIO_Volume d1, VIO_Volume m1, VIO_Real *centroid)
{

  VIO_Real
    sx,sy,sz,si,
    true_value,
    tx,ty,tz;
  int
    i,j,k,
    count[VIO_MAX_DIMENSIONS];
    

  sx = 0.0;                        /* init centroid vars */
  sy = 0.0;
  sz = 0.0;
  si = 0.0;

  get_volume_sizes(d1, count);
                                /* loop over all voxels */
  for(i=0; i<count[0]; i++)
    for(j=0; j<count[1]; j++)
      for(k=0; k<count[2]; k++) {
        
        convert_3D_voxel_to_world(d1, i,j,k, &tx, &ty, &tz);
        
              if (point_not_masked(m1, tx, ty, tz)) {
          
          GET_VALUE_3D( true_value, d1, i,j,k);
          
          sx +=  tx * true_value;
          sy +=  ty * true_value;
          sz +=  tz * true_value;
          si += true_value;
          
        } 
        
      }
                                /* calc centroids */
  if (si!=0.0) {
    centroid[0] = sx/ si;
    centroid[1] = sy/ si;
    centroid[2] = sz/ si;
    return(TRUE);
  }
  else {
    return(FALSE);
  }
}

VIO_BOOL vol_stat(VIO_Volume d1, VIO_Volume m1, 
                 VIO_Real *mean, 
                 VIO_Real *std, 
                 VIO_Real *minimum, 
                 VIO_Real *maximum)
{

  VIO_Real
      tx,ty,tz,
      true_value,
      sum, sum2,
      min,max;
  
  int
    i,j,k, counter,
    count[VIO_MAX_DIMENSIONS];
    

  sum = 0.0;
  sum2 = 0.0;
  counter = 0;
  min = DBL_MAX;
  max = -DBL_MAX;
  
  
  get_volume_sizes(d1, count);
                                /* loop over all voxels, tally stats */
  for(i=0; i<count[0]; i++)
    for(j=0; j<count[1]; j++)
      for(k=0; k<count[2]; k++) {
        
        convert_3D_voxel_to_world(d1, i,j,k, &tx, &ty, &tz);
        
              if (point_not_masked(m1, tx, ty, tz)) {
          
          GET_VALUE_3D( true_value, d1, i,j,k);
          
          sum  += true_value;
          sum2 += true_value*true_value;
          
          if (true_value < min) 
              min = true_value;
          if (true_value > max)
              max = true_value;

          counter++;
          
        } 
        
      }
                                /* compute stats */
  if (counter>0) {

      *mean = (VIO_Real) (sum / counter);
      if (counter>1) 
      {
          *std  = (VIO_Real) sqrt((sum2*counter - sum*sum)/(counter * (counter-1)));
      }
      else
          *std = 0.0;      

      *minimum = min;
      *maximum = max;      
      
      return(TRUE);
  }
  else {
    return(FALSE);
  }
}



char *prog_name;

main(int argc, char *argv[])
{
  VIO_Volume vol,mask;
  VIO_Real mean,std,min,max,cog[3];
  double step[3];

  prog_name = argv[0];

  vol = mask = NULL;

  if (argc<2 || argc>3) {
    (void)printf("usage: volume_stat volume.mnc [mask.mnc]\n");
    exit(EXIT_FAILURE);
  }

  input_volume(argv[1], 3, (char **)NULL, NC_UNSPECIFIED, FALSE, 0.0,0.0,
               TRUE, &vol, (minc_input_options *)NULL);

  if (argc==3) {
    input_volume(argv[2], 3, (char **)NULL, NC_UNSPECIFIED, FALSE, 0.0,0.0,
                 TRUE, &mask, (minc_input_options *)NULL);    
  }

  if (vol_stat(vol, mask, &mean, &std, &min, &max )) {
    (void)print ("mean.std.min.max %f %f %f %f\n",mean,std,min,max);
    exit(EXIT_SUCCESS);
    
  }
  else {
    (void)print ("No stats can be calculated\n");
    exit(EXIT_FAILURE);
  }

/*
  if (vol_cog(vol, mask, cog )) {
    (void)print ("%f %f %f\n",cog[0],cog[1],cog[2]);
    exit(EXIT_SUCCESS);
    
  }
  else {
    (void)print ("No COG can be calculated\n");
    exit(EXIT_FAILURE);
  }
*/

}
