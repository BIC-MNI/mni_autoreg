#include <internal_volume_io.h>

static char *default_dim_names[N_DIMENSIONS] =
   { MIzspace, MIyspace, MIxspace };


BOOLEAN vol_cog(Volume d1, Volume m1, float *centroid)
{

  Real
    sx,sy,sz,si,
    true_value,
    tx,ty,tz;
  int
    i,j,k,
    count[MAX_DIMENSIONS];
    

  sx = 0.0;			/* init centroid vars */
  sy = 0.0;
  sz = 0.0;
  si = 0.0;

  get_volume_sizes(d1, count);
				/* loop over all voxels */
  for_less(i,0,count[0])
    for_less(j,0,count[1])
      for_less(k,0,count[2]) {
	
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



char *prog_name;

main(int argc, char *argv[])
{
  Volume vol,mask;
  float cog[3];
  double step[3];

  prog_name = argv[0];

  vol = mask = NULL;

  if (argc<2 || argc>3) {
    (void)printf("usage: volume_cog volume.mnc [mask.mnc]\n");
    exit(EXIT_FAILURE);
  }

  input_volume(argv[1], 3, (char **)NULL, NC_UNSPECIFIED, FALSE, 0.0,0.0,
	       TRUE, &vol, (minc_input_options *)NULL);

  if (argc==3) {
    input_volume(argv[2], 3, (char **)NULL, NC_UNSPECIFIED, FALSE, 0.0,0.0,
		 TRUE, &mask, (minc_input_options *)NULL);    
  }

  if (vol_cog(vol, mask, cog )) {
    (void)print ("%f %f %f\n",cog[0],cog[1],cog[2]);
    exit(EXIT_SUCCESS);
    
  }
  else {
    (void)print ("No COG can be calculated\n");
    exit(EXIT_FAILURE);
  }

}
