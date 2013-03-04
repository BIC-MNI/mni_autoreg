#include <volume_io.h>
#include <config.h>
#include <Proglib.h>


static char *default_dim_names[VIO_N_DIMENSIONS] =
   { MIzspace, MIyspace, MIxspace };

#define LONGISH_CONST 8
#define BASELIM_CONST 3
#define LOWER_VOLUME_LIMIT 15000
#define UPPER_VOLUME_LIMIT 40000

int    clobber,debug;
VIO_Real   base_limit,
       length_limit, 
       vol_thres[2];

char   mask_name[1024];

  static ArgvInfo argTable[] = {
  {"-length_limit",      ARGV_FLOAT, (char *) 0, (char *) &length_limit,
       "Length ratio limit for 'tubeness' (long/short axis)."},
  {"-base_limit",      ARGV_FLOAT, (char *) 0, (char *) &length_limit,
       "Base ratio limit for 'tubeness' (mid/short axis)."},
  {"-volume_limits",      ARGV_FLOAT, (char *) 2, (char *)vol_thres,
       "Threshold limits (lower upper) for 'tubeness'."},
  {"-mask_name", ARGV_STRING, (char *) 0, (char *) mask_name,
       "Name of (optional) output mask volume."},
  {"-debug",     ARGV_CONSTANT, (char *) TRUE, (char *) &debug,
       "Write out debug info (default = no debug)."},
  {"-clobber",     ARGV_CONSTANT, (char *) TRUE, (char *) &clobber,
       "Overwrite existing file (default = no clobber)."},
  {"-version", ARGV_FUNC, (char *) print_version_info, (char *)MNI_AUTOREG_LONG_VERSION,
       "Print out version info and exit."},
  {NULL, ARGV_END, NULL, NULL, NULL}
  };
  
void remove_voxels(VIO_Volume d1, int number)
{
    int 
        v, sizes[3], i,j,k;
    VIO_Real 
        zero, voxel_value, obj_value;

    get_volume_sizes(d1, sizes);

    obj_value = CONVERT_VALUE_TO_VOXEL(d1, number);
    obj_value = FLOOR( obj_value );
    zero = CONVERT_VALUE_TO_VOXEL(d1, 0.0);
    zero = CONVERT_VOXEL_TO_VALUE(d1, zero);
    
    for(i=0; i<sizes[0]; i++)
        for(j=0; j<sizes[1]; j++)
            for(k=0; k<sizes[2]; k++) 
            {
                GET_VOXEL_3D(voxel_value, d1, i,j,k);
                if (voxel_value == obj_value) 
                    SET_VOXEL_3D(d1, i,j,k, zero);
            }
    
    
}

VIO_BOOL vol_to_cov(VIO_Volume d1, int number, float *vol, float *centroid, float **covar)
{

  VIO_Real
      obj_value,
    tx,ty,tz,
    voxel_value;
  int
    i,j,k,r,c,s,
    limits[3];

  float
    t,
    sxx,syy,szz,
    sxy,syz,sxz,
    sx,sy,sz,si; 
  VIO_Real
    thickness[3];
  int
    sizes[3];

  VIO_Real
      dx,dy,dz,
      zero,
      true_value,
      position[3];


  zero = CONVERT_VALUE_TO_VOXEL(d1, 0.0);
  zero = CONVERT_VOXEL_TO_VALUE(d1, zero);

  obj_value = CONVERT_VALUE_TO_VOXEL(d1, number);
  obj_value = FLOOR( obj_value );
  

  get_volume_separations(d1, thickness);
  get_volume_sizes(d1, sizes);


                                /* calculate centroids first */
  sx = 0.0;
  sy = 0.0;
  sz = 0.0;
  si = 0.0;
  
  for(i=0; i<sizes[0]; i++)
      for(j=0; j<sizes[1]; j++)
          for(k=0; k<sizes[2]; k++) 
          {
              GET_VOXEL_3D(voxel_value, d1, i,j,k);
              if (voxel_value == obj_value) {
                  
                  convert_3D_voxel_to_world(d1, (VIO_Real)i, (VIO_Real)j, (VIO_Real)k, 
                                            &tx, &ty, &tz); 
                  sx +=  tx;
                  sy +=  ty;
                  sz +=  tz;
                  si += 1.0;
                  
              }
          }
  
  *vol = si * thickness[0]*thickness[1]*thickness[2];
  

  if (si!=0.0) {
      sx = sx / si;
      sy = sy / si;
      sz = sz / si;
      
      centroid[1] = sx ;
      centroid[2] = sy ;
      centroid[3] = sz ;
      
      sxx = syy = szz = 0.0;
      sxy = syz = sxz = 0.0;
      
      /* now calculate variances and co-variances */
      for(i=0; i<sizes[0]; i++) {          
          for(j=0; j<sizes[1]; j++) {               
              for(k=0; k<sizes[2]; k++)   {
                  GET_VOXEL_3D(voxel_value, d1, i,j,k);
                  if (voxel_value == obj_value) {
                      
                      convert_3D_voxel_to_world(d1, (VIO_Real)i, (VIO_Real)j, (VIO_Real)k, 
                                                &tx, &ty, &tz); 
                      
                      dx = tx-sx;
                      dy = ty-sy;
                      dz = tz-sz;

                      
                      sxx += dx * dx;
                      syy += dy * dy;
                      szz += dz * dz;
                      sxy += dx * dy;
                      syz += dy * dz;
                      sxz += dx * dz;
                  }
                  
                  
              }
          }
      }
      
      
      covar[1][1] = sxx/si; covar[1][2] = sxy/si; covar[1][3] = sxz/si;
      covar[2][1] = sxy/si; covar[2][2] = syy/si; covar[2][3] = syz/si;
      covar[3][1] = sxz/si; covar[3][2] = syz/si; covar[3][3] = szz/si;
      
      return(TRUE);
      
  }
  else {
      if (debug) print ("no vol for object: %f %d\n",obj_value, number);
      
      return(FALSE);
  }
}

/* use the volume and the relative lengths of the principal axis to
   decide whether the object may be a tube. Recall that the principal
   axis are stored in columns 1,2 and 3 of the praxes[4][4] matrix. */

VIO_BOOL this_is_a_tube( float vol, float *centroid, float **praxes) 
{
    VIO_Real 
        max_len, min_len, lengths[3], tmp_len, base;
    VIO_BOOL
        base_val,
        longish,
        right_vol;
    int i, tmp_i;
    

                                /* check volume */
    right_vol =  (vol > vol_thres[0] && vol < vol_thres[1]);


                                /* check longishness  */
    for(i=0; i<3; i++) 
    {
        lengths[i] = sqrt(praxes[1][i+1]*praxes[1][i+1] +
                          praxes[2][i+1]*praxes[2][i+1] +
                          praxes[3][i+1]*praxes[3][i+1]);
        
    }

    base =  lengths[1] * lengths[2] ;
    
    if (base > 0.0) 
    {
         if (debug) print ("vol = %f lengths = %f,%f,%f; base = %f, len ratio = %f \n", 
                           vol, lengths[0] , lengths[1], lengths[2], 
                           lengths[1]  / lengths[2], lengths[0]  / lengths[2]);
        
        longish =  (lengths[0]  / lengths[2] ) > length_limit;
        base_val = (lengths[1]  / lengths[2] ) < base_limit;
    }
    else 
    {
        return(FALSE);
    }
    
    return( base_val && longish && right_vol );    
}

float angle_from( float x, float y) 
{
    float angle;
    
    angle = 0;

    if (y!=0.0) 
    {
        angle = 180 * atan( x/ y) / M_PI;        
    }
    else
        angle = 0.0;
    

    return(angle);    
}

void write_out_tube_options (float *centroid, float **praxes) 
{
    VIO_Real 
        xrot, yrot, zrot, max_len, min_len, lengths[3];
    VIO_BOOL
        longish,
        right_vol;
    int j,i,max_ind;
    float  *ang;
    

    ALLOC  (ang,4);
    
                                /* get principal axis lengths */
    for(i=0; i<3; i++) 
    {
        lengths[i] = sqrt(praxes[1][i+1]*praxes[1][i+1] +
                          praxes[2][i+1]*praxes[2][i+1] +
                          praxes[3][i+1]*praxes[3][i+1]);
    }

                        /* find index of principal axis */
    if ( praxes[1][1] > praxes[2][1])
        max_ind = 1; 
    else
        max_ind = 2;
    if (praxes[max_ind][1] < praxes[3][1])
        max_ind = 3;
    
    xrot = yrot = zrot = 0.0;
    
    switch (max_ind) 
    {
    case 1:                        /* x axis */
        yrot = angle_from( -praxes[3][1], praxes[1][1]);
        zrot = angle_from( praxes[2][1], praxes[1][1]);
        print ("model= LR_tube.mnc ");
        break;
    case 2:                        /* y axis */
        xrot = angle_from( praxes[3][1], praxes[2][1]);
        zrot = angle_from( -praxes[1][1], praxes[2][1]);
        print ("model= AP_tube.mnc ");
        break;
    case 3:                        /* z axis */
        xrot = angle_from( -praxes[2][1], praxes[3][1]);
        yrot = angle_from( praxes[1][1], praxes[3][1]);
        print ("model= VERT_tube.mnc ");
        break;
    }
    

    print ("-translation %f %f %f ",centroid[1],centroid[2],centroid[3]);
    print ("-rotation    %f %f %f \n",xrot, yrot, zrot);
    
    
    max_len = MAX3( lengths[0], lengths[1], lengths[2] );
    min_len = MIN3( lengths[0], lengths[1], lengths[2] );
    
    FREE(ang);
    
}


char *prog_name;
main(int argc, char *argv[])
{
  VIO_Volume vol;
  float **cov, **cov2, **praxes, *cog, volume_of_object;
  int i,j,obj_number;
  VIO_Real number_of_objects;
  
  char 
      *history;

  *mask_name = NULL;
  length_limit = LONGISH_CONST;
  base_limit   = BASELIM_CONST ;
  
  vol_thres[0] = LOWER_VOLUME_LIMIT;
  vol_thres[1] = UPPER_VOLUME_LIMIT;
  clobber = FALSE;
  debug = FALSE;

  
  if (ParseArgv(&argc, argv, argTable, 0) || (argc!=2)) {
      (void) fprintf(stderr, "Usage: %s [options] <file.mnc>\n",
                     argv[0]);
      exit(EXIT_FAILURE);
   }

  history = time_stamp(argc, argv);

  prog_name = argv[0];
  input_volume(argv[1],3,default_dim_names /*(char **)NULL*/, NC_BYTE, FALSE, 0.0,0.0,
               TRUE, &vol, (minc_input_options *)NULL);


  ALLOC2D(cov,4,4);
  ALLOC2D(cov2,4,4);
  ALLOC2D(praxes,4,4);
  ALLOC(cog,4);

  number_of_objects = get_volume_real_max (vol);
  
  for(obj_number=1; obj_number<(number_of_objects+0.1; obj_number++)) 
  {
      
      if (vol_to_cov(vol, obj_number, &volume_of_object, cog, cov )) 
      {
          for(i=1; i<4; i++)
              for(j=1; j<4; j++)
                  cov2[i][j] = cov[i][j];
          
          /* get the principal axis corresponding to the covariance
             matrix, note that the principal axes are in columns! */
          cov_to_praxes(3, cov2, praxes);

          /* now test the principal axis to judge 'tube-ness' */

          if (debug) {

              print ("\nStarting object %d\n", obj_number);
              
              if (this_is_a_tube(volume_of_object, cog, praxes )) 
                  (void)print ("object %d is a tube\n",obj_number);
              else
                  (void)print ("object %d is not considered a tube\n",obj_number);

              (void)print ("vol: %f \n",volume_of_object);
              (void)print ("cog: %f %f %f\n",cog[1],cog[2],cog[3]);
              
              print ("covar:\n");
              for(i=1; i<4; i++) {
                  for(j=1; j<4; j++) 
                      print ("%f ", cov[i][j] );
                  print ("\n");              
              }
              
              print ("praxes:\nc1: ");
              for(i=1; i<4; i++) {
                  for(j=1; j<4; j++) 
                      print ("%f ", praxes[j][i] );
                  if (i<3) 
                      print ("\nc%d: ",i+1);
                  else
                      print ("\n");              
              }
          }
          
          if (this_is_a_tube(volume_of_object, cog, praxes )) 
          {          

              print ("Object %d: ",obj_number);             

              write_out_tube_options( cog, praxes );
              
          }
          else 
          {
              remove_voxels(vol, obj_number);
          }
          
          
      }
      else 
      {
          (void)print ("Error in vol_to_cov for object %d\n", obj_number);
          remove_voxels(vol, obj_number);

      }
  }
  
  if (*mask_name != NULL) 
  {
      if (file_exists(mask_name) && !clobber) 
      {
          
          (void)fprintf (stderr,"File %s exists, use -clobber to overwrite.\n");
          exit( 1 );
      }
      
      if (output_modified_volume(mask_name, NC_UNSPECIFIED, FALSE, 
                                 0.0, 0.0, vol, argv[1], history, 
                                 (minc_output_options *)NULL) != OK) 
      {
          print_error_and_line_num("problems writing output mask...",__FILE__, __LINE__);
      }
  }

  FREE2D(cov);
  FREE2D(cov2);
  FREE2D(praxes);
  FREE(cog);

}
