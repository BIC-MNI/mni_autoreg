/*
------------------------------ MNI Header ----------------------------------
@NAME       : extras.c
@DESCRIPTION: this file contains a number of extra non-deformation specific
               procedures.               
@CREATED    : Mon Nov  3, 1997 , Louis Collins
@MODIFIED   : not yet!
@VERSION    : $Id: extras.c,v 1.9 2006-11-30 09:07:32 rotor Exp $
#-----------------------------------------------------------------------------
*/


#include <config.h>                
#include <float.h>
#include <volume_io.h>        
#include <time.h>

#include "local_macros.h"

void report_time(long start_time, VIO_STR text) 
{

  VIO_Real 
    time_total;
  long 
    present_time;
  VIO_STR 
    time_total_string;
    
  present_time = time(NULL);
  
  time_total = (VIO_Real)(present_time - start_time);
  time_total_string = format_time("%g %s", time_total);
  
  print ("\n%s : %s", text, time_total_string);
  
  if (time_total > 120)
    print (" (%d seconds)\n",(int)time_total);
  else
    print ("\n");
  
  
  delete_string(time_total_string);
}


void init_the_volume_to_zero(VIO_Volume volume)
{
    int             v0, v1, v2, v3, v4;
    VIO_Real            zero;
    int    sizes[VIO_MAX_DIMENSIONS];

    get_volume_sizes(volume, sizes);

    /* figure out what 0 is */
    zero = CONVERT_VALUE_TO_VOXEL(volume,0.0);
    
    for(v0=sizes[0]; v0--; ){
       for(v1=sizes[1]; v1--; ){
          for(v2=sizes[2]; v2--; ){
             for(v3=sizes[3]; v3--; ){
                set_volume_voxel_value( volume, v0, v1, v2, v3, 0, zero );
                }
             }
          }
       }
}

VIO_Real get_volume_maximum_real_value(VIO_Volume volume)
{
    int             v0, v1, v2, v3, v4;
    VIO_Real            val,max;
    int    sizes[VIO_MAX_DIMENSIONS];
    
    get_volume_sizes(volume, sizes);
    
    max = -DBL_MAX;
    
    for(v0=sizes[0]; v0--; ){
       for(v1=sizes[1]; v1--; ){
          for(v2=sizes[2]; v2--; ){
             for(v3=sizes[3]; v3--; ){
                val = get_volume_real_value( volume, v0, v1, v2, v3, 0 );
                if(val > max){
                   max = val;
                   }
                }
             }
          }
       }

    return(max);
}


/*************************************************************************/
/* debug procedure called to save the deformation at each iterative step */

void save_data(char *basename, int i, int j,
                      VIO_General_transform *transform)
{

  VIO_Status 
    status;
  VIO_STR 
    comments,name;
  FILE *file;

  ALLOC(comments, 512);
  ALLOC(name, 512);
  (void)sprintf(comments,"step %d of %d of the non-linear estimation",i,j);
  
  (void)sprintf(name,"%s%d",basename,i);
  status = open_file_with_default_suffix(name,
                                         get_default_transform_file_suffix(),
                                         WRITE_FILE, ASCII_FORMAT, &file );
  
  if( status == OK )
    status = output_transform(file,
                              basename,
                              &i,
                              comments,
                              transform);
  
  if( status == OK )
    status = close_file( file );
  
  if (status!=OK)
    print ("Error saving %s%d\n",basename,i);

  FREE(name);
  FREE(comments);
}
