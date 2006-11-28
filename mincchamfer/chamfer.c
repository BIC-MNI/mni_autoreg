/* ----------------------------- MNI Header -----------------------------------
@NAME       : chamfer.c
@DESCRIPTION: routines for computing a chamfer distance transform from a 
              binary mask volume.
@CREATED    : Nov 2, 1998 - louis
@MODIFIED   : 
---------------------------------------------------------------------------- */

#include <config.h>
#include <volume_io.h>

#define  MIN( x, y )  ( ((x) <= (y)) ? (x) : (y) )

extern int verbose;
extern int debug;

static void build_mask(VIO_Volume vol, VIO_Real mask_f[3][3][3], VIO_Real mask_b[3][3][3]);

/* ----------------------------- MNI Header -----------------------------------
@NAME       :  compute_chamfer
@INPUT/OUTPUT: chamfer
                   The original input volume will be destroyed and replaced
                   with the resulting chamfer distance volume.  The chamfer
                   will contains 0's where the mask was and values > 0.0
                   for all other voxels, where the value is an estimate of
                   the distance to the the nearest voxel of the mask.
@RETURNS    : ERROR if error, OK otherwise
@DESCRIPTION: Uses an idea from georges, who got it from claire, 
              who got it from Borgefors
@GLOBALS    : 
@CALLS      : 
@CREATED    : Nov 2, 1998 Louis
@MODIFIED   : 
---------------------------------------------------------------------------- */
VIO_Status compute_chamfer(VIO_Volume chamfer, VIO_Real max_val)
{

   VIO_Real
      mask_f[3][3][3],
      mask_b[3][3][3],
      zero, min,
      vox_min, vox_max,
      val, val2;
   int
      sizes[VIO_MAX_DIMENSIONS],
      i,j,k,
      ind0,ind1,ind2;
   
   VIO_progress_struct 
      progress;
   
   get_volume_sizes(chamfer, sizes);

   get_volume_voxel_range(chamfer, &vox_min, &vox_max);

   zero = CONVERT_VALUE_TO_VOXEL(chamfer,0.0);
   
   /* init chamfer to be  binary valued with 0.0 on the object,
      and infinity (or vox_max) elsewhere */
   
   if (debug) print ("initing chamfer vol (%d %d %d)\n",sizes[0],sizes[1],sizes[2]);
   for(ind0=0; ind0<sizes[0]; ind0++) {
      for(ind1=0; ind1<sizes[1]; ind1++) {
         for(ind2=0; ind2<sizes[2]; ind2++) {
            
            GET_VOXEL_3D(val, chamfer, ind0, ind1, ind2);
            if (val == zero) {
               SET_VOXEL_3D(chamfer, ind0, ind1, ind2, vox_max ); 
            }            
            else {
               SET_VOXEL_3D(chamfer, ind0, ind1, ind2, vox_min); 
            }
            
         }
      }
   }

   if (debug) print ("building mask\n");

   build_mask(chamfer, mask_f, mask_b);

   set_volume_real_range(chamfer, 0.0, max_val);
   zero = CONVERT_VALUE_TO_VOXEL(chamfer,0.0);

   if (verbose) initialize_progress_report( &progress, TRUE, sizes[0], 
                                           "forward pass");
   
   for(ind0=1; ind0<sizes[0]-1; ind0++) {
      for(ind1=1; ind1<sizes[1]-1; ind1++) {
         for(ind2=1; ind2<sizes[2]-1; ind2++) {

            GET_VALUE_3D(val, chamfer, ind0, ind1, ind2);
            
            if (val != zero) { /* then apply forward mask */

               min = val;

               for(i=-1; i<=+1; i++) {
                  for(j=-1; j<=+1; j++) {
                     for(k=-1; k<=+1; k++) {
                        
                        GET_VALUE_3D(val, chamfer, i+ind0, j+ind1, k+ind2);
                        
                        val2 = val + mask_f[i+1][j+1][k+1];
                        min = MIN (min, val2);
                        
                     }
                  }
               }
                           
               min = convert_value_to_voxel(chamfer, min);

               SET_VOXEL_3D(chamfer, ind0, ind1, ind2, min );              
               
            } /* if val != 0.0 */
            
         } /* ind2 */
      } /* ind1 */
      if (verbose) update_progress_report( &progress, ind0+1 );
   } /* ind0 */
   
   if (verbose) terminate_progress_report( &progress );


   if (verbose) initialize_progress_report( &progress, TRUE, sizes[0], 
                                           "reverse pass");
       
   for(ind0=sizes[0]-2; ind0>=1; ind0--) {
      for(ind1=sizes[1]-2; ind1>=1; ind1--) {
         for(ind2=sizes[2]-2; ind2>=1; ind2--) {
            
            GET_VALUE_3D(val, chamfer, ind0, ind1, ind2);
            
            if (val != zero) { /* then apply backwardmask */
                                
               min = val;

               for(i=-1; i<=1; i++) {
                  for(j=-1; j<=+1; j++) {
                     for(k=-1; k<=+1; k++) {
                        
                        GET_VALUE_3D(val, chamfer, i+ind0, j+ind1, k+ind2);
                        
                        val2 = val + mask_b[i+1][j+1][k+1];
                        min = MIN (min, val2);
                        
                     }
                  }
               }
                           
               min = convert_value_to_voxel(chamfer, min);

               SET_VOXEL_3D(chamfer, ind0, ind1, ind2, min );              
               
            } /* if val != 0.0 */
            
         } /* ind2 */
      } /* ind1 */
      if (verbose) update_progress_report( &progress, ind0+1 );
   } /* ind0 */
   
   if (verbose) terminate_progress_report( &progress );
   
   return (OK);
   
}


static void build_mask(VIO_Volume vol, 
                        VIO_Real mask_f[3][3][3], 
                        VIO_Real mask_b[3][3][3])
{
   int i,j,k;
   VIO_Real 
      steps[VIO_MAX_DIMENSIONS],
      d1;

   get_volume_separations(vol, steps);

   for(i=-1; i<=1; i++)
      for(j=-1; j<=1; j++)
         for(k=-1; k<=1; k++) {
            d1 = sqrt(i*i*steps[0]*steps[0] + j*j*steps[1]*steps[1]);
            mask_f[i+1][j+1][k+1] = sqrt( d1*d1 + k*k*steps[2]*steps[2] );
            mask_b[i+1][j+1][k+1] = mask_f[i+1][j+1][k+1];
         }
   
   mask_f[1][1][1] = mask_b[1][1][1] = 0.0;

   for(j=-1; j<=1; j++)
      for(k=-1; k<=1; k++) {
         mask_f[2][j+1][k+1] = 1000.0;
         mask_b[0][j+1][k+1] = 1000.0;
      }
   
   mask_f[1][1][2] = 1000.0;
   mask_f[1][2][0] = 1000.0;
   mask_f[1][2][1] = 1000.0;
   mask_f[1][2][2] = 1000.0;

   mask_b[1][1][0] = 1000.0;
   mask_b[1][0][0] = 1000.0;
   mask_b[1][0][1] = 1000.0;
   mask_b[1][0][2] = 1000.0;


   if (debug) {
      print ("%20s %35s\n","forward","reverse");

      for(i=2; i>=0; i--) {
         print ("slice %d\n",i);
         for(j=0; j<=2; j++) {
            for(k=0; k<=2; k++) 
               print ("%10.4f ",mask_f[i][j][k]);
            print ("       ");
            for(k=0; k<=2; k++) 
               print ("%10.4f ",mask_b[i][j][k]);
            print ("\n");
         }
         print ("\n");
      }
   }
}





