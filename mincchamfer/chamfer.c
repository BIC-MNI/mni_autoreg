/* ----------------------------- MNI Header -----------------------------------
@NAME       : chamfer.c
@DESCRIPTION: routines for computing a chamfer distance transform from a 
              binary mask volume.
@CREATED    : Nov 2, 1998 - louis
@MODIFIED   : 
---------------------------------------------------------------------------- */


#include <internal_volume_io.h>

extern int verbose;
extern int debug;

private void build_mask(Volume vol, Real mask_f[3][3][3], Real mask_b[3][3][3]);

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
public Status compute_chamfer(Volume chamfer) 
{

   Real
      mask_f[3][3][3],
      mask_b[3][3][3],
      steps[MAX_DIMENSIONS],
      zero, min, max_val, edge,
      vox_min, vox_max,
      val, val2;
   int
      sizes[MAX_DIMENSIONS],
      i,j,k,
      ind0,ind1,ind2;
   
   progress_struct 
      progress;
   
   Status 
      stat;

   stat = OK;
   
   get_volume_sizes(chamfer, sizes);

   get_volume_voxel_range(chamfer, &vox_min, &vox_max);

   zero = CONVERT_VALUE_TO_VOXEL(chamfer,0.0);
   
   /* init chamfer to be  binary valued with 0.0 on the object,
      and infinity (or vox_max) elsewhere */
   
   if (debug) print ("initing chamfer vol (%d %d %d)\n",sizes[0],sizes[1],sizes[2]);
   for_less (ind0, 0, sizes[0]) {
      for_less (ind1, 0, sizes[1]) {
         for_less (ind2, 0, sizes[2]) {
            
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

                                /* figure out the maximum possible distance */
   get_volume_separations(chamfer, steps);
   edge = sqrt(sizes[0]*steps[0]*sizes[0]*steps[0] + 
               sizes[1]*steps[1]*sizes[1]*steps[1]);

   max_val = sqrt(sizes[2]*steps[2]*sizes[2]*steps[2] + edge*edge);

   max_val = 254.0;
   set_volume_real_range(chamfer, 0.0, max_val);
   
   zero =CONVERT_VALUE_TO_VOXEL(chamfer,0.0);

   if (verbose) initialize_progress_report( &progress, TRUE, sizes[0], 
                                           "forward pass");
   
   for_less (ind0, 1, sizes[0]-1) {
      for_less (ind1, 1, sizes[1]-1) {
         for_less (ind2, 1, sizes[2]-1) {

            GET_VALUE_3D(val, chamfer, ind0, ind1, ind2);
            
            if (val != zero) { /* then apply forward mask */

/*
if (ind0==40 && ind1==32 && ind2==32) {
   print ("Hi!\n");
}
*/                            
               min = val;

               for_inclusive(i, -1, +1) {
                  for_inclusive(j, -1, +1) {
                     for_inclusive(k, -1, +1) {
                        
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

   return (OK);

   

   if (verbose) initialize_progress_report( &progress, TRUE, sizes[0], 
                                           "reverse pass");
       
   for_down (ind0,  sizes[0]-2, 1) {
      for_down (ind1,  sizes[1]-2, 1) {
         for_down (ind2,  sizes[2]-2, 1) {
            
            GET_VALUE_3D(val, chamfer, ind0, ind1, ind2);
            
            if (val != zero) { /* then apply backwardmask */
                                
               min = val;

               for_inclusive(i, -1, 1) {
                  for_inclusive(j, -1, +1) {
                     for_inclusive(k, -1, +1) {
                        
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


private void build_mask(Volume vol, 
                        Real mask_f[3][3][3], 
                        Real mask_b[3][3][3])
{
   int i,j,k;
   Real 
      steps[MAX_DIMENSIONS],
      d1, d2;

   get_volume_separations(vol, steps);

   for_inclusive(i,-1,1)
      for_inclusive(j,-1,1)
         for_inclusive(k,-1,1) {
            d1 = sqrt(i*i*steps[0]*steps[0] + j*j*steps[1]*steps[1]);
            mask_f[i+1][j+1][k+1] = sqrt( d1*d1 + k*k*steps[2]*steps[2] );
            mask_b[i+1][j+1][k+1] = mask_f[i+1][j+1][k+1];
         }
   
   mask_f[1][1][1] = mask_b[1][1][1] = 0.0;

   for_inclusive(j,-1,1)
      for_inclusive(k,-1,1) {
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

      for_down(i,2,0) {
         print ("slice %d\n",i);
         for_inclusive(j,0,2) {
            for_inclusive(k,0,2) 
               print ("%10.4f ",mask_f[i][j][k]);
            print ("       ");
            for_inclusive(k,0,2) 
               print ("%10.4f ",mask_b[i][j][k]);
            print ("\n");
         }
         print ("\n");
      }
   }
}





