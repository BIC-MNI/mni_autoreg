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

/* ----------------------------- MNI Header -----------------------------------
@NAME       :  compute_chamfer
@INPUT      : orig         original mask volume - will be destroyed by call
                   I assume that the volume is binary valued - 1 for the
                   mask, and 0 for the background.

@OUTPUT     : chamfer      resulting chamfer distance volume.  This volume
                   will contains 1's where the mask was and values > 1 for
                   all other voxels, where the value is proportional to the
                   distance of that voxel to the nearest mask region.
@RETURNS    : ERROR if error, OK otherwise
@DESCRIPTION: This routine dilates each mask region, adding the dilated
              shell to the mask volume, thus estimating the distance for
              each voxel in the volume to the edge of the masked regions.
              
              Each dilation can be thought of as an onion skin wrapped
              around each region in the mask.  Since each dilation (skin)
              is numbered, and this number is ctored in the volume, we have
              a direct look-up of the distance from any voxel to the edge
              of the object represented in the masked region.
@GLOBALS    : 
@CALLS      : 
@CREATED    : Nov 2, 1998 Louis
@MODIFIED   : 
---------------------------------------------------------------------------- */
public Status compute_chamfer(Volume orig, Volume  chamfer) 
{

   Real
      zero, val, val2;
   int
      pass,                       /* this is the number of dilation passes */
      sizes[3],
      done,
      found,
      i,j,k,
      ind0,ind1,ind2;
   
   progress_struct 
      progress;
   
   Status 
      stat;

   Volume
      workvol;

   stat = OK;
   
   get_volume_sizes(orig, sizes);
   zero =CONVERT_VALUE_TO_VOXEL(orig,0.0);
   
   
   /* make working volume */
   workvol = copy_volume(orig);
   
   /* init orig vol to be binary valued and
      init chamfer to 0.0 distance every where */
   
   if (debug) print ("initing data and chamfer vols (%d %d %d)\n",sizes[0],sizes[1],sizes[2]);
   for_less (ind0, 0, sizes[0]) {
      for_less (ind1, 0, sizes[1]) {
         for_less (ind2, 0, sizes[2]) {
            
                                /* ensure binary {0,1} vol for orig and working vol */
            GET_VOXEL_3D(val, orig, ind0, ind1, ind2);
            if (val == zero) {
               SET_VOXEL_3D(orig,    ind0, ind1, ind2, 0.0 ); 
               SET_VOXEL_3D(workvol, ind0, ind1, ind2, 0.0 ); 
            }            
            else {
               SET_VOXEL_3D(orig,    ind0, ind1, ind2, 1.0 ); 
               SET_VOXEL_3D(workvol, ind0, ind1, ind2, 1.0 ); 
            }
            
                                /* init chamfer to dist=0.0 */
            SET_VOXEL_3D(chamfer, ind0, ind1, ind2, 0.0 ); 
            
         }
      }
   }

   set_volume_real_range(orig,0.0,255.0);  /* note that these values won't work for large vols */
   set_volume_voxel_range(orig,0.0,255.0);
   
   set_volume_real_range(workvol,0.0,255.0);
   set_volume_voxel_range(workvol,0.0,255.0);
   
   set_volume_real_range(chamfer,0.0,255.0);
   set_volume_voxel_range(chamfer,0.0,255.0);
   
   zero =CONVERT_VALUE_TO_VOXEL(orig,0.0);

   /* now loop until all voxels in the mask are set */

   pass = 1; done = FALSE;
   
   if (verbose) print ("computing distances:\n");
   while (!done && (pass < 255 )) {
      print ("Pass %d\n",pass);
      
      done = TRUE;              /* act as if there will be no voxels to label */

      if (verbose) initialize_progress_report( &progress, TRUE, sizes[0], "computing");
      
      for_less (ind0, 1, sizes[0]-1) {
         for_less (ind1, 1, sizes[1]-1) {
            for_less (ind2, 1, sizes[2]-1) {
               
               GET_VOXEL_3D(val, orig, ind0, ind1, ind2);
         
               if (val == zero) {

                  done = FALSE; /* we found one, and will need at least one more pass */

                                /* increment the distance for this voxel */
                  SET_VOXEL_3D(chamfer, ind0, ind1, ind2, pass ); 
                  
                                /* look and see if we are near the object
                                   in the mask, and if so then dilate */
                  found = FALSE;
                  for_inclusive(i, ind0-1, ind0+1) {
                     for_inclusive(j, ind1-1, ind1+1) {
                        for_inclusive(k, ind2-1, ind2+1) {
                           
                           GET_VOXEL_3D(val2, orig, i, j, k);
                           if (val2 != zero && !(i==ind0 && 
                                                 j==ind1 && 
                                                 k==ind2)
                               ) {
                              found=TRUE;	
                              /* break;*/  
                           }
                        }
                     }
                  }
                  
                  if (found==TRUE)
                     SET_VOXEL_3D(workvol, ind0, ind1, ind2, 1.0 );  	    

               } /* if val == 0.0 */

            } /* ind2 */
         } /* ind1 */
         if (verbose) update_progress_report( &progress, ind0+1 );
      } /* ind0 */

      if (verbose) terminate_progress_report( &progress );
   
      for_less (ind0, 0, sizes[0]) {
         for_less (ind1, 0, sizes[1]) {
            for_less (ind2, 0, sizes[2]) {
               
               GET_VOXEL_3D(val, workvol, ind0, ind1, ind2);
               SET_VOXEL_3D(orig, ind0, ind1, ind2, val ); 
               
            }
         }
      }
      
      if (debug && done) print ("We're done after pass %d\n",pass);
      pass++;
      
   } /* while !done */
   
   delete_volume(workvol);

   return (OK);
   
}
