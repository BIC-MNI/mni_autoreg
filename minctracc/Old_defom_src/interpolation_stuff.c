/*--------------------------------------------------------*/

find_N_closest(x,y,z)
     float
        x,y,z;
{
   int i;
   APOINT a_point;


   /* initialize list of closest neighbours... */
   for ( i = 0; i < N_NEIGHBOURS; i++ ) minima[i] = FAR_AWAY;
   
   for ( i = 0; i < N_NEIGHBOURS; i++ ) {
      closest_neighbour[i].xyz[0] = FAR_AWAY;
      closest_neighbour[i].xyz[1] = FAR_AWAY;
      closest_neighbour[i].xyz[2] = FAR_AWAY;
      closest_neighbour[i].err = FAR_AWAY;
   }

   a_point.xyz[0] = x;
   a_point.xyz[1] = y;
   a_point.xyz[2] = z;



   (void) find_closest( a_point.xyz );

}

interpolate_displacement(x,y,z,dx,dy,dz)
     float
        *x,*y,*z;
     DATA
        *dx,*dy,*dz;
{
   int 
      i;
   Voxel_value_type c;
   float
      tmpxyz[3],
      tx,ty,tz,
      tot_x,tot_y,tot_z;
   
   tot_x = 0.0;
   tot_y = 0.0;
   tot_z = 0.0;
   
   for ( i = 0; i < N_NEIGHBOURS; i++ ) {


      my_f_to_v(closest_neighbour[i].xyz[0],closest_neighbour[i].xyz[1],closest_neighbour[i].xyz[2],
                &(tmpxyz[0]),&(tmpxyz[1]),&(tmpxyz[2]),dx);
      c = U_interpolation(dx,tmpxyz,'f');
      tx= get_float_DATA(c, dx);
      
      my_f_to_v(closest_neighbour[i].xyz[0],closest_neighbour[i].xyz[1],closest_neighbour[i].xyz[2],
                &(tmpxyz[0]),&(tmpxyz[1]),&(tmpxyz[2]),dy);
      c = U_interpolation(dy,tmpxyz,'f');
      ty = get_float_DATA(c, dy);
      
      my_f_to_v(closest_neighbour[i].xyz[0],closest_neighbour[i].xyz[1],closest_neighbour[i].xyz[2],
                &(tmpxyz[0]),&(tmpxyz[1]),&(tmpxyz[2]),dz);
      c = U_interpolation(dz,tmpxyz,'f');
      tz = get_float_DATA(c, dz);

      tot_x += tx;
      tot_y += ty;
      tot_z += tz;
         

      printf ("%3d: %7.3f (%7.3f,%7.3f,%7.3f) %7.3f %7.3f %7.3f\n",i,minima[i],
        closest_neighbour[i].xyz[0],closest_neighbour[i].xyz[1],closest_neighbour[i].xyz[2],
        tx,ty,tz);

   }

   *x = tot_x / (float)N_NEIGHBOURS;
   *y = tot_y / (float)N_NEIGHBOURS;
   *z = tot_z / (float)N_NEIGHBOURS;
}


/* given r,c and dist<dist of last point in list,
   insert the new point in the proper place in the closest list. */
void 
insert_closest(closest, nn, r,c,dist)
     float
        **closest;
     int 
        nn,  r,c;
     float 
        dist;
{
   int
      i,j;
   
   for (i=1; i<=nn; ++i)                /* find position in the list */
      if (dist<closest[i][3]) 
         break;

   for (j=nn-1; j>=i; --j)        {        /* make space in the list */
      closest[j+1][1] = closest[j][1];
      closest[j+1][2] = closest[j][2];
      closest[j+1][3] = closest[j][3];
   }

   closest[i][1] = (float)c;            /* insert into list */
   closest[i][2] = (float)r;
   closest[i][3] = dist;
}


/* given row,col and r,c,  
   check to see is r,c is an orig point, and if so,
      calc dist of row,col to r,c and insert r,c,dist into the list */
void         
check_closest(closest, nn, row,col, r,c,s, d)
     float
        **closest;
     int 
        nn, row,col ,r,c,  s;
     DATA 
        *d;
{
   int 
      i,dr,dc;
   Voxel_value_type
      t;
   float
      dist;

   t = get_voxel_value(c,r,s, d);

   if (t==0) 
      return;

   dr = (row-r)*d->pixel_size_row;
   dc = (col-c)*d->pixel_size_col;
   
   dist = dr*dr + dc*dc;

   if (dist<closest[nn][3])
      insert_closest(closest, nn, r,c,dist);
}

void
get_value_from_nn_closest(closest, nn, cx,cy,s, dx,dy)
     float
        **closest;
     int 
        nn,s;
     Voxel_value_type
        *cx,*cy;
     DATA
        *dx,*dy;
{
   Voxel_value_type
      c_dx, c_dy;
   float
      val_x,val_y,
      frac, total;
   int 
      i;

   total = 0.0;
   for (i=1; i<=nn; i++) {
      closest[i][3] = sqrtf(closest[i][3]);
      total += closest[i][3];
   }

   val_x = 0.0;
   val_y = 0.0;

   for (i=1; i<=nn; i++) {

/*
      c_dx = dx->voxels + (s*dx->slice_size) + 
         (ROUND(closest[i][1])*dx->cols + ROUND(closest[i][2]))*dx->bytes_per_voxel;
      c_dy = dy->voxels + (s*dy->slice_size) +
         (ROUND(closest[i][1])*dy->cols + ROUND(closest[i][2]))*dy->bytes_per_voxel;
*/
      c_dx = get_voxel_value(ROUND(closest[i][1]),ROUND(closest[i][2]),s, dx);
      c_dy = get_voxel_value(ROUND(closest[i][1]),ROUND(closest[i][2]),s, dy);
      
      frac = closest[i][3]/total;
      val_x += frac*(float)(c_dx);
      val_y += frac*(float)(c_dy);
   }

   *cx = (Voxel_value_type) ROUND(val_x);
   *cy = (Voxel_value_type) ROUND(val_y);
   
}

void interpolate_from_nn_closest (row,col,cx,cy, nn, slice, dx,dy,dz)
     int
        row,col,nn,slice;
     Voxel_value_type
        *cx,*cy;
     DATA
        *dx,*dy,*dz;
{
   Voxel_value_type
      c_dx,
      c_dy,
      c_dz;
   int
      i,j,r1,r2,c1,c2,r,c,step;
   float      
      **closest;  /* row,col,dist */
 
   c_dz = get_voxel_value(col, row, slice, dz);
             /*dz->voxels + (slice*dz->slice_size) + (row*dz->cols + col)*dz->bytes_per_voxel;*/

   if (c_dz>0) {
      *cx = get_voxel_value(col, row, slice, dx);
      *cy = get_voxel_value(col, row, slice, dy);
      return;
   }
   
   closest = matrix( 1,nn,1,3);

   for (i=1; i<=nn; ++i) {        /* init the closest list */
      closest[i][1] = 9999.0;          /* col */
      closest[i][2] = 9999.0;          /* row */
      closest[i][3] = 99999999.0; /* dist */
   }

   for(step=1; step<dx->rows; ++step) {
      
      r1 = row-step; r2 = row+step;
      c1 = col-step; c2 = col+step;
       
      if (r1<0)          r1 = 0;
      if (r2>dx->rows-1) r2 = dx->rows-1;
      if (c1<0)          c1 = 0;
      if (c2>dx->cols-1) c2 = dx->cols-1;

                                /* check a ring of pixels around row,col */
                                /* and place values in closest list      */
      for (r=r1; r<=r2; ++r)                
         check_closest(closest, nn, row,col,r,c1, slice, dz);
      for (r=r1; r<=r2; ++r)
         check_closest(closest, nn, row,col,r,c2, slice, dz);
      for (c=c1+1; c<c2; ++c)
         check_closest(closest, nn, row,col,r1,c, slice, dz);
      for (c=c1+1; c<c2; ++c)
         check_closest(closest, nn, row,col,r2,c, slice, dz);


      /* end when the farthest of the closest list is closer than closest pixel checked.*/
      if (closest[nn][3]< step*MIN(dx->pixel_size_row,dx->pixel_size_col))
         break;
   }

   get_value_from_nn_closest(closest, nn, cx,cy, slice, dx,dy);

   free_matrix(closest,      1,nn,1,3);
}


/* 
  this routine will do 2D interpolation of the vector field

  the deform_dz will be set to zero.
      deform_dx and deform_dy will return with interpolated values.

*/
void interpolate_deform_slice(deform_dx,deform_dy,deform_dz)
     DATA 
        *deform_dx,*deform_dy,*deform_dz;
{
   Voxel_value_type
      cx,cy,
      c_dx, c_dy, c_dz,
      t_dx, t_dy, t_dz;
    int 
      nn,origs,total,resamps, i,j,k,
      row,col,slice;

   origs = 0;
   total = 0;
   resamps = 0;

  
/*
   count how many have to get done, and use the dz volume to flag them.
*/
   
   slice = 0; 
   for (row = 0; row < deform_dx->rows; row++) {
            
      for (col = 0; col < deform_dx->cols; col++) {
         
         c_dx = get_voxel_value(col, row, slice, deform_dx);
         c_dy = get_voxel_value(col, row, slice, deform_dy);

         total++;

         if ((c_dx != 0) || (c_dy != 0) ) {
            origs++;
            put_voxel_value( (Voxel_value_type)1, col, row, slice, deform_dz);
         } 
         else
            put_voxel_value( (Voxel_value_type)0, col, row, slice, deform_dz);

      }
   }

   printf ("there are %d original samples out of %d points.\n",origs,total);

/*
  now interpolate them, using the nearest neighbour strategy
*/
   if (total<1 || origs<1) {
      printf ("there are no samples.\n");
      return;
   }


   if (origs<N_NEIGHBOURS) {
      printf ("there are not enough samples (%d<%d). Resetting NN=%d\n",total,N_NEIGHBOURS,N_NEIGHBOURS);
      nn = total;
   }
   else
      nn = N_NEIGHBOURS;

   for (row = 0; row < deform_dx->rows; row++) {
      
      for (col = 0; col < deform_dx->cols; col++) {
         
         c_dz = get_voxel_value(col, row, slice, deform_dz);

         if (!(c_dz>0)) {
            
            resamps++;
            
            interpolate_from_nn_closest (row,col,&cx,&cy, nn, 0, deform_dx,deform_dy,deform_dz);

            put_voxel_value(cx, col, row, slice, deform_dx);
            put_voxel_value(cy, col, row, slice, deform_dy);

         }
      }
   }

/*
  reset the dz volume

   memset(deform_dz->voxels, 0, deform_dz->slice_size*deform_dz->slices);
*/
   
   printf ("resampled %d out of %d points.\n",resamps,total);


}




/* this routine will be used to interpolate the missing samples in the 
   vector deformation field.

   in the estimation of the deformation field, only voxels that had a gradient
   magnitude greater than 10% maximum have been calculated.  This leaves out
   a number of voxels, where no estimation has been done.  

   In this routine, the deformation value for these voxels will be estimated
   from its N_NEIGHBOURS nearest neighbors.
*/
void interpolate_deform_volume(deform_dx,deform_dy,deform_dz,labels,ndim)
     DATA 
        *deform_dx,*deform_dy,*deform_dz,*labels;
     int ndim;
{
   Voxel_value_type
     label,
     c_dx,
     c_dy,
     c_dz,
     t_dx,
     t_dy,
     t_dz;
   int 
     count, i,j,k,
     row,col,slice;
   APOINT
      *the_list,
     *the_neighbors;
   float
     tempxyz[3],
     xpos,ypos,zpos,
     xpos2,ypos2,zpos2,
     xpf,ypf,zpf,
      x_displacement,
     y_displacement,
     z_displacement;
   DATA
     *tmp_dx, *tmp_dy, *tmp_dz;
   int zlimit;
   
   VIO_progress_struct progress;
   
   VIO_Status status;
   

   if (ndim==2) {
     interpolate_deform_slice(deform_dx,deform_dy,deform_dz,labels);
     return; 
   }
   
   
   
   
   /* count how many true samples will be on the list.--------------------------------- */
   
   count = 0;
   
   if (ndim < 3)
     zlimit = 1;
   else
     zlimit = deform_dx->slices;
   
   for (slice = 0; slice < zlimit ; slice++) {
     for (row = 0; row < deform_dx->rows; row++) {
       for (col = 0; col < deform_dx->cols; col++) {
         
         label = get_voxel_value(col, row, slice, labels);

         if (label == MAPPING_FOUND)
           count++;
         
       }
      }
   }
   
   
   /* build the list of true samples  ------------------------------------------------ */

   ALLOC1(status, list, count, APOINT);
   count = 0;                        
   the_list = list;

   for (slice = 0; slice < zlimit; slice++) {
     for (row = 0; row < deform_dx->rows; row++) {
       for (col = 0; col < deform_dx->cols; col++) {
         
         label = get_voxel_value(col, row, slice, labels);
         
         if (label == MAPPING_FOUND) {
           
                 count++;
           
           my_v_to_f((float)(col), (float)(row), (float)(slice),
                     &(the_list->xyz[0]),&(the_list->xyz[1]),&(the_list->xyz[2]),
                     deform_dx);
           
           the_list++;
           
         }
         
       }
     }
   }
   

   /* make the seach data structure -------------------------------------------------*/

   qcount = count;

   divide_volume( list, count ,ndim );

/*
   for (i=0; i<MESH_SIZE; i++)
      for (j=0; j<MESH_SIZE; j++)
         for (k=0; k<MESH_SIZE; k++)
            printf ("%3d %3d %3d -> %3d\n",i,j,k,z_size[i][j][k]);

printf("after divide cubes\n");
*/
   /* go through the volume, one last time, finding the nearest neighbors, and interpolating---*/

   count = 0;

   initialize_progress_report( &progress, FALSE, deform_dx->slices*deform_dx->rows+1, 
                              "Estimating deformations");

   for (slice = 0; slice < zlimit; slice++)
      for (row = 0; row < deform_dx->rows; row++) {
         for (col = 0; col < deform_dx->cols; col++) {
                 

           label = get_voxel_value(col, row, slice, labels);
           
           if (label != MAPPING_FOUND) {
             
             count++;
             
             my_v_to_f((float)(col), (float)(row), (float)(slice), 
                       &xpf, &ypf, &zpf, deform_dx);
             
/* 1
printf ("for point %12.4f, %12.4f, %12.4f\n",xpf,ypf,zpf);
*/

             find_N_closest(xpf,ypf,zpf);
             
             interpolate_displacement(&x_displacement,&y_displacement,&z_displacement,
                                      deform_dx,deform_dy,deform_dz);
             

/* 1
printf ("and the displacement is: %12.5f, %12.5f, %12.5f\n\n",
        x_displacement,y_displacement,z_displacement);
*/
             c_dx = get_voxel_value_DATA(x_displacement, deform_dx);  
             put_voxel_value(c_dx, col, row, slice, deform_dx);
             c_dy = get_voxel_value_DATA(y_displacement, deform_dy);
             put_voxel_value(c_dy, col, row, slice, deform_dy);
             c_dz = get_voxel_value_DATA(z_displacement, deform_dz);
             put_voxel_value(c_dz, col, row, slice, deform_dz);
           }
           
         }
         update_progress_report( &progress, slice*deform_dx->rows + row +1 );
       }

   terminate_progress_report(&progress);

   printf ("resampled  %d out of %d.\n",count, deform_dx->slices*deform_dx->rows*deform_dx->cols);
}






/* 
   Divide a list of points into MESH_SIZE regions or subvolumes
   by sorting x,y,z coordinates... 

   After sorting, the space occupied by the points wAPOINTl have been partitioned into
   a cube with MESH_SIZE units per side: each one being the same shape and having equal volume...
   kind of like a Rubik's cube.
*/

void divide_volume( APOINT *lst, int lst_size, int ndim )
{
   long int x_size, y_size, incx, incy, incz;
   register i, j, k;
   float min_x, min_y, min_z;
   float max_x, max_y, max_z;

   float rangex, rangey, rangez;

/* find range of coords...  */
   min_x = min_y = min_z = 9999.0;
   max_x = max_y = max_z = -9999.0;
   for ( i = 0; i < lst_size; i++ ) {
      if ( (lst+i)->xyz[0] > max_x ) max_x = (lst+i)->xyz[0];
      else if ( (lst+i)->xyz[0] < min_x ) min_x = (lst+i)->xyz[0];
      if ( (lst+i)->xyz[1] > max_y ) max_y = (lst+i)->xyz[1];
      else if ( (lst+i)->xyz[1] < min_y ) min_y = (lst+i)->xyz[1];
      if ( (lst+i)->xyz[2] > max_z ) max_z = (lst+i)->xyz[2];
      else if ( (lst+i)->xyz[2] < min_z ) min_z = (lst+i)->xyz[2];
   }

/*
printf (" in divide cubes: range %f %f, %f %f, %f %f\n",
        min_x,max_x,min_y,max_y,min_z,max_z);
*/
   rangex = max_x - min_x;
   rangey = max_y - min_y;
   rangez = max_z - min_z;

   min_x -= rangex/(float)(MESH_SIZE-4);
   max_x += rangex/(float)(MESH_SIZE-4);
   min_y -= rangey/(float)(MESH_SIZE-4);
   max_y += rangey/(float)(MESH_SIZE-4);
   min_z -= rangez/(float)(MESH_SIZE-4);
   max_z += rangez/(float)(MESH_SIZE-4);

   for ( i = 0; i < MESH_SIZE; i++ ) {
      xdiv[i] = (max_x-min_x)*(float)(i+1)/(float)MESH_SIZE + min_x;
      ydiv[i] = (max_y-min_y)*(float)(i+1)/(float)MESH_SIZE + min_y;
      zdiv[i] = (max_z-min_z)*(float)(i+1)/(float)MESH_SIZE + min_z;

/*printf(" div: %7.2f %7.2f %7.2f \n",xdiv[i],ydiv[i],zdiv[i]);*/

   }

/* sort the points in memory with smallest X's at the top of the lst... */
   if ( lst_size > 1 ) 
      (void) qsort( (char *)lst, (size_t)qcount, (size_t)sizeof(APOINT), comparex );

/*
for (i=0; i<lst_size; i++)
   printf("%3d: %f %f %f\n",i,(lst+i)->xyz[0],(lst+i)->xyz[1],(lst+i)->xyz[2]);
*/

   incx = 0;
 
   for ( i = 0; i < MESH_SIZE; i++ ) {
      x_size = incx;
      while( (lst+incx)->xyz[0] <= xdiv[i] && incx < qcount ) ++incx;
      x_size = incx - x_size;
      if ( x_size > 1 ) (void) qsort((char *)(lst+incx-x_size), (size_t)x_size, 
                                     (size_t)sizeof(APOINT), comparey );
      incy = 0;
      for ( j = 0; j < MESH_SIZE; j++ ) {
         y_size = incy;
         while( (lst+incx-x_size+incy)->xyz[1] <= ydiv[j] && incy < x_size ) ++incy;
         y_size = incy - y_size;
         if ( y_size > 1 ) (void) qsort((char *)(lst+incx-x_size+incy-y_size), 
                                        (size_t)y_size, (size_t)sizeof(APOINT), comparez );
         incz = 0;
         for ( k = 0; k < MESH_SIZE; k++ ) {
            list_ptr[i][j][k] = incx - x_size + incy - y_size + incz;
            z_size[i][j][k] = incz;
            while( (lst+incx-x_size+incy-y_size+incz)->xyz[2] <= zdiv[k] && incz < y_size ) ++incz;
            z_size[i][j][k] = incz - z_size[i][j][k];
            /*printf("%2d %2d %2d --> %d\n",i,j,k,z_size[i][j][k]);*/
         }
      }
   }

}


int comparex( APOINT *elem1, APOINT *elem2 )
{
   if ( elem1->xyz[0] > elem2->xyz[0] ) return( 1 );
   else if ( elem1->xyz[0] < elem2->xyz[0] ) return( -1 );
   else return( 0 );
}


int comparey( APOINT *elem1, APOINT *elem2 )
{
   if ( elem1->xyz[1] > elem2->xyz[1] ) return( 1 );
   else if ( elem1->xyz[1] < elem2->xyz[1] ) return( -1 );
   else return( 0 );
}


int comparez( APOINT *elem1, APOINT *elem2 )
{
   if ( elem1->xyz[2] > elem2->xyz[2] ) return( 1 );
   else if ( elem1->xyz[2] < elem2->xyz[2] ) return( -1 );
   else return( 0 );
}



void find_closest( float a[3] )
{
   register i, j, k, m;
   int save_r1, save_r2, save_r3, r1_a, r1_b, r2_a, r2_b, r3_a, r3_b;
   int search_level, r1_min, r1_max, r2_min, r2_max, r3_min, r3_max;
   float best, this_best, b[3];

/*printf("      in find closest, a=%f,%f,%f\n",a[0],a[1],a[2]);*/


   (void) which_region( a );
   save_r1 = r1;
   save_r2 = r2;
   save_r3 = r3;

   best = FAR_AWAY;

/* find distance to closest point within sub-volume... */
   best = min_dist( a, (list+list_ptr[r1][r2][r3]), z_size[r1][r2][r3] );

/*printf("      in find closest, best=%f\n",best);*/

/* define boundaries of box already searched... */

   r1_min = r1;
   r1_max = r1;
   r2_min = r2;
   r2_max = r2;
   r3_min = r3;
   r3_max = r3;

/* search... */

   search_level = 0;

/*printf("      in find closest, r1,r2,r3 = %d %d %d\n",r1,r2,r3);*/

   while( minima[N_NEIGHBOURS-1] == FAR_AWAY ) {
      
      ++search_level;

/*printf("      in find closest, s .. r1,r2,r3 = %d   %d %d %d\n",search_level,r1,r2,r3);*/

      r1_a = save_r1 - search_level;
      r2_a = save_r2 - search_level;
      r3_a = save_r3 - search_level;
      r1_b = save_r1 + search_level;
      r2_b = save_r2 + search_level;
      r3_b = save_r3 + search_level;

      if ( r1_a < 0 ) r1_a = 0;
      if ( r1_b > (MESH_SIZE-1) ) r1_b = MESH_SIZE - 1;
      if ( r2_a < 0 ) r2_a = 0;
      if ( r2_b > (MESH_SIZE-1) ) r2_b = MESH_SIZE - 1;
      if ( r3_a < 0 ) r3_a = 0;
      if ( r3_b > (MESH_SIZE-1) ) r3_b = MESH_SIZE - 1;

      for ( i = r1_a; i <= r1_b; i++ ) {
         for ( j = r2_a; j <= r2_b; j++ ) {
            for ( k = r3_a; k <= r3_b; k++ ) {
               if ( (i>r1_a&&i<r1_b) && (j>r2_a&&j<r2_b) && (k>r3_a&&k<r3_b) ) continue;
               this_best = min_dist( a, (list+list_ptr[i][j][k]), z_size[i][j][k] );

/*
printf("      %d %d %d in find closest, this_best,best=%f %f\n",i,j,k,this_best,best);
for (m=0;m<N_NEIGHBOURS;++m) printf ("%7.3f ",minima[m]);
printf("\n");
printf("in vol = %d\n",z_size[i][j][k]);
*/
               if ( this_best < best ) {
                  best = this_best;
               }
            }
         }
      }

/* define boundaries of box already searched... */

      r1_min = r1_a;
      r1_max = r1_b;
      r2_min = r2_a;
      r2_max = r2_b;
      r3_min = r3_a;
      r3_max = r3_b;
   }

/* increase allowable volume to be searched to prevent errors due to sub-volume boundary effects... */

   b[0] = a[0] - best;
   b[1] = a[1] - best;
   b[2] = a[2] - best;
   
   (void) which_region( b );
   r1_a = r1;
   r2_a = r2;
   r3_a = r3;

   b[0] = a[0] + best;
   b[1] = a[1] + best;
   b[2] = a[2] + best;
   
   (void) which_region( b );
   r1_b = r1;
   r2_b = r2;
   r3_b = r3;

   for ( i = r1_a; i <= r1_b; i++ ) {
      for ( j = r2_a; j <= r2_b; j++ ) {
         for ( k = r3_a; k <= r3_b; k++ ) {
            if ( (i>=r1_min&&i<=r1_max) && (j>=r2_min&&j<=r2_max) && (k>=r3_min&&k<=r3_max) ) continue;
            this_best = min_dist( a, (list+list_ptr[i][j][k]), z_size[i][j][k] );
            if ( this_best < best ) {
               best = this_best;
            }
         }
      }
   }
   
   return;
}


/* determine which region a point is in... */

void which_region( float xyz[3] )
{
   register i;

   r1 = r2 = r3 = -1;

   for ( i = (MESH_SIZE-1); i >= 0; i-- ) {
      if ( xyz[0] <= xdiv[i] ) r1 = i;
      if ( xyz[1] <= ydiv[i] ) r2 = i;
      if ( xyz[2] <= zdiv[i] ) r3 = i;
   }

   if ( r1 == -1 ) r1 = (MESH_SIZE-1);
   if ( r2 == -1 ) r2 = (MESH_SIZE-1);
   if ( r3 == -1 ) r3 = (MESH_SIZE-1);

   return;
}

/* 
   Find the closest points in the given sub-region...

   return: worst of the closest neighbours...
*/



float min_dist( float xyz[3], APOINT *search_list, long int num_in_search_list )
{
   long int i;
   register j, k;
   long int save_i[N_NEIGHBOURS];
   float dx, dy, dz, dist;

   for ( i = 0; i < num_in_search_list; i++ ) {
      dx = xyz[0] - (search_list+i)->xyz[0];
      dy = xyz[1] - (search_list+i)->xyz[1];
      dz = xyz[2] - (search_list+i)->xyz[2];
      dist = dx*dx + dy*dy + dz*dz;
      if ( dist <= 0.000001 ) continue;   /* don't consider self as closest point... */
      for ( j = 0; j < N_NEIGHBOURS; j++ ) {
         if ( minima[N_NEIGHBOURS-j-1] < dist ) break;
      }
      for ( k = 0; k < (j-1); k++ ) {
         minima[N_NEIGHBOURS-k-1] = minima[N_NEIGHBOURS-k-2];
         closest_neighbour[N_NEIGHBOURS-k-1].xyz[0] = closest_neighbour[N_NEIGHBOURS-k-2].xyz[0];
         closest_neighbour[N_NEIGHBOURS-k-1].xyz[1] = closest_neighbour[N_NEIGHBOURS-k-2].xyz[1];
         closest_neighbour[N_NEIGHBOURS-k-1].xyz[2] = closest_neighbour[N_NEIGHBOURS-k-2].xyz[2];
         closest_neighbour[N_NEIGHBOURS-k-1].err = closest_neighbour[N_NEIGHBOURS-k-2].err;
      }
      if ( j > 0 ) {
         minima[N_NEIGHBOURS-j] = dist;
         closest_neighbour[N_NEIGHBOURS-j].xyz[0] = (search_list+i)->xyz[0];
         closest_neighbour[N_NEIGHBOURS-j].xyz[1] = (search_list+i)->xyz[1];
         closest_neighbour[N_NEIGHBOURS-j].xyz[2] = (search_list+i)->xyz[2];
         closest_neighbour[N_NEIGHBOURS-j].err = (search_list+i)->err;
      }
      else if ( (dz = dz*dz) > minima[N_NEIGHBOURS-1] ) break;
   }
   
   return( minima[N_NEIGHBOURS-1] );
}




