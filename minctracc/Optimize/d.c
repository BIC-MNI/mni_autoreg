

public void get_average_warp_of_neighbours(Volume dx, Volume dy, Volume dz, 
					   int i,int j,int k,
					   Real *mx, Real *my, Real *mz)
{
  int 
    count, sizes[3];
  Real 
    ri,rj,rk,
    voxel,
    val_x, val_y, val_z,
    tx, ty, tz;
  
  get_volume_sizes(dx, sizes);

  *mx = 0.0; *my = 0.0; *mz = 0.0;
  ri = (Real)i; rj = (Real)j; rk = (Real)k;
  count = 0;

  if (sizes[0]>1) {
    convert_3D_voxel_to_world(dx, ri+1, rj, rk, &tx, &ty, &tz);
    GET_VOXEL_3D(voxel, dx, i+1,j,k); val_x = CONVERT_VOXEL_TO_VALUE(dx , voxel ); 
    GET_VOXEL_3D(voxel, dy, i+1,j,k); val_y = CONVERT_VOXEL_TO_VALUE(dy , voxel ); 
    GET_VOXEL_3D(voxel, dz, i+1,j,k); val_z = CONVERT_VOXEL_TO_VALUE(dz , voxel ); 
    *mx += tx+val_x; *my += ty+val_y; *mz+= tz+val_z;
    
    convert_3D_voxel_to_world(dx, ri-1, rj, rk, &tx, &ty, &tz);
    GET_VOXEL_3D(voxel, dx, i-1,j,k); val_x = CONVERT_VOXEL_TO_VALUE(dx , voxel ); 
    GET_VOXEL_3D(voxel, dy, i-1,j,k); val_y = CONVERT_VOXEL_TO_VALUE(dy , voxel ); 
    GET_VOXEL_3D(voxel, dz, i-1,j,k); val_z = CONVERT_VOXEL_TO_VALUE(dz , voxel ); 
    *mx += tx+val_x; *my += ty+val_y; *mz+= tz+val_z;
    count += 2;
  }

  if (sizes[1]>1) {
    convert_3D_voxel_to_world(dy, ri, rj+1, rk, &tx, &ty, &tz);
    GET_VOXEL_3D(voxel, dx, i,j+1,k); val_x = CONVERT_VOXEL_TO_VALUE(dx , voxel ); 
    GET_VOXEL_3D(voxel, dy, i,j+1,k); val_y = CONVERT_VOXEL_TO_VALUE(dy , voxel ); 
    GET_VOXEL_3D(voxel, dz, i,j+1,k); val_z = CONVERT_VOXEL_TO_VALUE(dz , voxel ); 
    *mx += tx+val_x; *my += ty+val_y; *mz+= tz+val_z;
    
    convert_3D_voxel_to_world(dy, ri, rj-1, rk, &tx, &ty, &tz);
    GET_VOXEL_3D(voxel, dx, i,j-1,k); val_x = CONVERT_VOXEL_TO_VALUE(dx , voxel ); 
    GET_VOXEL_3D(voxel, dy, i,j-1,k); val_y = CONVERT_VOXEL_TO_VALUE(dy , voxel ); 
    GET_VOXEL_3D(voxel, dz, i,j-1,k); val_z = CONVERT_VOXEL_TO_VALUE(dz , voxel ); 
    *mx += tx+val_x; *my += ty+val_y; *mz+= tz+val_z;
    count += 2;
  }

  if (sizes[2]>1) {
    convert_3D_voxel_to_world(dz, ri, rj, rk+1, &tx, &ty, &tz);
    GET_VOXEL_3D(voxel, dx, i,j,k+1); val_x = CONVERT_VOXEL_TO_VALUE(dx , voxel ); 
    GET_VOXEL_3D(voxel, dy, i,j,k+1); val_y = CONVERT_VOXEL_TO_VALUE(dy , voxel ); 
    GET_VOXEL_3D(voxel, dz, i,j,k+1); val_z = CONVERT_VOXEL_TO_VALUE(dz , voxel ); 
    *mx += tx+val_x; *my += ty+val_y; *mz+= tz+val_z;
    
    convert_3D_voxel_to_world(dz, ri, rj, rk-1, &tx, &ty, &tz);
    GET_VOXEL_3D(voxel, dx, i,j,k-1); val_x = CONVERT_VOXEL_TO_VALUE(dx , voxel ); 
    GET_VOXEL_3D(voxel, dy, i,j,k-1); val_y = CONVERT_VOXEL_TO_VALUE(dy , voxel ); 
    GET_VOXEL_3D(voxel, dz, i,j,k-1); val_z = CONVERT_VOXEL_TO_VALUE(dz , voxel ); 
    *mx += tx+val_x; *my += ty+val_y; *mz+= tz+val_z;
    count += 2;
  }
  
  *mx /= count; *my /= count; *mz /= count; /* average warp of neighbours, in target space  */
}

				/* add additional to current, return answer in additional */

public void add_additional_warp_to_current(Volume adx, Volume ady, Volume adz,
					   Volume dx, Volume dy, Volume dz,
					   Real weight)
{
  int
    loop_start[3], loop_end[3],
    i,j,k,sizes[3];
  Real 
    value, value2, voxel;

  get_volume_sizes(adx, sizes);

  for_less(i,0,3) {
    if (sizes[i]>3) {
      loop_start[i] = 1;
      loop_end[i] = sizes[i]-1;
    }
    else {
      loop_start[i]=0;
      loop_end[i] = sizes[i];
    }
  }


  for_less(i,loop_start[0],loop_end[0]) {
    for_less(j,loop_start[1],loop_end[1]) {
      for_less(k,loop_start[2],loop_end[2]){
	
	
	GET_VOXEL_3D(voxel, adx, i,j,k); value = CONVERT_VOXEL_TO_VALUE( adx, voxel ); 
	GET_VOXEL_3D(voxel, dx,  i,j,k); value2 = CONVERT_VOXEL_TO_VALUE( dx, voxel ); value2 += value*weight;
	voxel = CONVERT_VALUE_TO_VOXEL(adx, value2);
	SET_VOXEL_3D(adx, i,j,k, voxel); 

	GET_VOXEL_3D(voxel, ady, i,j,k); value = CONVERT_VOXEL_TO_VALUE( ady, voxel ); 
	GET_VOXEL_3D(voxel, dy,  i,j,k); value2 = CONVERT_VOXEL_TO_VALUE( dy, voxel ); value2 += value*weight;
	voxel = CONVERT_VALUE_TO_VOXEL(ady, value2);
	SET_VOXEL_3D(ady, i,j,k, voxel); 
		
	GET_VOXEL_3D(voxel, adz, i,j,k); value = CONVERT_VOXEL_TO_VALUE( adz, voxel ); 
	GET_VOXEL_3D(voxel, dz,  i,j,k); value2 = CONVERT_VOXEL_TO_VALUE( dz, voxel ); value2 += value*weight;
	voxel = CONVERT_VALUE_TO_VOXEL(adz, value2);
	SET_VOXEL_3D(adz, i,j,k, voxel); 
		
      }
    }
  }
  
}
					    

public void clamp_warp_deriv(Volume dx, Volume dy, Volume dz)
{
  int clamp_needed, clamped_once, i,j,k,sizes[3];
  Real steps[3], diff, voxel, value1, value2, old2 ;
  progress_struct
    progress;

  get_volume_sizes(dx, sizes);
  get_volume_separations(dx, steps);
  
  initialize_progress_report( &progress, FALSE, sizes[0]+sizes[1]+sizes[2] + 1,
			     "Deriv check" );

  clamped_once = FALSE;

  if (sizes[0]>2)
  for_less(j,0,sizes[1]) {
    for_less(k,0,sizes[2]){      

      
      do {
	
	clamp_needed = FALSE;
	
	GET_VOXEL_3D(voxel, dz, 0 ,j,k); value1 = CONVERT_VOXEL_TO_VALUE( dz, voxel ); 
	
	for_less(i,1,sizes[0]) {
	  GET_VOXEL_3D(voxel, dz, i,j,k);   value2 = CONVERT_VOXEL_TO_VALUE( dz, voxel ); 
	  old2 = value2;
	  
	  diff = value1 - (value2+steps[0]);
	  
	  if (diff>0.0) {
	    clamp_needed = TRUE; clamped_once = TRUE;
	    
	    value1 -= diff * DERIV_FRAC;
	    value2 += diff * DERIV_FRAC;
	    
	    voxel = CONVERT_VALUE_TO_VOXEL( dz, value1); SET_VOXEL_3D(dz, i-1,j,k, voxel);
	    voxel = CONVERT_VALUE_TO_VOXEL( dz, value2); SET_VOXEL_3D(dz, i,j,k, voxel);
	  }
	  
	  value1 = old2;
	}

      } while (clamp_needed); 
    }

    update_progress_report( &progress, j+1 );
  }


  if (sizes[1]>2)
  for_less(k,0,sizes[2]){      
    for_less(i,0,sizes[0]) {

      
      do { 
	
      clamp_needed = FALSE;

	GET_VOXEL_3D(voxel, dy, i ,0,k); value1 = CONVERT_VOXEL_TO_VALUE( dy, voxel ); 
	
	for_less(j,1,sizes[1]) {
	  GET_VOXEL_3D(voxel, dy, i,j,k);   value2 = CONVERT_VOXEL_TO_VALUE( dy, voxel ); 
	  old2 = value2;
	  
	  diff = value1 - (value2+steps[1]);
	  
	  if (diff>0.0) {
	    clamp_needed = TRUE; clamped_once = TRUE;

	    value1 -= diff * DERIV_FRAC;
	    value2 += diff * DERIV_FRAC;
	    
	    voxel = CONVERT_VALUE_TO_VOXEL( dy, value1); SET_VOXEL_3D(dy, i,j-1,k, voxel);
	    voxel = CONVERT_VALUE_TO_VOXEL( dy, value2); SET_VOXEL_3D(dy, i,j,k, voxel);
	  }
	  
	  value1 = old2;
	}

      } while (clamp_needed); 
    }
    update_progress_report( &progress, sizes[1]+k+1 );
  }


  if (sizes[2]>2)
  for_less(i,0,sizes[0]) {
    for_less(j,0,sizes[1]) {

      
      do {
	
      clamp_needed = FALSE;

	GET_VOXEL_3D(voxel, dx, i ,j,0); value1 = CONVERT_VOXEL_TO_VALUE( dx, voxel ); 
	
	for_less(k,0,sizes[2]){      
	  GET_VOXEL_3D(voxel, dx, i,j,k);   value2 = CONVERT_VOXEL_TO_VALUE( dx, voxel ); 
	  old2 = value2;
	  
	  diff = value1 - (value2+steps[2]);
	  
	  if (diff>0.0) {
	    clamp_needed = TRUE; clamped_once = TRUE;

	    value1 -= diff * DERIV_FRAC;
	    value2 += diff * DERIV_FRAC;
	    
	    voxel = CONVERT_VALUE_TO_VOXEL( dx, value1); SET_VOXEL_3D(dx, i,j,k-1, voxel);
	    voxel = CONVERT_VALUE_TO_VOXEL( dx, value2); SET_VOXEL_3D(dx, i,j,k, voxel);
	  }
	  
	  value1 = old2;
	}

      } while (clamp_needed);
    }
    update_progress_report( &progress, sizes[2]+sizes[1]+i+1 );

  }


  if (clamped_once)
    print ("Clamped once!\n");

}


/*******************************************************************************/
/*  procedure: smooth_the_warp

    desc: this procedure will smooth the current warp stored in warp_d? and 
          return the smoothed warp in smooth_d?

    meth: smoothing is accomplished by averaging the 6 neighbour of each node
          with the value at that node.

	  new_val = FRAC1*old_val + (1-FRAC1)*sum_6_neighbours(val);
    
*/
public void smooth_the_warp(Volume smooth_dx, Volume smooth_dy, Volume smooth_dz, 
			    Volume warp_dx, Volume warp_dy, Volume warp_dz,
			    Volume warp_mag, Real thres) 
{
  
  int 
    i,j,k,loop_start[3], loop_end[3],sizes[3];
  Real 
    mag,
    wx, wy, wz, mx, my, mz, tx,ty,tz,
    voxel;
  progress_struct
    progress;

  get_volume_sizes(warp_dx, sizes);
  
  for_less(i,0,3) {
    if (sizes[i]>3) {
      loop_start[i] = 1;
      loop_end[i] = sizes[i]-1;
    }
    else {
      loop_start[i]=0;
      loop_end[i] = sizes[i];
    }
  }


  initialize_progress_report( &progress, FALSE, 
			     (loop_end[0]-loop_start[0])*(loop_end[1]-loop_start[1]) + 1,
			     "Smoothing deformations" );

  for_less(i,loop_start[0],loop_end[0]) {
    for_less(j,loop_start[1],loop_end[1]) {
      for_less(k,loop_start[2],loop_end[2]){

	GET_VOXEL_3D(voxel, warp_mag, i,j,k); 
	mag = CONVERT_VOXEL_TO_VALUE(warp_mag , voxel ); 
	if (mag >= thres) {
	  
	  convert_3D_voxel_to_world(warp_dx, i, j, k, &tx, &ty, &tz);
	  
	  GET_VOXEL_3D(voxel, warp_dx, i,j,k); 
	  wx = CONVERT_VOXEL_TO_VALUE( warp_dx, voxel ); 
	  GET_VOXEL_3D(voxel, warp_dy, i,j,k); 
	  wy = CONVERT_VOXEL_TO_VALUE( warp_dy, voxel ); 
	  GET_VOXEL_3D(voxel, warp_dz, i,j,k); 
	  wz = CONVERT_VOXEL_TO_VALUE( warp_dz, voxel ); 
	  
	  tx += wx; ty += wy; tz += wz;
	  
	  get_average_warp_of_neighbours(warp_dx, warp_dy, warp_dz, 
					 i,j,k,
					 &mx, &my, &mz);
	  
	  wx += (1.0 - FRAC1)*(mx - tx);
	  wy += (1.0 - FRAC1)*(my - ty);
	  wz += (1.0 - FRAC1)*(mz - tz);
	  
	  voxel = CONVERT_VALUE_TO_VOXEL(smooth_dx, wx);
	  SET_VOXEL_3D(smooth_dx, i,j,k, voxel); 
	  voxel = CONVERT_VALUE_TO_VOXEL(smooth_dy, wy);
	  SET_VOXEL_3D(smooth_dy, i,j,k, voxel); 
	  voxel = CONVERT_VALUE_TO_VOXEL(smooth_dz, wz);
	  SET_VOXEL_3D(smooth_dz, i,j,k, voxel); 
	}
	else {
	  GET_VOXEL_3D(voxel, warp_dx, i,j,k); SET_VOXEL_3D(smooth_dx, i,j,k, voxel); 
	  GET_VOXEL_3D(voxel, warp_dy, i,j,k); SET_VOXEL_3D(smooth_dy, i,j,k, voxel); 
	  GET_VOXEL_3D(voxel, warp_dz, i,j,k); SET_VOXEL_3D(smooth_dz, i,j,k, voxel); 
	}
	  
      }
      update_progress_report( &progress,
			     (loop_end[1]-loop_start[1])*(i-loop_start[0])+(j-loop_start[1])+1  );

    }
    terminate_progress_report( &progress );
  }
}
