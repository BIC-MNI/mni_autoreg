/*
# ------------------------------ MNI Header ----------------------------------
#@NAME       : vox_space.c
#@DESCRIPTION: routines to map global lattice information from
               world coordinates to voxel coordinates.

               These routines are used only for linear registration (with
               all objective functions except mutual information.

               It would be desirable to have similar functionality (and the
               equivalent speed up) for the non-linear code.

#@CREATED    : Wed Jun 25, 1997, Louis Collins
#@MODIFIED   : not yet!
#@VERSION    : $Id: vox_space.c,v 1.5 2002-03-26 14:15:46 stever Exp $
#----------------------------------------------------------------------------- */

#include <volume_io/internal_volume_io.h>
#include "constants.h"
#include "arg_data.h"
#include "vox_space.h"
#include "local_macros.h"

public Voxel_space_struct* new_voxel_space_struct(void) {

   Voxel_space_struct *vox_space;

   ALLOC(vox_space, 1);
   
   ALLOC(vox_space->voxel_to_voxel_space, 1);

   create_linear_transform(vox_space->voxel_to_voxel_space, (Transform *)NULL);

   return(vox_space);
}

public void delete_voxel_space_struct( Voxel_space_struct *vox_space) {

   delete_general_transform(vox_space->voxel_to_voxel_space);
   FREE(vox_space->voxel_to_voxel_space);

   FREE(vox_space);

}

public void build_reorder_matrix_vox2xyz(General_transform *trans, Volume volume) {
   
   int axis;
   Transform *lin;

   lin = get_linear_transform_ptr(trans);

   axis = volume->spatial_axes[X];
   if( axis >= 0 ) {
      Transform_elem( *lin, 0, axis ) = 1.0;
      Transform_elem( *lin, 1, axis ) = 0.0;
      Transform_elem( *lin, 2, axis ) = 0.0;
   }
   axis = volume->spatial_axes[Y];
   if( axis >= 0 ) {
      Transform_elem( *lin, 0, axis ) = 0.0;
      Transform_elem( *lin, 1, axis ) = 1.0;
      Transform_elem( *lin, 2, axis ) = 0.0;
   }
   axis = volume->spatial_axes[Z];
   if( axis >= 0 ) {
      Transform_elem( *lin, 0, axis ) = 0.0;
      Transform_elem( *lin, 1, axis ) = 0.0;
      Transform_elem( *lin, 2, axis ) = 1.0;
   }
   
}

public void build_reorder_matrix_xyz2vox(General_transform *trans, Volume volume) {
   
   int axis;
   Transform *lin;

   lin = get_linear_transform_ptr(trans);

   axis = volume->spatial_axes[X];
   if( axis >= 0 ) {
      Transform_elem( *lin, axis, 0 ) = 1.0;
      Transform_elem( *lin, axis, 1 ) = 0.0;
      Transform_elem( *lin, axis, 2 ) = 0.0;
   }
   axis = volume->spatial_axes[Y];
   if( axis >= 0 ) {
      Transform_elem( *lin, axis, 0 ) = 0.0;
      Transform_elem( *lin, axis, 1 ) = 1.0;
      Transform_elem( *lin, axis, 2 ) = 0.0;
   }
   axis = volume->spatial_axes[Z];
   if( axis >= 0 ) {
      Transform_elem( *lin, axis, 0 ) = 0.0;
      Transform_elem( *lin, axis, 1 ) = 0.0;
      Transform_elem( *lin, axis, 2 ) = 1.0;
   }
   
}

public void get_into_voxel_space(Arg_Data *globals,
                                 Voxel_space_struct *vox,
                                 Volume v1, Volume v2) {
   Transform 
      *lin;
   General_transform 
      *reorder,
      *w2v;
   Real tx,ty,tz;
   PointR pnt,tmp_pt;
   Real 
      voxel_vector[MAX_DIMENSIONS],
      s_voxel_xyz[MAX_DIMENSIONS],
      s_voxel[MAX_DIMENSIONS],
      s_world[N_DIMENSIONS],
      t_voxel[MAX_DIMENSIONS],
      t_world[N_DIMENSIONS];
   int i;
                                /* take care of the starting coordinate */
   convert_3D_world_to_voxel(v1,
                             globals->start[X], globals->start[Y], globals->start[Z],
                             &vox->start[0],    &vox->start[1],    &vox->start[2]);

                                /* take care of the directions required to step
                                   through the volume */
   for_less(i,0,3) {
      convert_world_vector_to_voxel(v1,
                                    RVector_x(globals->directions[i]),
                                    RVector_y(globals->directions[i]),
                                    RVector_z(globals->directions[i]),
                                    voxel_vector);
      fill_Vector(vox->directions[i], voxel_vector[0],voxel_vector[1],voxel_vector[2]);
   }

                                /* take care of the
                                     voxel-world * world-world * world-voxel
                                   transformation */
   ALLOC(w2v,1);
   ALLOC(reorder,1);
   create_linear_transform(reorder, (Transform *)NULL); /* build identity */
   

   create_inverse_general_transform(get_voxel_to_world_transform( v2 ),
                                    w2v);

   
                                /* go from voxel to xyz space */
   build_reorder_matrix_vox2xyz(reorder, v1);
   concat_general_transforms(reorder,
                             get_voxel_to_world_transform( v1 ), 
                             vox->voxel_to_voxel_space );
   concat_general_transforms(vox->voxel_to_voxel_space,
                             globals->trans_info.transformation,
                             vox->voxel_to_voxel_space );
   concat_general_transforms(vox->voxel_to_voxel_space, w2v, 
                             vox->voxel_to_voxel_space);
   build_reorder_matrix_xyz2vox(reorder, v2);
   concat_general_transforms(vox->voxel_to_voxel_space, reorder, 
                             vox->voxel_to_voxel_space);
   FREE(w2v);
   FREE(reorder);



   if (FALSE && globals->flags.debug) {

      lin = get_linear_transform_ptr(globals->trans_info.transformation);

      print ("global trans:\n");
      print ("    %9.3f %9.3f %9.3f %9.3f\n", 
             Transform_elem( *lin, 0, 0 ),
             Transform_elem( *lin, 0, 1 ),
             Transform_elem( *lin, 0, 2 ),
             Transform_elem( *lin, 0, 3 ));
      print ("    %9.3f %9.3f %9.3f %9.3f\n", 
             Transform_elem( *lin, 1, 0 ),
             Transform_elem( *lin, 1, 1 ),
             Transform_elem( *lin, 1, 2 ),
             Transform_elem( *lin, 1, 3 ));
      print ("    %9.3f %9.3f %9.3f %9.3f\n", 
             Transform_elem( *lin, 2, 0 ),
             Transform_elem( *lin, 2, 1 ),
             Transform_elem( *lin, 2, 2 ),
             Transform_elem( *lin, 2, 3 ));

      print ("start:  w- %9.3f %9.3f %9.3f -> v- %9.3f %9.3f %9.3f\n",
             globals->start[X], globals->start[Y], globals->start[Z],
             vox->start[0],     vox->start[1],     vox->start[2]);

      for_less(i,0,3) {
         print ("dir[%d]: w- %9.3f %9.3f %9.3f -> v- %9.3f %9.3f %9.3f\n", 
                i,
                RVector_x(globals->directions[i]),
                RVector_y(globals->directions[i]),
                RVector_z(globals->directions[i]),
                RVector_x(vox->directions[i]),
                RVector_y(vox->directions[i]),
                RVector_z(vox->directions[i]));
      }

                                /* define a point in the source volume */
      s_voxel_xyz[0] = 10.0; 
      s_voxel_xyz[1] = 15.0; 
      s_voxel_xyz[2] = 20.0; 
      s_voxel_xyz[3] = 0.0; 
      s_voxel_xyz[4] = 0.0; 

      reorder_xyz_to_voxel( v1, s_voxel_xyz, s_voxel );
      convert_voxel_to_world(v1,s_voxel,&s_world[0],&s_world[1],&s_world[2]);
      fill_Point(pnt, s_world[0], s_world[1], s_world[2]);

      print ("source: w- %9.3f %9.3f %9.3f -> v- %9.3f %9.3f %9.3f\n",
             s_world[0], s_world[1], s_world[2],
             s_voxel[0], s_voxel[1], s_voxel[2]);


      DO_TRANSFORM(tmp_pt, globals->trans_info.transformation, pnt);

      convert_3D_world_to_voxel(v2, Point_x(tmp_pt), Point_y(tmp_pt), Point_z(tmp_pt), 
                                &tx, &ty, &tz);

      print ("classical method:\n");
      print ("target: v- %9.3f %9.3f %9.3f -> w- %9.3f %9.3f %9.3f\n",
             tx, ty, tz,
             Point_x(tmp_pt), Point_y(tmp_pt), Point_z(tmp_pt));


      general_transform_point(vox->voxel_to_voxel_space,
                              s_voxel[0],  s_voxel[1],  s_voxel[2],
                              &t_voxel[0], &t_voxel[1], &t_voxel[2]);

      convert_3D_voxel_to_world(v2, 
                                t_voxel[0],  t_voxel[1],  t_voxel[2],
                                &t_world[0], &t_world[1], &t_world[2]);
                                

      print ("direct method:\n");
      
      print ("target: v- %9.3f %9.3f %9.3f -> w- %9.3f %9.3f %9.3f\n",
             t_voxel[0], t_voxel[1], t_voxel[2],
             t_world[0], t_world[1], t_world[2]);
             

   }
}

public  void  my_homogenous_transform_point(
    Transform  *transform,
    Real       x,
    Real       y,
    Real       z,
    Real       w,
    Real       *x_trans,
    Real       *y_trans,
    Real       *z_trans )
{
    Real       w_trans;

    *x_trans = Transform_elem(*transform,0,0) * x +
               Transform_elem(*transform,0,1) * y +
               Transform_elem(*transform,0,2) * z +
               Transform_elem(*transform,0,3) * w;

    *y_trans = Transform_elem(*transform,1,0) * x +
               Transform_elem(*transform,1,1) * y +
               Transform_elem(*transform,1,2) * z +
               Transform_elem(*transform,1,3) * w;

    *z_trans = Transform_elem(*transform,2,0) * x +
               Transform_elem(*transform,2,1) * y +
               Transform_elem(*transform,2,2) * z +
               Transform_elem(*transform,2,3) * w;

    w_trans =  Transform_elem(*transform,3,0) * x +
               Transform_elem(*transform,3,1) * y +
               Transform_elem(*transform,3,2) * z +
               Transform_elem(*transform,3,3) * w;

    if( w_trans != 0.0 && w_trans != 1.0 )
    {
        *x_trans /= w_trans;
        *y_trans /= w_trans;
        *z_trans /= w_trans;
    }
}
