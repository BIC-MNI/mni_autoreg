/* ----------------------------- MNI Header -----------------------------------
@NAME       : def_to_mag
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : status
@DESCRIPTION: Program to compute the magnitude of a deformation volume.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Mon Sep 30 11:35:52 EDT 1996    LC
@MODIFIED   : $Log: def_to_mag.c,v $
@MODIFIED   : Revision 1.3  2004-02-12 05:54:05  rotor
@MODIFIED   :  * removed public/private defs
@MODIFIED   :
@MODIFIED   : Revision 1.2  2002/03/26 14:15:28  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 1.1  1999/10/25 19:52:06  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :

@COPYRIGHT  :
              Copyright 1996 Louis Collins, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Extra_progs/def_to_mag.c,v 1.3 2004-02-12 05:54:05 rotor Exp $";
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <volume_io/internal_volume_io.h>
#include <Proglib.h>

/* Constants */
#ifndef TRUE
#  define TRUE 1
#  define FALSE 0
#endif

				/* type of job to compute */
#define JOB_UNDEF 0
#define JOB_SCALE 1
#define JOB_MAG   2

				/* type of neighbourhood to use */
#define NEIGHBOUR_6     1
#define NEIGHBOUR_333   2
#define NEIGHBOUR_555   3

static char *my_ZYX_dim_names[] = { MIzspace, MIyspace, MIxspace };

void print_usage_and_exit(char *pname);

void get_volume_XYZV_indices(Volume data, int xyzv[]);

BOOLEAN get_average_scale_from_neighbours(General_transform *trans,
						 int voxel[],
						 int avg_type,
						 Real *scale);

BOOLEAN get_magnitude(General_transform *trans,
                             int voxel[],
                             Real *dist);


/* Main program */
char *prog_name;

Real real_range[2];
int
    neighbour_type,
    job_type,
    verbose, 
    clobber_flag,
    debug;


static ArgvInfo argTable[] = {
{"-range",      ARGV_FLOAT,   (char *) 2,     (char *) &real_range,
     "Specify real range for output volume."},
{"-scale",      ARGV_CONSTANT, (char *) JOB_SCALE, (char *) &job_type,
     "Compute local scale (default)."},
{"-6", ARGV_CONSTANT, (char *)NEIGHBOUR_6 , (char *) &neighbour_type,
     "Use only the 6 4-connected neighbours  (default)."},
{"-333", ARGV_CONSTANT, (char *)NEIGHBOUR_333 , (char *) &neighbour_type,
     "Use the 3x3x3  neighbours."},
{"-555", ARGV_CONSTANT, (char *)NEIGHBOUR_555 , (char *) &neighbour_type,
     "Use the 5x5x5  neighbours."},
{"-magnitude",  ARGV_CONSTANT, (char *) JOB_MAG, (char *) &job_type,
     "Compute local def magnitude."},
{"-no_clobber", ARGV_CONSTANT, (char *) FALSE, (char *) &clobber_flag,
     "Do not overwrite output file (default)."},
{"-clobber",    ARGV_CONSTANT, (char *) TRUE,  (char *) &clobber_flag,
     "Overwrite output file."},
{"-verbose",    ARGV_CONSTANT, (char *) TRUE,     (char *) &verbose,
     "Write messages indicating progress (default)"},
{"-quiet",      ARGV_CONSTANT, (char *) FALSE,    (char *) &verbose,
     "Do not write log messages"},
{"-debug",      ARGV_CONSTANT, (char *) TRUE,  (char *) &debug,
     "Print out debug info."},
{NULL, ARGV_END, NULL, NULL, NULL}
};

int main(int argc, char *argv[])
{

   Volume
       output_vol,
       def_volume;
   General_transform 
       *grid_transform_ptr,def_field;
   
   Real
       value,
       voxel[MAX_DIMENSIONS],
       steps[MAX_DIMENSIONS],
       new_steps[MAX_DIMENSIONS],
       start[MAX_DIMENSIONS],
       new_start[MAX_DIMENSIONS];
   
   int
       parse_flag,
       i, trans_count,
       counter,
       ind[MAX_DIMENSIONS],
       count[MAX_DIMENSIONS],
       new_count[MAX_DIMENSIONS],
       new_xyzv[MAX_DIMENSIONS],
       xyzv[MAX_DIMENSIONS];
   
   char *infile,*outfile;
   
   Status status;
   progress_struct
       progress;

   
   verbose       = TRUE;	/* init some variables */
   clobber_flag  = FALSE;
   debug         = FALSE;
   job_type      = JOB_SCALE;
   real_range[0] =  0.0;
   real_range[1] =  3.0;
   neighbour_type= NEIGHBOUR_6;
   prog_name     = argv[0];

   
   /* Call ParseArgv to interpret all command line args (returns TRUE if error) */
   parse_flag = ParseArgv(&argc, argv, argTable, 0);

   /* Check remaining arguments */
   if (parse_flag || argc != 3) print_usage_and_exit(prog_name);

   infile   = argv[1];
   outfile  = argv[2];

   /* Check the writability of the output file */

   if (file_exists(outfile) && !clobber_flag)
   {
       (void) fprintf(stderr, "Error: file <%s> already exists,\nuse -clobber to overwrite.\n",
		      outfile);
       exit(EXIT_FAILURE);
       
   }
   
   /* Read the deformation field .xfm file */
   if (input_transform_file(infile, &def_field) != OK) {
      (void) fprintf(stderr, "%s: Error reading transform file %s\n",
                     argv[0], infile);
      exit(EXIT_FAILURE);
   }


   for_less(trans_count,0, get_n_concated_transforms(&def_field) ) {
       
       grid_transform_ptr = get_nth_general_transform(&def_field, trans_count );

       def_volume = (Volume)NULL;
       
       if (grid_transform_ptr->type == GRID_TRANSFORM) {
	   def_volume = grid_transform_ptr->displacement_volume;
       }
   }
   
   if (def_volume == (Volume)NULL) 
   {
       (void) fprintf(stderr, "Input tranformation <%s> does\n",  infile);
       (void) fprintf(stderr, "not have a GRID_TRANSFORM.\n");
       exit(EXIT_FAILURE);
   }
   
   
   get_volume_sizes(       def_volume, count);
   get_volume_separations( def_volume, steps);
   get_volume_XYZV_indices(def_volume, xyzv);
   
   output_vol = create_volume(3, my_ZYX_dim_names, NC_SHORT, TRUE, 0.0, 0.0);
   get_volume_XYZV_indices(output_vol, new_xyzv);

   for_less(i,0,3) {
     new_count[ new_xyzv[i] ] = count[ xyzv[i] ];
     new_steps[ new_xyzv[i] ] = steps[ xyzv[i] ];
   }
   for_less(i,0,MAX_DIMENSIONS) {
     start[i] = 0.0;
     voxel[i] = 0.0;
     new_start[i] = 0.0;
   }

   convert_voxel_to_world(def_volume,
			  voxel,
			  &start[0], &start[1], &start[2]);
   
   set_volume_sizes(output_vol, new_count);
   set_volume_separations(output_vol, new_steps);
   set_volume_translation(output_vol,  voxel, start);
   alloc_volume_data(output_vol);
   set_volume_real_range(output_vol, real_range[0], real_range[1]);

   for_less(i, 0, MAX_DIMENSIONS)
       ind[ xyzv[i] ] = 0;

				/* init the volume to the default value */
   switch (job_type) 
   {
   case JOB_SCALE:
       value = 1.0;
       break;
   case JOB_MAG:
       value = 0.0;
       break;
   case JOB_UNDEF:
   default:
       (void) fprintf(stderr,"No default value for job_type=%d\n",job_type);
       exit( EXIT_FAILURE );
   }
   
   
   for_less(ind[ xyzv[X] ], 0, count[  xyzv[X] ]) 
       for_less(ind[ xyzv[Y] ], 0, count[  xyzv[Y] ]) 
	   for_less(ind[ xyzv[Z] ], 0, count[  xyzv[Z] ]) {
	       
	       set_volume_real_value(output_vol,
				     ind[ xyzv[Z] ], ind[ xyzv[Y] ], ind[ xyzv[X] ], 0, 0,
				     value);
	       
	   }
   
				/* now compute the scale factor */
   
   if (verbose) initialize_progress_report(&progress, FALSE, 
                                           (count[xyzv[X]]-2)*(count[xyzv[Y]]-2)*(count[xyzv[Z]]-2) + 1,
                                           "Computing");

   counter = 0;
   
   for_less(ind[ xyzv[X] ], 1, count[  xyzv[X] ]-1) 
      for_less(ind[ xyzv[Y] ], 1, count[  xyzv[Y] ]-1) 
         for_less(ind[ xyzv[Z] ], 1, count[  xyzv[Z] ]-1) {
            
            
            switch (job_type) {
               case JOB_SCALE:
	       if (get_average_scale_from_neighbours(grid_transform_ptr,
						     ind, neighbour_type, 
                                                     &value)) {
                  set_volume_real_value(output_vol,
                                        ind[ xyzv[Z] ], ind[ xyzv[Y] ], ind[ xyzv[X] ], 
                                        0, 0,
                                        value);
               }
               break;
               case JOB_MAG:
	       if (get_magnitude(grid_transform_ptr,
                                 ind, &value)) {
                  set_volume_real_value(output_vol,
                                        ind[ xyzv[Z] ], ind[ xyzv[Y] ], ind[ xyzv[X] ], 
                                        0, 0,
                                        value);
               }
               break;
               case JOB_UNDEF:
               default:
               (void) fprintf(stderr,"No default value for job_type=%d\n",job_type);
               exit( EXIT_FAILURE );
            }
            counter++;
            
            if (verbose) 	update_progress_report( &progress, counter );
            
         }

   if (verbose) terminate_progress_report(&progress);

   
   
   if (output_volume(outfile, NC_UNSPECIFIED, FALSE, 0.0, 0.0,
		     output_vol,
		     NULL, (minc_output_options *)NULL) != OK) {
       
       (void) fprintf(stderr,"Cannot write %s\n",outfile);
       exit( EXIT_FAILURE );
       
   }

   delete_volume(output_vol);

   exit(EXIT_SUCCESS);
}

void print_usage_and_exit(char *pname) {

  (void) fprintf(stderr, "This program is used to extract local geometric information\n");
  (void) fprintf(stderr, "the GRID_TRANSFORM representation of a deformation field.\n\n");
  (void) fprintf(stderr, "Usage: %s [options] <input.xfm> <result.mnc>\n",
		 pname);
  exit(EXIT_FAILURE);

}

void get_volume_XYZV_indices(Volume data, int xyzv[])
{
  
  int 
    axis, i, vol_dims;
  char 
    **data_dim_names;

  vol_dims       = get_volume_n_dimensions(data);
  data_dim_names = get_volume_dimension_names(data);

  for_less(i,0,N_DIMENSIONS+1) xyzv[i] = -1;
  for_less(i,0,vol_dims) {
    if (convert_dim_name_to_spatial_axis(data_dim_names[i], &axis )) {
      xyzv[axis] = i; 
    } 
    else {     /* not a spatial axis */
      xyzv[Z+1] = i;
    }
  }
  delete_dimension_names(data,data_dim_names);

}


BOOLEAN get_average_scale_from_neighbours(General_transform *trans,
						 int voxel[],
						 int avg_type,
						 Real *scale)
{
  int
    start[MAX_DIMENSIONS],
    end[MAX_DIMENSIONS],
    count, i, 
    voxel2[MAX_DIMENSIONS],
    xyzv[MAX_DIMENSIONS],
    sizes[MAX_DIMENSIONS];
  Real 
      R_voxel[MAX_DIMENSIONS],
      mx, my, mz,
      dx,dy,dz,
      cx,cy,cz,
      cbx,cby,cbz,
      px,py,pz,
      dist,
      norm_dist,
      sum_dist,
      seps[MAX_DIMENSIONS],
      def_vector[N_DIMENSIONS];
  Volume volume;
  
  if (trans->type != GRID_TRANSFORM) {
    print_error_and_line_num("get_average_scale_from_neighbours not called with GRID_TRANSFORM",
			     __FILE__, __LINE__);
    return (FALSE);
  }

  volume = trans->displacement_volume;
  
  get_volume_sizes(volume, sizes);
  get_volume_separations(volume, seps);
  get_volume_XYZV_indices(volume, xyzv);
  count = 0;
  mx = 0.0; my = 0.0; mz = 0.0; /* assume no warp vector in volume */

				/* make sure voxel is in volume */
  for_less(i,0,3 ) {
    if (voxel[ xyzv[i]]<0 || voxel[ xyzv[i] ]>=sizes[ xyzv[i]] ) {
       return (FALSE);
    }
  }
  
  for_less(i,0,MAX_DIMENSIONS ) { /* copy the voxel position */
    voxel2[i] = voxel[i];
  }
  

  for_less(i,0, MAX_DIMENSIONS) R_voxel[i] = voxel[i];
  convert_voxel_to_world(volume, R_voxel, &cbx, &cby, &cbz);

  for_less(voxel2[ xyzv[Z+1] ], 0, sizes[ xyzv[Z+1] ]) {
      def_vector[voxel2[ xyzv[Z+1] ]] = 
	  get_volume_real_value(volume,
				voxel2[0],voxel2[1],voxel2[2],voxel2[3],voxel2[4]);
  }
  cx = cbx + def_vector[X];	     
  cy = cby + def_vector[Y];
  cz = cbz + def_vector[Z];
  
  switch (avg_type) {
  case 1:  {			/* get the 6 neighbours, 2 along
				   each of the spatial axes,  if they exist 
				   ie the 4-connected immediate neighbours only */

      sum_dist = 0.0;
      
      for_less(i, 0 , N_DIMENSIONS) {
	  
	  if ((voxel[ xyzv[i] ]+1) < sizes[ xyzv[i] ]) {
	      
	      voxel2[ xyzv[i] ] = voxel[ xyzv[i] ] + 1;
	      
	      for_less(voxel2[ xyzv[Z+1] ], 0, sizes[ xyzv[Z+1] ]) {
		  def_vector[voxel2[ xyzv[Z+1] ]] = 
		      get_volume_real_value(volume,
					    voxel2[0],voxel2[1],voxel2[2],voxel2[3],voxel2[4]);
	      }
	      
	      for_less(i,0, MAX_DIMENSIONS) R_voxel[i] = voxel2[i];
	      convert_voxel_to_world(volume, R_voxel, &px, &py, &pz);

	      norm_dist = sqrt ( (cbx-px)*(cbx-px) + (cby-py)*(cby-py) + (cbz-pz)*(cbz-pz) );

	      px += def_vector[X];	     
	      py += def_vector[Y];
	      pz += def_vector[Z];
	      	      
	      dist = sqrt ( (cx-px)*(cx-px) + (cy-py)*(cy-py) + (cz-pz)*(cz-pz) ) / 
		  norm_dist;
	      

	      voxel2[ xyzv[i] ] = voxel[ xyzv[i] ];
		  
	      mx += def_vector[X]; my += def_vector[Y]; mz += def_vector[Z];
	      ++count;

	      sum_dist += dist;
	  }

	  
	  if ((voxel[ xyzv[i] ]-1) >= 0) {
	      voxel2[ xyzv[i] ] = voxel[ xyzv[i] ] - 1;
	      
	      for_less(voxel2[ xyzv[Z+1] ], 0, sizes[ xyzv[Z+1] ]) {
		  def_vector[voxel2[ xyzv[Z+1] ]] = 
		      get_volume_real_value(volume,
					    voxel2[0],voxel2[1],voxel2[2],voxel2[3],voxel2[4]);
	      }
	      
	      for_less(i,0, MAX_DIMENSIONS) R_voxel[i] = voxel2[i];
	      convert_voxel_to_world(volume, R_voxel, &px, &py, &pz);

	      norm_dist = sqrt ( (cbx-px)*(cbx-px) + (cby-py)*(cby-py) + (cbz-pz)*(cbz-pz) );

	      px += def_vector[X];	     
	      py += def_vector[Y];
	      pz += def_vector[Z];
	      	      
	      dist = sqrt ( (cx-px)*(cx-px) + (cy-py)*(cy-py) + (cz-pz)*(cz-pz) ) / 
		  norm_dist;
	      
	      
	      voxel2[ xyzv[i] ] = voxel[ xyzv[i] ];
	      
	      mx += def_vector[X]; my += def_vector[Y]; mz += def_vector[Z];
	      ++count;
	      sum_dist += dist;
	  }

	  
      }
      break;
  }
  case 2: {			/* 3x3x3 ie the 26 8-connected
				   immediate neighbours only */
    
    for_less(i, 0 , N_DIMENSIONS) {
      start[i] = voxel[ xyzv[i] ] - 1;
      if (start[i]<0) start[i]=0;
      end[i] = voxel[ xyzv[i] ] + 1;
      if (end[i]>sizes[ xyzv[i] ]-1) end[i] = sizes[ xyzv[i] ]-1;
    }
    
    sum_dist = 0.0;    
    for_inclusive(voxel2[ xyzv[X] ],start[X], end[X])
      for_inclusive(voxel2[ xyzv[Y] ],start[Y], end[Y])
	for_inclusive(voxel2[ xyzv[Z] ],start[Z], end[Z]) {

	  if ((voxel2[ xyzv[X]] != voxel[ xyzv[X] ]) ||
	      (voxel2[ xyzv[Y]] != voxel[ xyzv[Y] ]) ||
	      (voxel2[ xyzv[Z]] != voxel[ xyzv[Z] ])) {
	    for_less(voxel2[ xyzv[Z+1] ], 0, sizes[ xyzv[Z+1] ]) {
	      def_vector[voxel2[ xyzv[Z+1] ]] = 
		get_volume_real_value(volume,
				      voxel2[0],voxel2[1],voxel2[2],voxel2[3],voxel2[4]);
	    }

	    for_less(i,0, MAX_DIMENSIONS) R_voxel[i] = voxel2[i];
	    convert_voxel_to_world(volume, R_voxel, &px, &py, &pz);

	    norm_dist = sqrt ( (cbx-px)*(cbx-px) + (cby-py)*(cby-py) + (cbz-pz)*(cbz-pz) );
	    
	    px += def_vector[X];	     
	    py += def_vector[Y];
	    pz += def_vector[Z];
	    
	    dist = sqrt ( (cx-px)*(cx-px) + (cy-py)*(cy-py) + (cz-pz)*(cz-pz) ) / 
		norm_dist;
	    	    
	    mx += def_vector[X]; my += def_vector[Y]; mz += def_vector[Z];
	    ++count;
	    sum_dist += dist;
	  }
	}
    break;
  }
  case 3: {			/* 5x5x5 */
    for_less(i, 0 , N_DIMENSIONS) {
      start[i] = voxel[ xyzv[i] ] - 2;
      if (start[i]<0) start[i]=0;
      end[i] = voxel[ xyzv[i] ] + 2;
      if (end[i]>sizes[ xyzv[i] ]-1) end[i] = sizes[ xyzv[i] ]-1;
    }
    
    sum_dist = 0.0;    
    for_inclusive(voxel2[ xyzv[X] ],start[X], end[X])
      for_inclusive(voxel2[ xyzv[Y] ],start[Y], end[Y])
	for_inclusive(voxel2[ xyzv[Z] ],start[Z], end[Z]) {
	  
	  if ((voxel2[ xyzv[X]] != voxel[ xyzv[X] ]) ||
	      (voxel2[ xyzv[Y]] != voxel[ xyzv[Y] ]) ||
	      (voxel2[ xyzv[Z]] != voxel[ xyzv[Z] ])) {

	    for_less(voxel2[ xyzv[Z+1] ], 0, sizes[ xyzv[Z+1] ]) {
	      def_vector[voxel2[ xyzv[Z+1] ]] = 
		get_volume_real_value(volume,
				      voxel2[0],voxel2[1],voxel2[2],voxel2[3],voxel2[4]);
	    }

	    for_less(i,0, MAX_DIMENSIONS) R_voxel[i] = voxel2[i];
	    convert_voxel_to_world(volume, R_voxel, &px, &py, &pz);

	    norm_dist = sqrt ( (cbx-px)*(cbx-px) + (cby-py)*(cby-py) + (cbz-pz)*(cbz-pz) );
	    
	    px += def_vector[X];	     
	    py += def_vector[Y];
	    pz += def_vector[Z];

	    
	    dist = sqrt ( (cx-px)*(cx-px) + (cy-py)*(cy-py) + (cz-pz)*(cz-pz) ) / 
		norm_dist;
	    
	    mx += def_vector[X]; my += def_vector[Y]; mz += def_vector[Z];
	    ++count;
	    sum_dist += dist;

	  }
	}
    break;
    
  }
  }    

  if (count>0) {		/* average scale value */
    *scale = sum_dist / count; 
    return(TRUE);
  }
  else {
    return(FALSE);
  }
}


BOOLEAN get_magnitude(General_transform *trans,
                             int voxel[],
                             Real *dist)
{
  int
    start[MAX_DIMENSIONS],
    end[MAX_DIMENSIONS],
    count, i, 
    voxel2[MAX_DIMENSIONS],
    xyzv[MAX_DIMENSIONS],
    sizes[MAX_DIMENSIONS];
  Real 
      seps[MAX_DIMENSIONS],
      def_vector[N_DIMENSIONS];
  Volume volume;
  
  if (trans->type != GRID_TRANSFORM) {
    print_error_and_line_num("get_magnitude not called with GRID_TRANSFORM",
			     __FILE__, __LINE__);
    return (FALSE);
  }

  volume = trans->displacement_volume;
  
  get_volume_sizes(volume, sizes);
  get_volume_separations(volume, seps);
  get_volume_XYZV_indices(volume, xyzv);
  count = 0;

				/* make sure voxel is in volume */
  for_less(i,0,3 ) {
    if (voxel[ xyzv[i]]<0 || voxel[ xyzv[i] ]>=sizes[ xyzv[i]] ) {
       return (FALSE);
    }
  }
  
  
  for_less(voxel[ xyzv[Z+1] ], 0, sizes[ xyzv[Z+1] ]) {
     def_vector[voxel[ xyzv[Z+1] ]] = 
        get_volume_real_value(volume,
                              voxel[0],voxel[1],voxel[2],voxel[3],voxel[4]);
  }
  *dist = sqrt ((def_vector[X]*def_vector[X]) + 
                (def_vector[Y]*def_vector[Y]) + 
                (def_vector[Z]*def_vector[Z]));
  return(TRUE);

}


