/* ----------------------------- MNI Header -----------------------------------
@NAME       : randomize_tags.c
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : status
@DESCRIPTION: program to add a random displacement to the 
              target tag list.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Tue Aug 23 15:32:30 EST 1994 - Louis
@MODIFIED   : $Log: randomize_tags.c,v $
@MODIFIED   : Revision 1.4  2006-11-29 09:09:31  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 1.3  2005/07/20 20:45:46  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 1.2  2002/03/26 14:15:32  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 1.1  1999/10/25 19:52:09  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :


---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="";
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <volume_io.h>

#include <sys/types.h>
#include <time.h>


time_t time(time_t *tloc);

char *prog_name;

VIO_Real gaussian_random_w_std(VIO_Real sigma);

/* Main program */

int main(int argc, char *argv[])
{
  union
    {
      long   l;
      char   c[4];
    } seedval;
  
  time_t t;
  char tmp;
  
  char 
    *in_tag, *out_tag, **labels;
  int 
    i,j,n_vols, n_points,
    *struc_ids, *pat_ids;
  VIO_Real
    std,**tags1, **tags2, 
    *weights;

  t = 2*time(NULL);
  seedval.l = t; 
  
  tmp = seedval.c[0]; seedval.c[0] = seedval.c[3]; seedval.c[3] = tmp; 
  tmp = seedval.c[1]; seedval.c[1] = seedval.c[2]; seedval.c[2] = tmp;
  
  srand48(seedval.l);
  
  if (argc!=4) {
    (void) fprintf(stderr, "Usage: randomize_tags in.tag out.tag std\n",
                   argv[0]);
      exit(EXIT_FAILURE);
  }

  prog_name = argv[0];
  in_tag    = argv[1];
  out_tag   = argv[2];
  std       = atof( argv[3] );

  input_tag_file(in_tag, 
                 &n_vols, &n_points, 
                 &tags1, &tags2, 
                 &weights, &struc_ids, &pat_ids, &labels);

  for(i=0; i<n_points; i++)
    for(j=0; j<3; j++)
      tags1[i][j] += gaussian_random_w_std(std);
  
  output_tag_file(out_tag, "Target (1st vol) points have been randomized",
                 n_vols, n_points, 
                 tags1, tags2, 
                 weights, struc_ids, pat_ids, labels);
                 
  FREE2D(tags1);
  FREE2D(tags2);
  FREE(weights);
  FREE(struc_ids);
  FREE(pat_ids);
  FREE2D(labels);

  exit(EXIT_SUCCESS);
}
