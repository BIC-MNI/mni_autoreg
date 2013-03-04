/* ----------------------------- MNI Header -----------------------------------
@NAME       : test_quad
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : status
@DESCRIPTION: Program to quad fitting routines.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Mon May  8 11:29:40 MET DST 1995  LC
@MODIFIED   : $Log: test_quad.c,v $
@MODIFIED   : Revision 1.5  2006-11-29 09:09:31  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 1.4  2005/07/20 20:45:47  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 1.3  2004/02/12 05:54:06  rotor
@MODIFIED   :  * removed public/private defs
@MODIFIED   :
@MODIFIED   : Revision 1.2  2002/03/26 14:15:34  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 1.1  1999/10/25 19:52:11  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :

@COPYRIGHT  :
              Copyright 1995 Louis Collins, McConnell Brain Imaging Centre, 
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
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Extra_progs/test_quad.c,v 1.5 2006-11-29 09:09:31 rotor Exp $";
#endif

#include <stdio.h>
#include <volume_io.h>

/* Constants */
#ifndef TRUE
#  define TRUE 1
#  define FALSE 0
#endif

VIO_BOOL return_3D_disp_from_quad_fit(VIO_Real r[3][3][3], /* the values used in the quad fit */
                                            VIO_Real *dispu, /* the displacements returned */
                                            VIO_Real *dispv, 
                                            VIO_Real *dispw);        


int main(int argc, char *argv[])
{
  float tmp[3][3][3];
  VIO_Real local_corr[3][3][3];
  VIO_Real dx,dy,dz;
  int i,j,k,flag;
  float test_val;

  while (TRUE) {
    print ("\ninput local corr data:\n");

    if (scanf ("%f %f %f   %f %f %f   %f %f %f",
           &(tmp[0][0][0]), &(tmp[1][0][0]), &(tmp[2][0][0]), 
           &(tmp[0][0][1]), &(tmp[1][0][1]), &(tmp[2][0][1]), 
           &(tmp[0][0][2]), &(tmp[1][0][2]), &(tmp[2][0][2])) != 9) {
      print ("error reading first line"); exit;
    }
    if (scanf ("%f %f %f   %f %f %f   %f %f %f",
           &(tmp[0][1][0]), &(tmp[1][1][0]), &(tmp[2][1][0]), 
           &(tmp[0][1][1]), &(tmp[1][1][1]), &(tmp[2][1][1]), 
           &(tmp[0][1][2]), &(tmp[1][1][2]), &(tmp[2][1][2])) != 9) {
      print ("error reading second line"); exit;
    }
    if (scanf ("%f %f %f   %f %f %f   %f %f %f",
           &(tmp[0][2][0]), &(tmp[1][2][0]), &(tmp[2][2][0]), 
           &(tmp[0][2][1]), &(tmp[1][2][1]), &(tmp[2][2][1]), 
           &(tmp[0][2][2]), &(tmp[1][2][2]), &(tmp[2][2][2])) != 9) {
      print ("error reading third line"); exit;
    }

    for(i=0; i<3; i++)
      for(j=0; j<3; j++)
        for(k=0; k<3; k++)
          local_corr[i][j][k] = (VIO_Real)tmp[i][j][k];

    flag = return_3D_disp_from_quad_fit(local_corr, &dx, &dy, &dz);

    if (flag) 
      print ("TRUE ");
    else
      print ("FALSE");
    print ("disp: %f %f %f\n\n", dx,dy,dz);



  }
}
