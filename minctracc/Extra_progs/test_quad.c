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
@MODIFIED   : Revision 1.1  1999-10-25 19:52:11  louis
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
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Extra_progs/test_quad.c,v 1.1 1999-10-25 19:52:11 louis Exp $";
#endif

#include <stdio.h>
#include <internal_volume_io.h>

/* Constants */
#ifndef TRUE
#  define TRUE 1
#  define FALSE 0
#endif
#ifndef public
#  define public
#  define private static
#endif

public BOOLEAN return_3D_disp_from_quad_fit(Real r[3][3][3], /* the values used in the quad fit */
					    Real *dispu, /* the displacements returned */
					    Real *dispv, 
					    Real *dispw);	


int main(int argc, char *argv[])
{
  float tmp[3][3][3];
  Real local_corr[3][3][3];
  Real dx,dy,dz;
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

    for_less(i,0,3)
      for_less(j,0,3)
	for_less(k,0,3)
	  local_corr[i][j][k] = (Real)tmp[i][j][k];

    flag = return_3D_disp_from_quad_fit(local_corr, &dx, &dy, &dz);

    if (flag) 
      print ("TRUE ");
    else
      print ("FALSE");
    print ("disp: %f %f %f\n\n", dx,dy,dz);



  }
}
