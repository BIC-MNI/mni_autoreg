/* ----------------------------- MNI Header -----------------------------------
@COPYRIGHT  :
              Copyright 1993 Louis Collins, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
@MODIFIED   : $Log: deform_support.c,v $
@MODIFIED   : Revision 1.1  1999-10-25 19:59:10  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Release_stubs/Attic/deform_support.c,v 1.1 1999-10-25 19:59:10 louis Exp $";
#endif

#include <internal_volume_io.h>		

private dummy_routine() {
  print ("Sorry, non-linear spatial registration not supported in\n");
  print ("this version of minctracc\n");

  return (OK);

}

