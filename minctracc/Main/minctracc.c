/* ----------------------------- MNI Header -----------------------------------
   @NAME       : minctracc.c
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

   @CREATED    : February 3, 1992 - louis collins
   @MODIFIED   : $Log: minctracc.c,v $
   @MODIFIED   : Revision 96.24  2011-02-24 20:02:35  louis
   @MODIFIED   : update for normalized mutual informatoin
   @MODIFIED   :
   @MODIFIED   : Revision 96.23  2010-04-01 04:49:16  rotor
   @MODIFIED   :  * fixed time.h include
   @MODIFIED   :
   @MODIFIED   : Revision 96.22  2009/05/26 18:03:07  claude
   @MODIFIED   : free more memory after usage
   @MODIFIED   :
   @MODIFIED   : Revision 96.21  2009/05/22 15:49:19  claude
   @MODIFIED   : fixed memory bug freeing initial transform
   @MODIFIED   :
   @MODIFIED   : Revision 96.20  2009/04/03 18:36:59  louis
   @MODIFIED   : made changes to use only DOUBLES for input source and model volumes, and for all estimation of deformation fields
   @MODIFIED   :
   @MODIFIED   : Revision 96.19  2009/03/13 19:51:31  claude
   @MODIFIED   : fixed bug in offsets for minctracc and free memory upon exit
   @MODIFIED   :
   @MODIFIED   : Revision 96.18  2008/10/08 15:17:49  louis
   @MODIFIED   : added -nmi option for linear normalized mutual information
   @MODIFIED   :
   @MODIFIED   : Revision 96.17  2006/11/30 09:07:31  rotor
   @MODIFIED   :  * many more changes for clean minc 2.0 build
   @MODIFIED   :
   @MODIFIED   : Revision 96.15  2006/06/04 07:02:35  rotor
   @MODIFIED   :  * Fixed 64 bit function pointer and ParseArgv problem with an enum for
   @MODIFIED   :       objective function type and interpolation type (thanks jason)
   @MODIFIED   :
   @MODIFIED   : Revision 96.14  2005/07/20 20:45:48  rotor
   @MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
   @MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
   @MODIFIED   :     * Removed old VOLUME_IO cruft #defines
   @MODIFIED   :     * Fixed up all Makefile.am's in subdirs
   @MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
   @MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
   @MODIFIED   :
   @MODIFIED   : Revision 96.13  2005/06/28 18:56:14  rotor
   @MODIFIED   :  * added masking for feature volumes (irina and patricia)
   @MODIFIED   :
   @MODIFIED   : Revision 96.12  2004/03/18 06:51:03  rotor
   @MODIFIED   :  * changed make_model from csh to sh
   @MODIFIED   :  * removed an extraneous printf from minctracc
   @MODIFIED   :
   @MODIFIED   : Revision 96.11  2004/02/12 05:54:22  rotor
   @MODIFIED   :  * removed /static defs
   @MODIFIED   :
   @MODIFIED   : Revision 96.10  2004/02/04 20:42:33  lenezet
   @MODIFIED   : *** empty log message ***
   @MODIFIED   :
   @MODIFIED   : Revision 96.9  2003/02/04 06:08:44  stever
   @MODIFIED   : Add support for correlation coefficient and sum-of-squared difference.
   @MODIFIED   :
   @MODIFIED   : Revision 96.8  2002/12/13 21:16:11  lenezet
   @MODIFIED   : nonlinear in 2D has changed. The option -2D-non-lin is no more necessary. The grid transform has been adapted to feet on the target volume whatever is size. The Optimization is done on the dimensions for which "count" is greater than 1.
   @MODIFIED   :
   @MODIFIED   : Revision 96.7  2002/11/20 21:38:31  lenezet
   @MODIFIED   :
   @MODIFIED   : Fix the code to take in consideration the direction cosines especially in the grid transform.
   @MODIFIED   : Add an option to choose the maximum expected deformation magnitude.
   @MODIFIED   :
   @MODIFIED   : Revision 96.6  2002/08/14 19:54:42  lenezet
   @MODIFIED   :  quaternion option added for the rotation
   @MODIFIED   :
   @MODIFIED   : Revision 96.5  2002/03/26 14:15:38  stever
   @MODIFIED   : Update includes to <volume_io/foo.h> style.
   @MODIFIED   :
   @MODIFIED   : Revision 96.4  2002/03/07 19:08:34  louis
   @MODIFIED   : Added -lattice_diameter as an optionto minctracc to account for a
   @MODIFIED   : problem with the automated calculation of the sub-lattice diameter.
   @MODIFIED   : It used to be step*3*2 - which was pretty big, when step = 8mm.
   @MODIFIED   :
   @MODIFIED   : Now, the sub lattice diameter can be input on the command line, and I
   @MODIFIED   : suggest a lattice size 3 times greater than the step size.
   @MODIFIED   :
   @MODIFIED   : If not on the command line, the default is = 24mm.
   @MODIFIED   :
   @MODIFIED   : Revision 96.3  2000/03/15 08:42:41  stever
   @MODIFIED   : Code cleanup: all functions prototyped (except ParseArgs.c), no useless declarations, etc
   @MODIFIED   :
   @MODIFIED   : Revision 96.2  2000/02/20 04:01:03  stever
   @MODIFIED   : * use new history_string() function to generate history strings
   @MODIFIED   :   when outputting MNI files (.mnc, .xfm)
   @MODIFIED   : * removed unused vax routines from Proglib
   @MODIFIED   : * tuned configure script; CPPFLAGS and LDFLAGS are now left alone,
   @MODIFIED   :   for the installer to use
   @MODIFIED   :
   @MODIFIED   : Revision 96.1  1999/10/25 19:52:19  louis
   @MODIFIED   : final checkin before switch to CVS
   @MODIFIED   :
 * Revision 96.0  1996/08/21  18:21:51  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:21:49  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:15:42  louis
 * Never released version 0.95
 *
 * Revision 1.16  1996/08/12  14:15:40  louis
 * Pre-release
 *
 * Revision 1.15  1995/09/28  13:24:26  collins
 * added multiple feature volume loading with get_feature_volumes()
 * to be able to correlate multiple features at the same time.
 *
 * Revision 1.14  1995/09/28  11:54:52  collins
 * working version, just prior to release 0.9 of mni_autoreg
 *
 * Revision 1.13  1995/02/22  08:56:06  collins
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.12  94/05/28  16:18:54  louis
 * working version before modification of non-linear optimiation
 * 
 * Revision 1.11  94/04/26  12:54:30  louis
 * updated with new versions of make_rots, extract2_parameters_from_matrix 
 * that include proper interpretation of skew.
 * 
 * Revision 1.10  94/04/06  11:48:43  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.9  94/02/21  16:35:51  louis
 * version before feb 22 changes
 * 
 * Revision 1.8  93/11/15  16:27:06  louis
 * working version, with new library, with RCS revision stuff,
 * before deformations included
 * 
 * Revision 1.7  93/11/15  13:12:10  louis
 * working version, before deform installation
 * 

Thu May 20 11:21:22 EST 1993 lc
             complete re-write to use minc library

Wed May 26 13:05:44 EST 1993 lc

        complete re-write of fit_vol to:
        - use minc library (libminc.a) to read and write .mnc files
        - use libmni.a
        - use dave's volume structures.
        - read and write .xfm files, for the transformations
        - allow different interpolation schemes on voxel values
        - allow different objective functions: xcorr, ssc, var ratio
        - allow mask volumes on both source and model
        - calculate the bounding box of both source and model,
          to map samples from smallest volume into larger one (for speed)

---------------------------------------------------------------------------- */

#ifndef lint
static char minctracc_rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/minctracc/Main/minctracc.c,v 96.24 2011-02-24 20:02:35 louis Exp $";
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /*HAVE_CONFIG_H*/

#include <float.h>
#include <volume_io.h>
#include <minctracc.h>
#include <globals.h>
#include <objectives.h>
#include "local_macros.h"

/*************************************************************************/
int main ( int argc, char* argv[] )
{
	return minctraccOldFashioned(argc,argv);
}

