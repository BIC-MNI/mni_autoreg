/* ----------------------------- MNI Header -----------------------------------
@NAME       : switch_obj_func.c
@INPUT      : *a1, s1, s3 and sample
@OUTPUT     : updated values for s1 and s3, updated pointed for *a1
@DESCRIPTION: this is a case statement that is to be included within the 
              procedure go_get_samples_with_offset(), inside the loop over
	      all nodes of the sublattice, just after interpolation of the 
	      value for 'sample'
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

@CREATED    : 
@MODIFIED   : $Log: switch_obj_func.c,v $
@MODIFIED   : Revision 96.2  1997-11-12 21:07:43  louis
@MODIFIED   : added support for chamfer distance local objective function
@MODIFIED   :
 * Revision 96.1  1997/11/03  20:05:41  louis
 * no changes
 *
 * Revision 96.0  1996/08/21  18:22:10  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 96.0  1996/08/21  18:22:10  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:22:06  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:16:03  louis
 * Never released version 0.95
 *
 * Revision 1.2  1996/08/12  14:15:59  louis
 * Pre-release
 *
 * Revision 1.1  1996/03/25  10:33:15  collins
 * Initial revision
 *
---------------------------------------------------------------------------- */

	switch (obj_func) {
	case NONLIN_XCORR:
	  s1 += *a1++ * sample; /* compute correlation */
	  s3 += sample * sample;
	  break;
	case NONLIN_DIFF:
	  tmp = *a1++ - sample;
	  if (tmp<0) tmp *= -1.0;
	  s1 += tmp;            /* add the difference */
	  break;
	case NONLIN_LABEL:
	  tmp = *a1++ - sample;
	  if (tmp<0) tmp *= -1.0;
	  if (tmp < 0.01)
	    s1 += 1.0;          /* count up similar labels */
	  break;
	case NONLIN_CHAMFER:
	  if (*a1++ > 0) {
             s1 += sample;         /* add the distance */
             number_of_nonzero_samples++;
          }
	  break;
	default:
	  print_error_and_line_num("Objective function %d not supported in go_get_samples_with_offset",__FILE__, __LINE__,obj_func);
	}
