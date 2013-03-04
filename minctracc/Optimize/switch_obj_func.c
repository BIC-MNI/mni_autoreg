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
@MODIFIED   : Revision 96.7  2009-04-03 18:36:59  louis
@MODIFIED   : made changes to use only DOUBLES for input source and model volumes, and for all estimation of deformation fields
@MODIFIED   :
@MODIFIED   : Revision 96.6  2006/11/29 09:09:34  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 96.5  2005/07/18 19:14:02  rotor
@MODIFIED   :  * Optimisations to code resulting in 30% speed increase for nonlinear fitting
@MODIFIED   :
@MODIFIED   : Revision 96.4  2005/06/28 18:56:18  rotor
@MODIFIED   :  * added masking for feature volumes (irina and patricia)
@MODIFIED   :
@MODIFIED   : Revision 96.3  2003/02/04 06:08:46  stever
@MODIFIED   : Add support for correlation coefficient and sum-of-squared difference.
@MODIFIED   :
@MODIFIED   : Revision 96.2  1997/11/12 21:07:43  louis
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
   
 case NONLIN_CORRCOEFF:		/* compute the correlation coefficient (placed first, because used most often) */
   s1 += *a1;
   s2 += sample;
   s3 += (*a1) * (*a1);
   s4 += sample * sample;
   s5 += *a1 * sample;
   number_of_nonzero_samples++;
   break;
   
 case NONLIN_XCORR:          /* compute correlation */
   s2 += (*a1) * (*a1);
   s1 += *a1 * sample; 
   s3 += sample * sample;
   break;
   
 case NONLIN_CHAMFER:
   if (*a1 > 0) {
     s1 += sample;         /* add the distance */
     number_of_nonzero_samples++;
   }
   break;
   
 case NONLIN_SQDIFF:		/* sompute the squared intensity difference */
   tmp = *a1 - sample;
   s1 += tmp*tmp;
   number_of_nonzero_samples++;
   break;
   
 case NONLIN_DIFF:           /* compute the sample-to-sample difference */
   tmp = *a1 - sample;
   if (tmp<0){
     tmp *= -1.0;
   }
   s1 += tmp;            
   number_of_nonzero_samples++;
   break;
   
 case NONLIN_LABEL:		/* sum up the number of similar labels */
   tmp = *a1 - sample;
   if (tmp<0){
     tmp *= -1.0;
   }
   if (tmp < 0.01){
     s1 += 1.0;          /* count similar labels */
   }
   number_of_nonzero_samples++;
   break;
   
 default:
   print_error_and_line_num("Objective function %d not supported in go_get_samples_with_offset",__FILE__, __LINE__,obj_func);
 }
