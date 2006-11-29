/* ----------------------------- MNI Header -----------------------------------
@NAME       : segment_table.h
@DESCRIPTION: structures and prototypes for segment table manipulation
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

@CREATED    : Fri Jun 18 11:03:43 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */

typedef struct Segment_Table_Struct Segment_Table;

                                /* this function will return the group value */
typedef int (*Segment_Function)(int value, Segment_Table *table);


struct Segment_Table_Struct
{
  int min;                        /* minimum voxel value */
  int max;                        /* maximum voxel value */
  int groups;                        /* number of groups in table */
  int *table;                        /* list of look up values for segmentation */
  Segment_Function segment;     /* name of the function used to apply the segmentation */
};








/* ------------------------ prototypes for segment table manipulation ---------- */

VIO_BOOL build_segment_table(Segment_Table **table, VIO_Volume d1, int groups);

VIO_BOOL free_segment_table(Segment_Table *table);

int     get_segment_LUT_value(int value, Segment_Table *table);


