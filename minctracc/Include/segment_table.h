/* ----------------------------- MNI Header -----------------------------------
@NAME       : segment_table.h
@DESCRIPTION: structures and prototypes for segment table manipulation
@CREATED    : Fri Jun 18 11:03:43 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */

typedef struct Segment_Table_Struct Segment_Table;

				/* this function will return the group value */
typedef int (*Segment_Function)(int value, Segment_Table *table);


struct Segment_Table_Struct
{
  int min;			/* minimum voxel value */
  int max;			/* maximum voxel value */
  int groups;			/* number of groups in table */
  int *table;			/* list of look up values for segmentation */
  Segment_Function segment;     /* name of the function used to apply the segmentation */
};








/* ------------------------ prototypes for segment table manipulation ---------- */

public BOOLEAN build_segment_table(Segment_Table **table, Volume d1, int groups);

public BOOLEAN free_segment_table(Segment_Table *table);

public int     get_segment_LUT_value(int value, Segment_Table *table);


