#include <volume_io.h>
#include <print_error.h>

#include "segment_table.h"

public BOOLEAN build_segment_table(Segment_Table **s_table, Volume d1, int groups)
{

  float
    frac;
  int
    i,min,max;

  Segment_Table *st;
  int *it;

  switch (d1->data_type) {
  case   UNSIGNED_BYTE:
    min = 0; max = (1<<8)-1;
    break;
  case  SIGNED_BYTE:
    min = -(1<<7); max = (1<<7)-1;
    break;
  case  UNSIGNED_SHORT:
    min = 0; max = (1<<16)-1;
    break;
  case  SIGNED_SHORT:
    min = -(1<<15); max = (1<<15)-1;
    break;
  default:
    print_error_and_line_num("Currently an unsupported data type (%d).", __FILE__, __LINE__, 
		d1->data_type);
    return(FALSE);
  }

  ALLOC( st, 1);

  *s_table = st;
  if (*s_table != NULL) {
    (*s_table)->min = min;
    (*s_table)->max = max;
    (*s_table)->groups = groups;
    (*s_table)->segment = get_segment_LUT_value;
    

    ALLOC(it,max-min+1);      /* allocate the table space! */
    it -= min;

    (*s_table)->table = it;
    
    for_inclusive(i, min, max) {	        /* build the look up table */
      frac = 0.5 + ((float)(groups)-0.00001 )* (float)(i-min)/(float)(max-min);
      (*s_table)->table[i] = ROUND( frac );
    }
    
    return(TRUE);
  }
  else
    return(FALSE);
}


public BOOLEAN free_segment_table(Segment_Table *s_table)
{

  int *p;

  p = s_table->table;
  p += s_table->min;

  FREE(p);
  FREE(s_table);
  return(TRUE);

}

public int     get_segment_LUT_value(int value, Segment_Table *s_table)
{

  if ((value>=s_table->min) && (value<=s_table->max))
    return(s_table->table[value]);
  else
    return(0);

}

