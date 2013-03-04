#include <volume_io.h>
#include <Proglib.h>

#include "segment_table.h"

VIO_BOOL build_segment_table(Segment_Table **s_table, VIO_Volume d1, int groups)
{

  float
    frac;
  Data_types
    data_type;
  int
    i,min,max;

  Segment_Table *st;
  int *it;

  data_type = get_volume_data_type (d1);
  switch (data_type) {
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
                data_type);
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
    
    for(i=min; i<=max; i++) {                /* build the look up table */
      frac = 0.5 + ((float)(groups)-0.00001 )* (float)(i-min)/(float)(max-min);
      (*s_table)->table[i] = ROUND( frac );
    }
    
    return(TRUE);
  }
  else
    return(FALSE);
}


VIO_BOOL free_segment_table(Segment_Table *s_table)
{

  int *p;

  p = s_table->table;
  p += s_table->min;

  FREE(p);
  FREE(s_table);
  return(TRUE);

}

int     get_segment_LUT_value(int value, Segment_Table *s_table)
{

  if ((value>=s_table->min) && (value<=s_table->max))
    return(s_table->table[value]);
  else
    return(0);

}

