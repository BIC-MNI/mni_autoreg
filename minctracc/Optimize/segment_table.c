#include <def_mni.h>
#include <recipes.h>

#include <def_segment_table.h>

public Boolean build_segment_table(Segment_Table **s_table, Volume d1, int groups)
{

  float
    frac;
  int
    i,min,max;

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
    print_error("Currently an unsupported data type (%d).", __FILE__, __LINE__, 
		d1->data_type,0,0,0,0);
    return(FALSE);
  }

  ALLOC( *s_table, 1);

  if (*s_table != NULL) {
    (*s_table)->min = min;
    (*s_table)->max = max;
    (*s_table)->groups = groups;
    (*s_table)->segment = get_segment_LUT_value;
    
    (*s_table)->table = ivector(min,max);	/* allocate the table space! */
    
    for_inclusive(i, min, max) {	        /* build the look up table */
      frac = 0.5 + ((float)(groups)-0.00001 )* (float)(i-min)/(float)(max-min);
      (*s_table)->table[i] = ROUND( frac );
    }
    
    return(TRUE);
  }
  else
    return(FALSE);
}


public Boolean free_segment_table(Segment_Table *s_table)
{

  free_ivector(s_table->table, s_table->min, s_table->max);
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

