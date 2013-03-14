#include <volume_io.h>
#include "arg_data.h"

void initializeArgs(Arg_Data *args);
VIO_General_transform* minctracc( VIO_Volume source, VIO_Volume target, VIO_Volume sourceMask, VIO_Volume targetMask, VIO_General_transform *initialXFM, int iterations, Arg_Data *args);

