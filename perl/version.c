/*
 * Program to write out the build information.
 */

#include "config.h"

int main (void)
{
   printf ("VERSION=%s\n", MNI_AUTOREG_VERSION);
   printf ("COMPILE_TIME=%s\n", MNI_AUTOREG_COMPILE_TIME);
   printf ("COMPILE_DATE=%s\n", MNI_AUTOREG_COMPILE_DATE);
   printf ("COMPILE_USER=%s\n", MNI_AUTOREG_COMPILE_USER);
   printf ("COMPILE_HOST=%s\n", MNI_AUTOREG_COMPILE_HOST);
   printf ("COMPILE_SYSTEM=%s\n", MNI_AUTOREG_COMPILE_SYSTEM);
   printf ("LONG_VERSION=%s\n", MNI_AUTOREG_LONG_VERSION);
   exit (0);
}
