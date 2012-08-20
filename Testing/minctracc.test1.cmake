#! /bin/sh
set -e

minctracc -identity object1_dxyz.mnc object2_dxyz.mnc \
     -est_center -debug -simplex 10 -lsq6 -step 8 8 8 \
     -clobber output.test1.xfm

     
param2xfm -rotation -4 7 10 -translation  5 2 -6 -clobber ideal.test1.xfm

#TODO: check why tolerance have to be 0.05 
if ! cmpxfm -linear_tolerance 0.05 -translation_tolerance 0.05 output.test1.xfm ideal.test1.xfm; then
  echo >&2 $0 failed: minctracc produced incorrect results.
  exit 1
fi
