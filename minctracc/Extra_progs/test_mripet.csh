#!/bin/csh

minctracc pcvelican_025943_46.mnc velican_marc_hres.mnc -est_center -est_translation \
	-vr -thresh 20 -tol 0.01 out.xfm -step 8 8 8 -groups 20 -debug -lsq6 -clobber

xfminvert out.xfm inv.xfm

mincresample velican_marc_hres.mnc tmp.mnc -like pcvelican_025943_46.mnc -transformation inv.xfm

register tmp.mnc pcvelican_025943_46.mnc


