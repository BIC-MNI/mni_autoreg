INCLUDE_DIRECTORIES(Include)

ADD_LIBRARY(minctracc_files STATIC
   Files/read_data_files.c)


ADD_LIBRARY(minctracc_optimize STATIC
  Optimize/def_obj_functions.c 
  Optimize/extras.c 
  Optimize/sub_lattice.c
  Optimize/amoeba.c 
  Optimize/vox_space.c 
  Optimize/objectives.c 
  Optimize/optimize.c 
  Optimize/segment_table.c 
  Optimize/deform_support.c 
  Optimize/super_sample_def.c 
  Optimize/my_grid_support.c 
  Optimize/obj_fn_mutual_info.c 
  Optimize/do_nonlinear.c
)

# this file produces different results when compiled with optimization
SET_SOURCE_FILES_PROPERTIES(
Numerical/matrix_basics.c 
Numerical/rotmat_to_ang.c 

Numerical/cov_to_praxes.c 
Numerical/init_params.c 
Numerical/quaternion.c
Numerical/make_rots.c 
Numerical/default_def.c 
Numerical/quad_max_fit.c 
Numerical/stats.c
PROPERTIES COMPILE_FLAGS -O0)


ADD_LIBRARY(minctracc_numerical STATIC
  Numerical/cov_to_praxes.c 
  Numerical/init_params.c 
  Numerical/quaternion.c
  Numerical/make_rots.c 
  Numerical/matrix_basics.c 
  Numerical/rotmat_to_ang.c 
  Numerical/default_def.c 
  Numerical/quad_max_fit.c 
  Numerical/stats.c
)

ADD_LIBRARY(minctracc_volume  STATIC
  Volume/init_lattice.c 
  Volume/interpolation.c 
  Volume/volume_functions.c
)

ADD_EXECUTABLE(minctracc 
  Main/make_matlab_data_file.c 
  Main/minctracc.c 
)

TARGET_LINK_LIBRARIES(minctracc
 minctracc_files 
 minctracc_optimize 
 minctracc_numerical
 minctracc_volume
 Proglib)

ADD_EXECUTABLE(check_scale  Extra_progs/check_scale.c)

TARGET_LINK_LIBRARIES(check_scale
 minctracc_files 
 minctracc_optimize 
 minctracc_numerical
 minctracc_volume
 Proglib
)

ADD_EXECUTABLE(param2xfm    Extra_progs/param2xfm.c)

TARGET_LINK_LIBRARIES(param2xfm
 minctracc_files 
 minctracc_optimize 
 minctracc_numerical
 minctracc_volume
 Proglib
)

ADD_EXECUTABLE(rand_param    Extra_progs/rand_param.c)

TARGET_LINK_LIBRARIES(rand_param
 minctracc_files 
 minctracc_optimize 
 minctracc_numerical
 minctracc_volume
 Proglib
)

ADD_EXECUTABLE(volume_cog   Extra_progs/volume_cog.c)

TARGET_LINK_LIBRARIES(volume_cog
 minctracc_files 
 minctracc_optimize 
 minctracc_numerical
 minctracc_volume
 Proglib)

ADD_EXECUTABLE(xfm2param    Extra_progs/xfm2param.c)

TARGET_LINK_LIBRARIES(xfm2param
 minctracc_files 
 minctracc_optimize 
 minctracc_numerical
 minctracc_volume
 Proglib)

ADD_EXECUTABLE(zscore_vol    Extra_progs/zscore_vol.c)

TARGET_LINK_LIBRARIES(zscore_vol
 minctracc_files 
 minctracc_optimize 
 minctracc_numerical
 minctracc_volume
 Proglib)

ADD_EXECUTABLE(crispify     Extra_progs/crispify.c)
ADD_EXECUTABLE(xcorr_vol    Extra_progs/xcorr_vol.c)
ADD_EXECUTABLE(cmpxfm       Extra_progs/cmpxfm.c)
ADD_EXECUTABLE(make_lvv_vol Extra_progs/make_lvv_vol.c)

TARGET_LINK_LIBRARIES(make_lvv_vol
 Proglib
 )

INSTALL(TARGETS 
 minctracc 
 check_scale 
 crispify 
 param2xfm 
 volume_cog 
 rand_param 
 xcorr_vol 
 xfm2param 
 zscore_vol
 cmpxfm
 make_lvv_vol
 RUNTIME DESTINATION bin)
