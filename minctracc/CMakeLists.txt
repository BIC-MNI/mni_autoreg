
INCLUDE_DIRECTORIES(Include)

IF(LIBLBFGS_FOUND)
  INCLUDE_DIRECTORIES(Include ${LIBLBFGS_INCLUDE_DIR})
  LINK_DIRECTORIES(${LIBLBFGS_LIBRARY_DIR})
ENDIF(LIBLBFGS_FOUND)

SET ( MINCTRACC_FILES
  Files/read_data_files.c
)

SET ( MINCTRACC_OPTIMIZE
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

SET (MINCTRACC_NUMERICAL
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

SET (MINCTRACC_VOLUME
  Volume/init_lattice.c 
  Volume/interpolation.c 
  Volume/volume_functions.c
)

SET (MINCTRACC_PROGLIB
  ../Proglib/get_history.c
  ../Proglib/print_error.c
  ../Proglib/print_version.c
)

SET (MINCTRACC_MAIN
  Main/minctracclib.c
  Main/make_matlab_data_file.c
)

SET ( LIB_MINCTRACC_HEADERS
  Include/minctracc_point_vector.h
  Include/minctracc_arg_data.h
  Include/libminctracc.h
)


SET ( MINCTRACC_HEADERS
  Include/libminctracc.h
  Include/amoeba.h
  Include/constants.h
  Include/cov_to_praxes.h
  Include/deform_support.h
  Include/extras.h
  Include/globals.h
  Include/init_lattice.h
  Include/interpolation.h
  Include/local_macros.h
  Include/make_rots.h
  Include/matrix_basics.h
  Include/minctracc.h
  Include/objectives.h
  Include/quad_max_fit.h
  Include/quaternion.h
  Include/rotmat_to_ang.h
  Include/segment_table.h
  Include/stats.h
  Include/sub_lattice.h
  Include/super_sample_def.h
  Include/vox_space.h
  ../Proglib/Proglib.h
  ${LIB_MINCTRACC_HEADERS}
)

#ADD_LIBRARY(minctracc_files STATIC ${MINCTRACC_FILES})

#ADD_LIBRARY(minctracc_optimize STATIC ${MINCTRACC_OPTIMIZE})

#ADD_LIBRARY(minctracc_numerical STATIC ${MINCTRACC_NUMERICAL})

#ADD_LIBRARY(minctracc_volume  STATIC ${MINCTRACC_VOLUME})

ADD_LIBRARY(_minctracc ${LIBRARY_TYPE} 
  ${MINCTRACC_HEADERS} 
  ${MINCTRACC_FILES} 
  ${MINCTRACC_OPTIMIZE} 
  ${MINCTRACC_NUMERICAL} 
  ${MINCTRACC_VOLUME}
  ${MINCTRACC_PROGLIB}
  ${MINCTRACC_MAIN}
)

TARGET_LINK_LIBRARIES(
  _minctracc
  ${VOLUME_IO_LIBRARIES}
  ${LIBMINC_LIBRARIES}
)

IF(LIBLBFGS_FOUND)
  TARGET_LINK_LIBRARIES(
    _minctracc
    ${LIBLBFGS_LIBRARY}
  )
ENDIF(LIBLBFGS_FOUND)

ADD_EXECUTABLE(check_scale  Extra_progs/check_scale.c)

TARGET_LINK_LIBRARIES(check_scale
  _minctracc
)

ADD_EXECUTABLE(minctracc 
  Main/make_matlab_data_file.c 
  Main/minctracc.c 
)

TARGET_LINK_LIBRARIES(
  minctracc
  _minctracc
)

ADD_EXECUTABLE(param2xfm    Extra_progs/param2xfm.c)

TARGET_LINK_LIBRARIES(param2xfm
 _minctracc 
)


# ADD_EXECUTABLE(rand_param    Extra_progs/rand_param.c)
# 
# TARGET_LINK_LIBRARIES(rand_param
#  minctracc_files 
#  minctracc_optimize 
#  minctracc_numerical
#  minctracc_volume
#  Proglib
# )
# 
ADD_EXECUTABLE(volume_cog   Extra_progs/volume_cog.c)
 
TARGET_LINK_LIBRARIES(volume_cog
  _minctracc )
 
ADD_EXECUTABLE(xfm2param    Extra_progs/xfm2param.c)
 
TARGET_LINK_LIBRARIES(xfm2param
  _minctracc
  )
 
# ADD_EXECUTABLE(zscore_vol    Extra_progs/zscore_vol.c)
# 
# TARGET_LINK_LIBRARIES(zscore_vol
#  minctracc_files 
#  minctracc_optimize 
#  minctracc_numerical
#  minctracc_volume
#  Proglib)
# 
ADD_EXECUTABLE(crispify     Extra_progs/crispify.c)
ADD_EXECUTABLE(xcorr_vol    Extra_progs/xcorr_vol.c)
ADD_EXECUTABLE(cmpxfm       Extra_progs/cmpxfm.c)

# ADD_EXECUTABLE(make_lvv_vol Extra_progs/make_lvv_vol.c)
# 
# TARGET_LINK_LIBRARIES(make_lvv_vol
#  Proglib
#  )

INSTALL(TARGETS _minctracc ${LIBRARY_INSTALL} DESTINATION lib)
INSTALL(FILES   ${LIB_MINCTRACC_HEADERS} DESTINATION include)
INSTALL(TARGETS minctracc RUNTIME DESTINATION bin) 

INSTALL(TARGETS 
# minctracc-e 
 check_scale 
 crispify 
 param2xfm 
 volume_cog 
#  rand_param 
 xcorr_vol 
 xfm2param 
#  zscore_vol
  cmpxfm
#  make_lvv_vol
 RUNTIME DESTINATION bin)


INSTALL_MAN_PAGES( ${CMAKE_INSTALL_PREFIX}/man
  minctracc.man1
)
