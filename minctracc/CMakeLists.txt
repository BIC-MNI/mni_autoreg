INCLUDE_DIRECTORIES(Include ${LIBLBFGS_INCLUDE_DIR})

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

SET ( MINCTRACC_HEADERS
	Include/libminctracc.h
	Include/amoeba.h
	Include/arg_data.h
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
	Include/point_vector.h
	Include/quad_max_fit.h
	Include/quaternion.h
	Include/rotmat_to_ang.h
	Include/segment_table.h
	Include/stats.h
	Include/sub_lattice.h
	Include/super_sample_def.h
	Include/vox_space.h
	../Proglib/Proglib.h
)


ADD_LIBRARY(minctracc_files STATIC ${MINCTRACC_FILES})

ADD_LIBRARY(minctracc_optimize STATIC ${MINCTRACC_OPTIMIZE})

ADD_LIBRARY(minctracc_numerical STATIC ${MINCTRACC_NUMERICAL})

ADD_LIBRARY(minctracc_volume  STATIC ${MINCTRACC_VOLUME})

ADD_LIBRARY(minctracc ${LIBRARY_TYPE} 
	${MINCTRACC_HEADERS} 
	${MINCTRACC_FILES} 
	${MINCTRACC_OPTIMIZE} 
	${MINCTRACC_NUMERICAL} 
	${MINCTRACC_VOLUME}
	${MINCTRACC_PROGLIB}
	${MINCTRACC_MAIN}
)


LINK_LIBRARIES(
 minctracc_files 
 minctracc_optimize 
 minctracc_numerical
 minctracc_volume
 Proglib
)

TARGET_LINK_LIBRARIES(
	minctracc
	${VOLUME_IO_LIBRARIES}
	${LIBMINC_LIBRARIES}
	${LIBLBFGS_LIBRARIES}
)

#ADD_EXECUTABLE(minctracc-e 
 # Main/make_matlab_data_file.c 
#	Main/minctracc.c )

ADD_EXECUTABLE(check_scale Extra_progs/check_scale.c)
ADD_EXECUTABLE(crispify    Extra_progs/crispify.c)
ADD_EXECUTABLE(param2xfm   Extra_progs/param2xfm.c)
ADD_EXECUTABLE(volume_cog  Extra_progs/volume_cog.c)
ADD_EXECUTABLE(rand_param  Extra_progs/rand_param.c)
ADD_EXECUTABLE(xcorr_vol   Extra_progs/xcorr_vol.c)
ADD_EXECUTABLE(xfm2param   Extra_progs/xfm2param.c)
ADD_EXECUTABLE(zscore_vol  Extra_progs/zscore_vol.c)
ADD_EXECUTABLE(cmpxfm      Extra_progs/cmpxfm.c)

INSTALL(TARGETS minctracc ${LIBRARY_INSTALL} DESTINATION lib)
INSTALL(FILES ${MINCTRACC_HEADERS} DESTINATION include)

INSTALL(TARGETS 
# minctracc-e 
 check_scale 
 crispify 
 param2xfm 
 volume_cog 
 rand_param 
 xcorr_vol 
 xfm2param 
 zscore_vol
 cmpxfm
 RUNTIME DESTINATION bin)
