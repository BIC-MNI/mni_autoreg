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

ADD_LIBRARY(minctracc_numerical STATIC
	Numerical/cov_to_praxes.c 
	Numerical/init_params.c 
	Numerical/quaternion.c
	Numerical/make_rots.c 
	Numerical/matrix_basics.c 
	Numerical/rotmat_to_ang.c 
	Numerical/default_def.c 
	Numerical/quad_max_fit.c 
	Numerical/stats.c)

ADD_LIBRARY(minctracc_volume STATIC
	Volume/init_lattice.c 
	Volume/interpolation.c 
	Volume/volume_functions.c)


LINK_LIBRARIES(
 minctracc_files 
 minctracc_optimize 
 minctracc_numerical
 minctracc_volume
 Proglib
)

ADD_EXECUTABLE(minctracc 
  Main/make_matlab_data_file.c 
	Main/minctracc.c )

ADD_EXECUTABLE(check_scale Extra_progs/check_scale.c)
ADD_EXECUTABLE(crispify Extra_progs/crispify.c)
ADD_EXECUTABLE(param2xfm Extra_progs/param2xfm.c)
ADD_EXECUTABLE(volume_cog Extra_progs/volume_cog.c)
ADD_EXECUTABLE(rand_param Extra_progs/rand_param.c)
ADD_EXECUTABLE(xcorr_vol Extra_progs/xcorr_vol.c)
ADD_EXECUTABLE(xfm2param Extra_progs/xfm2param.c)
ADD_EXECUTABLE(zscore_vol Extra_progs/zscore_vol.c)

INSTALL(TARGETS 
 minctracc check_scale crispify param2xfm volume_cog rand_param xcorr_vol xfm2param zscore_vol
 RUNTIME DESTINATION bin)