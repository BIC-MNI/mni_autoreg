# make -k 'LIB_DIRS=-L/nil/louis/IRIX_6 -L/usr/local/lib' 'MACHINE_LIB=-lm' 'OPT=-g -DNONLINEAR_RELEASE' animal
# --------------------------------------------------------------------
# 

# Makefile for minctracc
#

DAVE_LIB = volume_io
include ../Makefile.include

EXTRA_PROGS = \
	Extra_progs/avg_voxel_vol.c \
	Extra_progs/change_range.c \
	Extra_progs/load_vol.c \
	Extra_progs/check_scale.c \
	Extra_progs/overlap_corr.c \
	Extra_progs/param2xfm.c \
	Extra_progs/rand_param.c \
	Extra_progs/xcorr_vol.c \
	Extra_progs/xfm2param.c \
	Extra_progs/xfminvert.c \
	Extra_progs/def_to_vols.c \
	Extra_progs/def_to_lines.c \
	Extra_progs/deform_object.c \
	Extra_progs/mst_to_lines.c \
	Extra_progs/find_tubes.c \
	Extra_progs/crispify.c \
	Extra_progs/test_cov.c \
	Extra_progs/test_cov_pr.c \
	Extra_progs/test_quad.c \
	Extra_progs/test_kernels.c \
	Extra_progs/test_prindir.c \
	Extra_progs/test_eigen.c \
	Extra_progs/test_lines.c \
	Extra_progs/test_deriv.c \
	Extra_progs/reversedef.c \
	Extra_progs/inbetween_def.c \
	Extra_progs/random_def.c \
	Extra_progs/random_vol.c \
	Extra_progs/randomize_volume.c \
	Extra_progs/volume_cog.c \
	Extra_progs/volume_stat.c \
	Extra_progs/zscore_vol.c \
	Extra_progs/make_zscore_vol.c

  PROG = minctracc

MANPAGES =minctracc.1

  SRCS =  \
	Files/read_data_files.c \
	Main/make_matlab_data_file.c \
	Main/minctracc.c \
	Numerical/cov_to_praxes.c \
	Numerical/init_params.c \
	Numerical/make_rots.c \
	Numerical/matrix_basics.c \
	Numerical/rotmat_to_ang.c \
	Numerical/default_def.c \
	Numerical/quad_max_fit.c \
	Numerical/stats.c \
	Optimize/def_obj_functions.c \
	Optimize/extras.c \
	Optimize/sub_lattice.c \
	Optimize/amoeba.c \
	Optimize/vox_space.c \
	Optimize/objectives.c \
	Optimize/optimize.c \
	Optimize/segment_table.c \
	Optimize/deform_support.c \
	Optimize/super_sample_def.c \
	Optimize/my_grid_support.c \
	Optimize/obj_fn_mutual_info.c \
	Volume/init_lattice.c \
	Volume/interpolation.c \
	Volume/volume_functions.c \
	Optimize/do_nonlinear.c 


 HEADER_SRCS =  Include/boolean.h \
		Include/cov_to_praxes.h \
		Include/def_segment_table.h \
		Include/def_tps.h \
		Include/globals.h \
		Include/make_rots.h \
		Include/matrix_basics.h \
		Include/minctracc.h \
		Include/def_tag_io.h \
		Include/rotmat_to_ang.h 

  OBJS = $(SRCS:.c=.o)

  LINT_LIST1 = 	$(OBJS:.o=.ln)     # list of lint files in program 1

   ALL_PROGS = $(PROG) volume_cog test_cov check_scale param2xfm xfm2param

   ALL_LINTS = lint_$(PROG) 

     OPT = -g3 -DNONLINEAR_RELEASE

$(PROG): $(OBJS) 
	$(CC) $(OPT) $(OBJS) -o $@ $(LDOPT)


animal: $(OBJS) 
	$(CC) $(OPT) $(OBJS) -o $@ $(LDOPT) 

tags: $(SRCS) $(HEADER_SRCS)
	etags  $(SRCS) $(HEADER_SRCS) \
	/nil/louis/Release/Volume_io/*/*.c \
	/nil/louis/Release/Volume_io/Include/*.h 

lint: lint_$(PROG)

lint_$(PROG): $(LINT_LIST1)
	$(LINT) $(LINTFLAGS) $(LINT_LIST1) $(LINKOPT)

DEF_TO_VOLS_OBJS = Extra_progs/def_to_vols.o
def_to_vols: $(DEF_TO_VOLS_OBJS)
	$(CC) $(OPT) $(DEF_TO_VOLS_OBJS) -o $@ $(LDOPT)

DEF_TO_SCALE_OBJS = Extra_progs/def_to_scale.o
def_to_scale: $(DEF_TO_SCALE_OBJS)
	$(CC) $(OPT) $(DEF_TO_SCALE_OBJS) -o $@ $(LDOPT)

DEF_TO_MAG_OBJS = Extra_progs/def_to_mag.o
def_to_mag: $(DEF_TO_MAG_OBJS)
	$(CC) $(OPT) $(DEF_TO_MAG_OBJS) -o $@ $(LDOPT)

TEST_QUAD_OBJS = Extra_progs/test_quad.o Numerical/quad_max_fit.o
test_quad: $(TEST_QUAD_OBJS)
	$(CC) $(OPT) $(TEST_QUAD_OBJS) -o $@ $(LDOPT)

TEST_PRINDIR_OBJS = Extra_progs/test_prindir.o Numerical/quad_max_fit.o \
	Numerical/cov_to_praxes.o
test_prindir: $(TEST_PRINDIR_OBJS)
	$(CC) $(OPT) $(TEST_PRINDIR_OBJS) -o $@ $(LDOPT)

TEST_EIGEN_OBJS = Extra_progs/test_eigen.o Numerical/quad_max_fit.o \
	Numerical/cov_to_praxes.o
test_eigen: $(TEST_EIGEN_OBJS)
	$(CC) $(OPT) $(TEST_EIGEN_OBJS) -o $@ $(LDOPT)

TEST_KERNEL = Extra_progs/test_kernels.o 
test_kernels: $(TEST_KERNEL)
	$(CC) $(OPT) $(TEST_KERNEL) -o $@ $(LDOPT)

CORRECT_SD_VOL = Extra_progs/correct_sd_vols.o 
correct_sd_vol: $(CORRECT_SD_VOL)
	$(CC) $(OPT) $(CORRECT_SD_VOL) -o $@ $(LDOPT)

MAKE_LVV_VOL_OBJS = Extra_progs/make_lvv_vol.o 
make_lvv_vol: $(MAKE_LVV_VOL_OBJS)
	$(CC) $(OPT) $(MAKE_LVV_VOL_OBJS) -o $@ $(LDOPT)

CRISP_OBJS = Extra_progs/crispify.o 
crispify: $(CRISP_OBJS)
	$(CC) $(OPT) $(CRISP_OBJS) -o $@ $(LDOPT)

DEF_LINES_OBJS = Extra_progs/def_to_lines.o 
def_to_lines: $(DEF_LINES_OBJS)
	$(CC) $(OPT) $(DEF_LINES_OBJS) -o $@ -lbicpl $(LDOPT) 

DEFORM_OBJS = Extra_progs/deform_object.o 
deform_object: $(DEFORM_OBJS)
	$(CC) $(OPT) $(DEFORM_OBJS) -o $@ -lbicpl $(LDOPT) 

SCAN_OBJECTS_TO_VOLUME_O = Extra_progs/scan_objects_to_volume.o
Extra_progs/scan_objects_to_volume.o:Extra_progs/scan_objects_to_volume.c
	$(CC) $(OPT) -c Extra_progs/scan_objects_to_volume.c  -o $@ \
	      -I/nil/louis/minc/myprogs -IInclude -I../Proglib \
	      -I/nil/david/Release/BIC_PL/Include \
	      -I/nil/david/Release/Volume_io/Include \
	      -I/usr/local/include -Dpublic= -Dprivate=static -DHAVE_CONFIG_H -c 
scan_objects_to_volume: $(SCAN_OBJECTS_TO_VOLUME_O)
	$(CC) $(OPT) $(SCAN_OBJECTS_TO_VOLUME_O) -o $@ -lbicpl -L../Proglib \
		-L/nil/david/Release/BIC_PL -L/nil/david/Release/Volume_io \
		-L/usr/local/lib -L/nil/david/Release/BIC_PL/IRIX_5 \
		-lmincprog -lvolume_io -lminc -lnetcdf -lc_s  -lm

RMS_DIST_LINES_O = Extra_progs/rms_dist_lines.o
Extra_progs/rms_dist_lines.o:Extra_progs/rms_dist_lines.c
	$(CC) $(OPT) -c Extra_progs/rms_dist_lines.c  -o $@ \
	      -I/nil/louis/minc/myprogs -IInclude -I../Proglib \
	      -I/nil/david/Release/BIC_PL/Include \
	      -I/nil/david/Release/Volume_io/Include \
	      -I/usr/local/include -Dpublic= -Dprivate=static -DHAVE_CONFIG_H -c 
rms_dist_lines: $(RMS_DIST_LINES_O)
	$(CC) $(OPT) $(RMS_DIST_LINES_O) -o $@ -lbicpl -L../Proglib \
		-L/nil/david/Release/BIC_PL -L/nil/david/Release/Volume_io \
		-L/usr/local/lib -L/nil/david/Release/BIC_PL/IRIX_5 \
		-lmincprog -lvolume_io -lminc -lnetcdf -lc_s  -lm

RMS_DIST_OBJECTS_O = Extra_progs/rms_dist_objects.o
Extra_progs/rms_dist_objects.o:Extra_progs/rms_dist_objects.c
	$(CC) $(OPT) -c Extra_progs/rms_dist_objects.c  -o $@ \
	      -I/nil/louis/minc/myprogs -IInclude -I../Proglib \
	      -I/nil/david/Release/BIC_PL/Include \
	      -I/nil/david/Release/Volume_io/Include \
	      -I/usr/local/include -Dpublic= -Dprivate=static -DHAVE_CONFIG_H -c 
rms_dist_objects: $(RMS_DIST_OBJECTS_O)
	$(CC) $(OPT) $(RMS_DIST_OBJECTS_O) -o $@ -lbicpl -L../Proglib \
		-L/nil/david/Release/BIC_PL -L/nil/david/Release/Volume_io \
		-L/usr/local/lib -L/nil/david/Release/BIC_PL/IRIX_5 \
		-lmincprog -lvolume_io -lminc -lnetcdf -lc_s  -lm

MST_LINES_OBJS = Extra_progs/mst_to_lines.o 
mst_to_lines: $(MST_LINES_OBJS)
	$(CC) $(OPT) $(MST_LINES_OBJS) -o $@ -lbicpl $(LDOPT) 

TEST_LINES_OBJS = Extra_progs/test_lines.o 
test_lines: $(TEST_LINES_OBJS)
	$(CC) $(OPT) $(TEST_LINES_OBJS) -o $@ -lbicpl $(LDOPT) 

TEST_DERIV_OBJS = Extra_progs/test_deriv.o Numerical/quad_max_fit.o
test_deriv: $(TEST_DERIV_OBJS)
	$(CC) $(OPT) $(TEST_DERIV_OBJS) -o $@ $(LDOPT)

FLIPOBJS = Extra_progs/reversedef.o Optimize/my_grid_support.o
reversedef: $(FLIPOBJS)
	$(CC) $(OPT) $(FLIPOBJS) -o $@ $(LDOPT)

INBETWEEN = Extra_progs/inbetween_def.o 
inbetween_def: $(INBETWEEN)
	$(CC) $(OPT) $(INBETWEEN) -o $@ $(LDOPT)


RANDOMOBJS = Extra_progs/random_def.o Numerical/random_numbers.o
random_def: $(RANDOMOBJS)
	$(CC) $(OPT) $(RANDOMOBJS) -o $@ $(LDOPT)


RANDOMVOLOBJS = Extra_progs/random_vol.o 
random_vol: $(RANDOMVOLOBJS)
	$(CC) $(OPT) $(RANDOMVOLOBJS) -o $@ $(LDOPT)

RANDOMIZEOBJS = Extra_progs/randomize_volume.o 
randomize_volume: $(RANDOMIZEOBJS)
	$(CC) $(OPT) $(RANDOMIZEOBJS) -o $@ $(LDOPT)


AVGVOXVALOBJS = Extra_progs/avg_voxel_vol.o 
avg_voxel_val: $(AVGVOXVALOBJS)
	$(CC) $(OPT) $(AVGVOXVALOBJS) -o $@ $(LDOPT)


PARAM_OBJS = Extra_progs/param_test.o Numerical/make_rots.o Numerical/matrix_basics.o Numerical/rotmat_to_ang.o 
param_test: $(PARAM_OBJS)
	$(CC) $(OPT) $(PARAM_OBJS) -o $@ $(LDOPT) 

XFM2PARAM_OBJS = \
	Extra_progs/xfm2param.o \
	Numerical/make_rots.o \
	Numerical/matrix_basics.o \
	Numerical/rotmat_to_ang.o \
	Volume/interpolation.o
xfm2param: $(XFM2PARAM_OBJS)
	$(CC) $(OPT) $(XFM2PARAM_OBJS) -o $@ $(LDOPT) 

PARAM2_OBJS = \
	Extra_progs/param2xfm.o \
	Numerical/make_rots.o \
	Numerical/matrix_basics.o \
	Numerical/rotmat_to_ang.o

param2xfm: $(PARAM2_OBJS)
	$(CC) $(OPT) $(PARAM2_OBJS) -o $@ $(LDOPT) 



CHECK_OBJS = \
	Extra_progs/check_scale.o \
	Numerical/make_rots.o \
	Numerical/matrix_basics.o \
	Numerical/rotmat_to_ang.o \
	Volume/interpolation.o

check_scale: $(CHECK_OBJS)
	$(CC) $(OPT) $(CHECK_OBJS) -o $@ $(LDOPT) 

xcorr_vol: Extra_progs/xcorr_vol.o
	$(CC) $(OPT) Extra_progs/xcorr_vol.o -o $@ $(LDOPT) 

overlap_corr: Extra_progs/overlap_corr.o
	$(CC) $(OPT) Extra_progs/overlap_corr.o -o $@ $(LDOPT) 

avg_voxel_vol: Extra_progs/avg_voxel_vol.o
	$(CC) $(OPT) Extra_progs/avg_voxel_vol.o -o $@ $(LDOPT) 


zscore_vol: Extra_progs/zscore_vol.o  Volume/interpolation.o
	$(CC) $(OPT) Extra_progs/zscore_vol.o Volume/interpolation.o -o $@ $(LDOPT)

MAKE_ZSCORE_OBJ = \
		Extra_progs/make_zscore_vol.o \
		Volume/interpolation.o \
		Volume/volume_functions.o
make_zscore_vol: $(MAKE_ZSCORE_OBJ)
	$(CC) $(OPT) $(MAKE_ZSCORE_OBJ) -o $@ $(LDOPT)

change_range: Extra_progs/change_range.o
	$(CC) $(OPT) Extra_progs/change_range.o -o $@ $(LDOPT)

load_vol: Extra_progs/load_vol.o
	$(CC) $(OPT) Extra_progs/load_vol.o -o $@ $(LDOPT)


RAND_PARAM_OBJS = \
	Extra_progs/rand_param.o \
	Numerical/make_rots.o \
	Numerical/matrix_basics.o \
	Numerical/rotmat_to_ang.o 

rand_param: $(RAND_PARAM_OBJS)
	$(CC) $(OPT) $(RAND_PARAM_OBJS) -o $@ $(LDOPT) 

randomize_tags: Extra_progs/randomize_tags.o Numerical/random_numbers.o
	$(CC) $(OPT) Extra_progs/randomize_tags.o Numerical/random_numbers.o -o $@ $(LDOPT) 

COV_OBJS = \
	Extra_progs/test_cov.o \
	Numerical/make_rots.o \
	Numerical/matrix_basics.o \
	Numerical/rotmat_to_ang.o \
	Numerical/cov_to_praxes.o \
	Volume/interpolation.o

test_cov: $(COV_OBJS)
	$(CC) $(OPT) $(COV_OBJS) -o $@ $(LDOPT) 

FIND_OBJS = \
	Extra_progs/find_tubes.o \
	Numerical/make_rots.o \
	Numerical/matrix_basics.o \
	Numerical/rotmat_to_ang.o \
	Numerical/cov_to_praxes.o

find_tubes: $(FIND_OBJS)
	$(CC) $(OPT) $(FIND_OBJS) -o $@ $(LDOPT) 

test_cov_pr: Extra_progs/test_cov_pr.o Numerical/cov_to_praxes.o  Numerical/quad_max_fit.o 
	$(CC) $(OPT) Extra_progs/test_cov_pr.o Numerical/quad_max_fit.o Numerical/cov_to_praxes.o -o $@ $(LDOPT)

volume_cog: Extra_progs/volume_cog.o Volume/interpolation.o
	$(CC) $(OPT) Extra_progs/volume_cog.o Volume/interpolation.o -o $@ $(LDOPT) 

volume_stat: Extra_progs/volume_stat.o Volume/interpolation.o
	$(CC) $(OPT) Extra_progs/volume_stat.o Volume/interpolation.o -o $@ $(LDOPT) 


test3d2.xfm: animal
	./animal Testing/object1_blur.mnc Testing/object2_blur.mnc $@ \
		-feature Testing/object1_dxyz.mnc Testing/object2_dxyz.mnc xcorr 0.7 \
		-nonlin -optical -step 4 4 4 -iter 1 -weight 0.3 -identity \
		-debug -clobber




test3d.xfm: animal
	./animal Testing/object1_blur.mnc Testing/object2_blur.mnc $@ -optical \
        -feature Testing/object1_chamfer.mnc Testing/object2_chamfer.mnc chamfer \
        -nonlin  -step 4 4 4 -iter 10 -weight 0.7 -stiff 0.7 -identity -sub_lattice 5\
		-debug -clobber ; \
        mincresample Testing/sulci1.mnc -like Testing/object2.mnc Testing/object1-like-2.mnc \
                -transformation test3d.xfm -clobber -nearest

include ../Makefile.bottom




# --------------------------------------------------------------------


#depends:

Main/minctracc.o: \
	Include/minctracc.h \
	Include/arg_data.h \
	Include/globals.h \
	Include/point_vector.h \
	Include/constants.h \
	Include/objectives.h \
	Include/local_macros.h \
	Main/measure_code.c \
	../Proglib/print_error.h 

Optimize/do_nonlinear.o: \
	Include/stats.h \
	Include/constants.h \
	Include/local_macros.h \
	Include/arg_data.h 


Numerical/stats.o: Include/stats.h

Optimize/sub_lattice.o: \
	Optimize/switch_obj_func.c
