#
# Makefile for minctracc and associated programs for linear
# registration
# 
# version 0.1e
# $Id: Makefile,v 1.3 2000-01-19 18:47:48 louis Exp $
# $Log: Makefile,v $
# Revision 1.3  2000-01-19 18:47:48  louis
# Updating the Makefile for nonlinear minctracc
#
# Revision 1.2  2000/01/19 18:18:03  louis
# Updated minctracc/Makefile to be the saem as in mni_autoreg_0.98f
#
# Revision 1.2  1995/09/11  12:50:01  collins
# removed Numerical/transformations.c from source list, since it is no
# longer needed.
#
# Revision 1.1  1995/08/29  11:21:12  collins
# Initial revision
#
#

include ../Makefile.include

EXTRA_PROGS = \
	Extra_progs/check_scale.c \
	Extra_progs/param2xfm.c \
	Extra_progs/rand_param.c \
	Extra_progs/xcorr_vol.c \
	Extra_progs/xfm2param.c \
	Extra_progs/test_cov.c \
	Extra_progs/volume_cog.c \
	Extra_progs/reverse_def.c \
	Extra_progs/zscore_vol.c

  PROG = minctracc

MANPAGES=minctracc.1

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

 HEADER_SRCS =  \
	Include/boolean.h \
	Include/cov_to_praxes.h \
	Include/def_segment_table.h \
	Include/def_tps.h \
	Include/globals.h \
	Include/make_rots.h \
	Include/matrix_basics.h \
	Include/minctracc.h \
	Include/def_tag_io.h \
	Include/rotmat_to_ang.h 

  OBJS = \
	$(SRCS:.c=.o)

  LINT_LIST1 = \
	$(PROG_OBJ:.o=.ln) $(EXTRA_OBJS:.o=.ln)


 ALL_PROGS = \
	$(PROG) \
	check_scale param2xfm volume_cog \
	rand_param xcorr_vol xfm2param zscore_vol

   ALL_LINTS = lint_$(PROG)


default: $(ALL_PROGS)

$(PROG_LIB_DIR)/lib$(PROG_LIB).a:
	cd $(PROG_LIB_DIR); make build

$(ALL_PROGS): $(PROG_LIB_DIR)/lib$(PROG_LIB).a

$(PROG): $(OBJS) $(HEADER_SRCS)
	$(CC) $(OPT)  $(OBJS) -o $@ $(LDOPT)

lint_$(PROG): $(LINT_LIST1)                # how to lint the executable source
	lint -u $(LINTFLAGS) $(LINT_LIST1)

CHECK_OBJS = \
	Extra_progs/check_scale.o \
	Numerical/make_rots.o \
	Numerical/matrix_basics.o \
	Numerical/rotmat_to_ang.o \
	Volume/interpolation.o

check_scale: $(CHECK_OBJS)
	$(CC) $(OPT) $(CHECK_OBJS) -o $@ $(LDOPT)

COV_OBJS = \
	Extra_progs/test_cov.o \
	Numerical/make_rots.o \
	Numerical/matrix_basics.o \
	Numerical/rotmat_to_ang.o \
	Numerical/cov_to_praxes.o \
	Volume/interpolation.o

test_cov: $(COV_OBJS)
	$(CC) $(OPT) $(COV_OBJS) -o $@ $(LDOPT)

xcorr_vol: Extra_progs/xcorr_vol.o
	$(CC) $(OPT) Extra_progs/xcorr_vol.o -o $@ $(LDOPT)

zscore_vol: Extra_progs/zscore_vol.o  Volume/interpolation.o
	$(CC) $(OPT) Extra_progs/zscore_vol.o Volume/interpolation.o -o $@ $(LDOPT)

PARAM2_OBJS = \
	Extra_progs/param2xfm.o \
	Numerical/make_rots.o \
	Numerical/matrix_basics.o \
	Numerical/rotmat_to_ang.o

param2xfm: $(PARAM2_OBJS)
	$(CC) $(OPT) $(PARAM2_OBJS) -o $@ $(LDOPT)

XFM2_OBJS = \
	Extra_progs/xfm2param.o \
	Numerical/make_rots.o \
	Numerical/matrix_basics.o \
	Numerical/rotmat_to_ang.o \
	Volume/interpolation.o

xfm2param: $(XFM2_OBJS)
	$(CC) $(OPT) $(XFM2_OBJS) -o $@ $(LDOPT)

RAND_PARAM_OBJS = \
	Extra_progs/rand_param.o \
	Numerical/make_rots.o \
	Numerical/matrix_basics.o \
	Numerical/rotmat_to_ang.o

rand_param: $(RAND_PARAM_OBJS)
	$(CC) $(OPT) $(RAND_PARAM_OBJS) -o $@ $(LDOPT)

volume_cog: Extra_progs/volume_cog.o Volume/interpolation.o
	$(CC) $(OPT) Extra_progs/volume_cog.o Volume/interpolation.o -o $@ $(LDOPT)


reversedef: Extra_progs/reversedef.o 
	$(CC) $(OPT) Extra_progs/reversedef.o -o $@ $(LDOPT)

include ../Makefile.bottom







