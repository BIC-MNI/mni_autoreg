SUBDIRS = Files Numerical Optimize Volume Main Extra_progs

# Automake does not handle sources in subdirs (yet).  We kludge our way
# around this by making a bunch of (noninstalled) libraries --- one per subdir ---
# and building the minctracc program in subdirectory Main.  The extra programs
# are likewise built in their subdirectory, Extra_progs.

EXTRA_DIST = $(man_MANS) $(docs) $(includes) $(matlab) $(notes)

man_MANS = minctracc.1

docs = \
	Doc/doc \
	Doc/minctracc.doc

includes = \
	Include/amoeba.h \
	Include/arg_data.h \
	Include/constants.h \
	Include/cov_to_praxes.h \
	Include/def_segment_table.h \
	Include/def_tag_io.h \
	Include/def_tps.h \
	Include/deform_support.h \
	Include/extras.h \
	Include/globals.h \
	Include/local_macros.h \
	Include/make_rots.h \
	Include/matrix_basics.h \
	Include/minctracc.h \
	Include/objectives.h \
	Include/point_vector.h \
	Include/quad_max_fit.h \
	Include/rotmat_to_ang.h \
	Include/segment_table.h \
	Include/stats.h \
	Include/sub_lattice.h \
	Include/vox_space.h

matlab = \
	Matlab/avg.prof \
	Matlab/blurdata.m \
	Matlab/boxdata.m \
	Matlab/brute.m \
	Matlab/comments \
	Matlab/comments.tex \
	Matlab/corrdat.m \
	Matlab/corrdata.m \
	Matlab/cubic.m \
	Matlab/cubic_p.m \
	Matlab/cubic_pp.m \
	Matlab/derivdat.m \
	Matlab/differ.m \
	Matlab/do_realwarp.m \
	Matlab/do_test.m \
	Matlab/do_warp.m \
	Matlab/do_wrp.old \
	Matlab/estderi2.m \
	Matlab/estderiv.m \
	Matlab/gold.m \
	Matlab/int_like.m \
	Matlab/inter.m \
	Matlab/inter_1p.m \
	Matlab/inter_2p.m \
	Matlab/inter_p.m \
	Matlab/make_fig.m \
	Matlab/makedat2.m \
	Matlab/makedata.m \
	Matlab/makedata2.m \
	Matlab/makeft.m \
	Matlab/makekern.m \
	Matlab/makeoffs.m \
	Matlab/match.m \
	Matlab/normdist.m \
	Matlab/op_node.m \
	Matlab/op_node.working \
	Matlab/op_node.working2 \
	Matlab/plot_res.m \
	Matlab/profile \
	Matlab/profile2 \
	Matlab/recover.m \
	Matlab/rectdist.m \
	Matlab/set_rang.m \
	Matlab/subj.prof \
	Matlab/test.m \
	Matlab/test_res.m \
	Matlab/testcorr.m \
	Matlab/todo \
	Matlab/transf.m \
	Matlab/warpdata.m

notes = \
	Notes/limited_fov

