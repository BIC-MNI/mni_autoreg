# --------------------------------------------------------------------
# 
# Makefile for mincbbox
#

DAVE_LIB = volume_io
include ../Makefile.include

       PROG  = mincbbox
    PROG_OBJ = $(PROG).o
  EXTRA_SRCS = 
 HEADER_SRCS =  $(PROG).h
  EXTRA_OBJS = $(EXTRA_SRCS:.c=.o)

LINT_LIST1 = $(PROG_OBJ:.o=.ln) $(EXTRA_OBJS:.o=.ln)        # list of lint files in program 1

ALL_PROGS = $(PROG) 
ALL_LINTS = lint_$(PROG)
RCS_LIST = $(PROG).c $(EXTRA_SRCS) $(HEADER_SRCS)


$(PROG): $(PROG_OBJ) $(EXTRA_OBJS)
	$(CC) $(OPT) $(PROG_OBJ) $(EXTRA_OBJS) -o $@ $(LDOPT)

$(PROG).o: $(PROG).h

lint_$(PROG): $(LINT_LIST1)                # how to lint the executable source
	lint -u $(LINTFLAGS) $(LINT_LIST1)

include ../Makefile.bottom



