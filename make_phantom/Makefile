# --------------------------------------------------------------------
# 
# Makefile for minc_phantom
#

include ../Makefile.include

       PROG = make_phantom
   PROG_OBJ = $(PROG).o
 EXTRA_SRCS = 
HEADER_SRCS = make_phantom.h 
 EXTRA_OBJS = $(EXTRA_SRCS:.c=.o) 


LINT_LIST = $(PROG_OBJ:.o=.ln) $(EXTRA_OBJS:.o=.ln)  

ALL_PROGS = $(PROG) 
ALL_LINTS = lint_$(PROG) 

$(PROG): $(PROG_OBJ) $(EXTRA_OBJS) $(HEADER_SRCS) 
	$(CC) $(OPT) $(PROG_OBJ) $(EXTRA_OBJS) -o $@ $(LDOPT)

lint_$(PROG): $(LINT_LIST)                # how to lint the executable source
	lint -u $(LINTFLAGS) $(LINT_LIST)

include ../Makefile.bottom

