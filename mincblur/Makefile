include ../Makefile.include

       PROG  = mincblur
    PROG_OBJ = $(PROG).o
  EXTRA_SRCS =  blur_support.c \
		gradient_volume.c \
		blur_volume.c \
		gradmag_volume.c \
		apodize_data.c \
		fft.c
 HEADER_SRCS =  blur_support.h   mincblur.h gradmag_volume.h
  EXTRA_OBJS = $(EXTRA_SRCS:.c=.o)

    MANPAGES = mincblur.1


LINT_LIST1 = $(PROG_OBJ:.o=.ln) $(EXTRA_OBJS:.o=.ln)        # list of lint files in program 1

ALL_PROGS = $(PROG) 
ALL_LINTS = lint_$(PROG) 

$(PROG): $(PROG_OBJ) $(EXTRA_OBJS) $(PROG).h $(HEADER_SRCS) 
	$(CC) $(OPT) $(PROG_OBJ) $(EXTRA_OBJS) -o $@ $(LDOPT)

lint_$(PROG): $(LINT_LIST1)            
	lint -u $(LINTFLAGS) $(LINT_LIST1)


include ../Makefile.bottom

