# --------------------------------------------------------------------
#
# MINC Makefile
#

include ../Makefile.include

# Object names
OBJS     = ParseArgv.o vax_conversions.o time_stamp.o print_error.o print_version.o
RM       = rm
RM_FLAGS = -f

# --------------------------------------------------------------------

CC_PROG_LIB= lib$(PROG_LIB).a

LINT_OBJS = $(OBJS:.o=.ln)
LINT_PROG_LIB = llib-l$(PROG_LIB).ln

# --------------------------------------------------------------------

ALL_PROGS = $(CC_PROG_LIB) $(NR_PROG_LIB)

#Dependency on Makefile
$(OBJS) $(LINT_OBJS) : Makefile

all: $(CC_PROG_LIB) $(NR_PROG_LIB) lint

# How to create library
$(CC_PROG_LIB): $(OBJS)
	$(RM) $(RM_FLAGS) $(CC_PROG_LIB)
	$(AR) $(AR_FLAGS) $(CC_PROG_LIB) $(OBJS)
	$(RANLIB) $(CC_PROG_LIB)


# Lint
lint : $(LINT_PROG_LIB)

# How to create lint library
$(LINT_PROG_LIB) : $(LINT_OBJS)
	lint -o $(PROG_LIB) -x -u $(LINTFLAGS) $(LINT_OBJS)


# Remove all derived files in this directory

-:

test: test.o
	$(CC) $@.o -o $@ $(LDOPT)

clean:
	\rm -f   $(OBJS) $(LINT_OBJS) $(CC_PROG_LIB) $(NR_PROG_LIB)

build: $(CC_PROG_LIB) $(NR_PROG_LIB)

install:
	@ echo Nothing to install here.
