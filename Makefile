# 
# Top level makefile for MNI AutoReg package
#

SHELL = /bin/sh

ALL_SUBDIRS = \
	Proglib \
	mincblur \
	minccrop \
	make_phantom \
	minctracc \
	perl

TARGET = $@

default :  
	@echo ""
	@echo "In order to make the MNI registration utilities:"
	@echo ""
	@echo "   0) you must have an ANSI-compliant C compiler"
	@echo "   1) you may have to build the NETCDF library"
	@echo "   2) you may have to build the Minc library"
	@echo "   3) you may have to build the stereotaxic model"
	@echo "   4) edit Makefile.include to set paths and machine type"
	@echo "   5) type:"
	@echo "        make build"
	@echo "        make test"
	@echo "        make install"
	@echo "        make installman"
	@echo ""
 


all clean build install: 
	@for dir in $(ALL_SUBDIRS); \
	   do if [ -d $$dir ]; \
	   then \
	      cd $$dir; \
	      echo "" ; echo "" ; \
	      echo "*************************************"; \
	      echo Making $(TARGET) in $$dir; \
	      echo "*************************************"; echo "" ; \
	      $(MAKE) $(TARGET); \
	      cd ..; \
	   fi; \
	 done

distclean: clean testclean
	\rm -f config.cache config.log config.status Makefile.include \
               config.h version.h

testclean:
	\rm -f Testing/*.mnc Testing/*.act Testing/*.xfm

test:
	cd Testing; \
	$(MAKE) $(TARGET); \
	cd ..

installman:
	cd minctracc; \
	$(MAKE) $(TARGET); \
	cd ..
	cd mincblur; \
	$(MAKE) $(TARGET); \
	cd ..
	cd perl; \
	$(MAKE) $(TARGET); \
	cd ..
