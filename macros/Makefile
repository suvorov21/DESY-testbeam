# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

RC     := root-config
ifeq ($(shell which $(RC) 2>&1 | sed -ne "s@.*/$(RC)@$(RC)@p"),$(RC))
MKARCH := $(wildcard $(shell $(RC) --etcdir)/Makefile.arch)
endif
ifneq ($(MKARCH),)
include $(MKARCH)
else
$(error ROOTSYS is not defined)
include $(ROOTSYS)/Makefile.arch
endif
-include ../MyConfig.mk
#------------------------------------------------------------------------------
##CXX           =
ObjSuf        = o
SrcSuf        = cxx
OutPutOpt     = -o
OutPutSuf     = exe
BINDIR        = ../bin/
OBJDIR 				=	../obj/

# add utils: styles, geometry, etc
SRC=$(SRC)

OBJS= $(OBJDIR)/$(notdir $(SRC:.$(SrcSuf)=.$(ObjSuf)))
MCR=  $(BINDIR)/$(notdir $(SRC:.$(SrcSuf)=.$(OutPutSuf)))
RLXX=$(LFLAGS) $(ROOTLIBS) -lMinuit -lTreePlayer

##turning on debuging
RLXX +=  -g

$(info Beginning of compilation)
$(info )
$(info Compiler : $(CXX) $(CXXFLAGS))
$(info Linker   : $(CXX) $(LDFLAGS))
$(info RLXX     : $(RLXX))
$(info )
$(info Sources          : $(SRC))
$(info Objects          : $(OBJS))
$(info Binary directory : $(BINDIR))
$(info Executables      : $(MCR))
$(info -------------------------------------)
$(info )

all: build $(MCR)

$(MCR): $(OBJDIR)/$(OBJS)
	$(LD) $^ $(OutPutOpt) $@.$(OutPutSuf) $(RLXX)
	rm $(OBJDIR)/*.$(ObjSuf)
	@echo "Compilation done"

$(OBJDIR)/$(OBJS): $(SRC)
	$(LD) -c $^ $(CXXFLAGS) -o $(OBJDIR)/$(OBJS)

.PHONY: build
build:
	@mkdir -p $(OBJDIR)
	@mkdir -p $(BINDIR)

.PHONY: clean
clean:
	@echo "Cleaning..."
	rm -f $(OBJDIR)/*.$(ObjSuf)
	rm -r $(BINDIR)/$(MCR)

