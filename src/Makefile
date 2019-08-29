ROOTCFLAGS 	= $(shell root-config --cflags)
ROOTLIBS   	= $(shell root-config --libs)
ROOTGLIBS  	= $(shell root-config --glibs)
ROOTLDFLAGS	= $(shell root-config --ldflags)

##CXX           =
ObjSuf        = o
SrcSuf        = cxx
OutPutOpt     = -o
OutPutSuf     = exe

BINDIR        = ../bin/
OBJDIR 				=	../obj/

DICT 					= eventdict

RLXX=$(LFLAGS) $(ROOTLIBS) -lMinuit -lTreePlayer
# add ROOT and warnings
CXXFLAGS += $(ROOTCFLAGS) -Wall
# for release
ifeq ($(TYPE), RELEASE)
	CXXFLAGS += -O3
	RLXX += -O3
else
	# add degub
	CXXFLAGS += -g
	RLXX += -g
endif
# add include folders
CXXFLAGS += -I./base/DataStorage/ -I./base/ -I./Reconstruction/ -I../utils/

BASESRC 		= 	base/
BASEDATASRC = 	base/DataStorage/
BASESEL 		=		Reconstruction/
#SRC 				= 	$(notdir $(wildcard $(BASESRC)/*.$(SrcSuf))) $(notdir $(wildcard $(BASESEL)/*.$(SrcSuf)))
SRC 				= 	$(notdir $(wildcard $(BASEDATASRC)/*.$(SrcSuf))) $(notdir $(wildcard $(BASESRC)/*.$(SrcSuf))) $(notdir $(wildcard $(BASESEL)/*.$(SrcSuf)))
OBJ 				= 	$(OBJDIR)$(DICT).$(ObjSuf) $(addprefix $(OBJDIR)/, $(SRC:.$(SrcSuf)=.$(ObjSuf)))

$(info Beginning of compilation)
$(info )
$(info Compiler : $(CXX) $(CXXFLAGS))
$(info Linker   : $(CXX) $(LDFLAGS))
$(info RLXX     : $(RLXX))
$(info )
$(info Base Sources     : $(SRC))
$(info Base Objects     : $(OBJ))
$(info Binary directory : $(BINDIR))
$(info -------------------------------------)
$(info )

all: build SpatialResol dEdx

# declare your analysis like this
SpatialResol: $(OBJ) $(OBJDIR)/SpatialResolAna.o
	$(CXX) $^ $(OutPutOpt) $(BINDIR)/$@.$(OutPutSuf) $(RLXX)
	@echo "Compilation done"
# end of your analysis declaration

dEdx: $(OBJ) $(OBJDIR)/dEdxAna.o
	$(CXX) $^ $(OutPutOpt) $(BINDIR)/$@.$(OutPutSuf) $(RLXX)
	@echo "Compilation done"

$(OBJDIR)$(DICT).$(ObjSuf):
	rootcint -v1 -f $(BINDIR)$(DICT).$(SrcSuf) -c $(addprefix ../src/, $(wildcard $(BASEDATASRC)/*.hxx)) ../src/$(BASEDATASRC)/linkDef.h
	$(CXX) -c $(OutPutOpt) $@ $(BINDIR)$(DICT).$(SrcSuf) $(CXXFLAGS)

$(OBJDIR)/%.o: */DataStorage/%.cxx
	$(CXX) -c $(OutPutOpt) $@ $< $(CXXFLAGS)

$(OBJDIR)/%.o: */%.cxx
	$(CXX) -c $(OutPutOpt) $@ $< $(CXXFLAGS)

.PHONY: build
build:
	@mkdir -p $(OBJDIR)
	@mkdir -p $(BINDIR)

.PHONY: clean
clean:
	@echo "Cleaning..."
	rm -rf $(OBJDIR)/*
	rm -rf $(BINDIR)/*

.PHONY: cleanObj
cleanObj:
	rm -f $(OBJDIR)/*.$(ObjSuf)
