
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

RLXX=$(LFLAGS) $(ROOTLIBS) -lMinuit -lTreePlayer -g
# add ROOT
CXXFLAGS += $(ROOTCFLAGS)
# add degub and warnings
CXXFLAGS += -g -Wall
# add include folders
CXXFLAGS += -I./base/ -I./Reconstruction/ -I../utils/

BASESRC 		= 	base/
BASESEL 		=		Reconstruction/
SRC 				= 	$(notdir $(wildcard $(BASESRC)/*.$(SrcSuf))) $(notdir $(wildcard $(BASESEL)/*.$(SrcSuf)))
OBJ 				= 	$(addprefix $(OBJDIR)/, $(SRC:.$(SrcSuf)=.$(ObjSuf)))

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

all: build SpatialResol dEdx cleanObj

# declare your analysis like this
SpatialResol: $(OBJ) $(OBJDIR)/SpatialResolAna.o
	$(CXX) $^ $(OutPutOpt) $(BINDIR)/$@.$(OutPutSuf) $(RLXX)
	@echo "Compilation done"
# end of your analysis declaration

dEdx: $(OBJ) $(OBJDIR)/dEdxAna.o
	$(CXX) $^ $(OutPutOpt) $(BINDIR)/$@.$(OutPutSuf) $(RLXX)
	@echo "Compilation done"

$(OBJDIR)/%.o: */%.cxx
	$(CXX) -c $(OutPutOpt) $@ $< $(CXXFLAGS)

.PHONY: build
build:
	@mkdir -p $(OBJDIR)
	@mkdir -p $(BINDIR)

.PHONY: clean
clean:
	@echo "Cleaning..."
	rm -f $(OBJDIR)/*.$(ObjSuf)
	rm -r $(BINDIR)/*.$(OutPutSuf)

.PHONY: cleanObj
cleanObj:
	rm -f $(OBJDIR)/*.$(ObjSuf)
