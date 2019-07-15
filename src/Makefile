
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
# add include folders
CXXFLAGS += $(ROOTCFLAGS) -I./base/ -I./Selection/ -I../utils/

BASESRC			= base/
BASESEL 		= Selection/
SOURCEFILES = $(wildcard $(BASESRC)/*.$(SrcSuf)) $(wildcard $(BASESEL)/*.$(SrcSuf))
OBJECTFILES = $(patsubst $(BASESRC)%.$(SrcSuf),$(OBJDIR)%.$(ObjSuf),$(SOURCEFILES))

$(info Beginning of compilation)
$(info )
$(info Compiler : $(CXX) $(CXXFLAGS))
$(info Linker   : $(CXX) $(LDFLAGS))
$(info RLXX     : $(RLXX))
$(info )
$(info Base Sources     : $(SOURCEFILES))
$(info Base Objects     : $(OBJECTFILES))
$(info Binary directory : $(BINDIR))
$(info -------------------------------------)
$(info )

all: build SpatialResol

# declare your analysis like this
SpatialResol: $(OBJECTFILES) $(OBJDIR)/SpatialResolAna.o
	$(CXX) $^ $(OutPutOpt) $(BINDIR)/$@.$(OutPutSuf) $(RLXX)
	rm $(OBJDIR)/*.$(ObjSuf)
	@echo "Compilation done"

$(OBJDIR)/%.o: SpatialResol/%.cxx SpatialResol/%.hxx
	$(CXX) -c $(OutPutOpt) $@ $< $(CXXFLAGS)
# end of your analysis declaration

$(OBJDIR)/%.o: base/%.cxx
	$(CXX) -c $(OutPutOpt) $@ $< $(CXXFLAGS)

$(OBJDIR)/%.o: Selection/%.cxx
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

