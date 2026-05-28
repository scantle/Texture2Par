# compiler
FC = gfortran

# compile flags
FCFLAGS =  -O2 -cpp -fbacktrace -ffree-line-length-none  -ffixed-line-length-none -fdec -fallow-argument-mismatch -fopenmp
#
#TMPDIR := .mod
#ifneq ($(TMPDIR),)
#  $(shell mkdir -p $(TMPDIR))
#  FCFLAGS+= -J $(TMPDIR)
#endif

# link flags
FLFLAGS = -fopenmp

# set src dir
VPATH =                                \
    ../src/                                \
    ../src/Libraries/flibs/src/strings/    \
    ../src/Libraries/mfUSGLib/Jan_23_2018/ \
    ../src/Libraries/mfUSGLib/             \
    ../                                   \

# source files and objects
flibs = \
    m_vstring.f90               \
    m_vstringlist.F90

#    pmmodule_usg.f90            \
#    ppcgu1.f                    \

mfusg = \
    mf_fpath.f90                \
    sparse.f                    \
    gmodules.f                  \
    xmdlib_2.f                  \
    xmd.f                       \
    lak_gag_sfr_modules.f       \
    glo2p3du1.f                 \
    gwf2QRT8u.f                 \
    gwf2str7u1.f                \
    gwf2drt8u.f                 \
    gwf2dpf1u1.f                \
    pcgu7.f                     \
    gwf2wel7u1.f                \
    tvmu2.f                     \
    gwf2hfb7u1.f                \
    glo2sgbu1.f                 \
    disu2gncn1.f                \
    tvmu1.f                     \
    gwf2ghb7u1.f                \
    gwf2chd7u1.f                \
    gwf2fhb7u1.f                \
    gwf2sub7u1.f                \
    gwf2mnw17.f                 \
    gwf2riv7u1.f                \
    gwf2mnw27.f                 \
    gwf2drn7u1.f                \
    glo2btnu1.f                 \
    gwf2rch8u1.f                \
    gwf2ets8u1.f                \
    gwt2dptu1.f                 \
    gwf2evt8u1.f                \
    density.f                   \
    gwt2bndsu1.f                \
    gwf2basu1.f                 \
    glo2sms-u1.f                \
    gwf2lak7u1.f                \
    disu2gncb1.f                \
    gwf2gag7u1.f                \
    gwf2bcf-lpf-u1.f            \
    glo2basu1.f                 \
    gwf2sfr7u1.f                \
    parutl7.f                   \
    cln2props1.f                \
    utl7u1.f                    \
    cln2basu1.f                 \
    mf_t2p.f90

t2p2 = \
    Str_Module.f90              \
    sgemv.f                     \
    time.f90                    \
    kdtree2.f90                 \
    Lists.f90                   \
    fpath.f90                   \
    FileIO.f90                  \
    global.f90                  \
    ErrorHandler.f90            \
    variogram.f90               \
    data_transformer.f90        \
    Sparse_Symmetric.f90        \
    Categories.f90              \
    Tools.f90                   \
    cellblocks.f90              \
    Grid.f90                    \
    DataSets.f90                \
    weights.f90                 \
    PilotPoints.f90             \
    Options.f90                 \
    write_fmodel.f90            \
    Interpolators.f90           \
    kriging.f90                 \
    sgsim.f90                   \
    FlowModel.f90               \
    ColMap.f90                  \
    subsidence_funcs.f90        \
    T2P.f90                     \
    version.f90                 \
    t2p_log.f90                 \
    Read_Datasets.f90           \
    Read_Main_Input.f90         \
    Texture2Par.f90

objflibs = $(addsuffix .o, $(flibs))
objother = $(addsuffix .o, $(mfusg) $(t2p2))
objects = $(objflibs) $(objother)
# display all targets
$(info $$objflibs = $(objflibs))
$(info $$objother = $(objother))

# program name
PROGRAM = ../bin/Texture2Par

all: $(PROGRAM)

# gnu.org/software/make/manual/html_node/Automatic-Variables.html
# % match any part of the target name
# $@ -- The file name of the target of the rule.In a pattern rule that has multiple targets, $@ is the name of whichever target caused the rule's recipe to be run.
# $% -- The target member name, when the target is an archive member.
# $< -- The name of the first prerequisite.
# $? -- The names of all the prerequisites that are newer than the target, with spaces between them.
# $^ -- The names of all the prerequisites, with spaces between them.
# $+ -- This is like ‘$^’, but prerequisites listed more than once are duplicated in the order they were listed in the makefile.
# $| -- The names of all the order-only prerequisites, with spaces between them.
# $* -- The stem with which an implicit rule matches

$(PROGRAM): $(objects)
	$(FC) -static $(FLFLAGS) -o $@ $^

$(objflibs): %.o: %
	$(FC) $(FCFLAGS) -c                 -o $@ $<

$(objother): %.o: %
	$(FC) $(FCFLAGS) -c -freal-4-real-8 -o $@ $<


clean:
	rm -rf *.mod
	rm -rf *.o
	rm -rf $(PROGRAM)*
	#rm -rf $(TMPDIR)\*