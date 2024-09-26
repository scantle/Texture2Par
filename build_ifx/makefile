# compiler
FC = ifx

# compile flags
FCFLAGS = -traceback -heap-arrays -O2 -fpe:0 -fpp -nologo -static -real-size:64
#
TMPDIR := .temp
ifneq ($(TMPDIR),)
  $(shell mkdir $(TMPDIR))
  MODSWITCH = -module:$(TMPDIR)
  FCFLAGS+= $(MODSWITCH)
endif

# link flags
FLFLAGS =

# set src dir
VPATH =                                \
    ../                                \
    ../Libraries/flibs/src/strings/    \
    ../Libraries/mfUSGLib/Jan_23_2018/ \
    ../Libraries/mfUSGLib/             \

# source files and objects
flibs = \
    m_vstring.F90               \
    m_vstringlist.F90

#    ppcgu1.f                    \
#    pmmodule_usg.f90            \

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
    T2P.f90                     \
    t2p_log.f90                 \
    Read_Datasets.f90           \
    Read_Main_Input.f90         \
    Texture2Par.f90

srcs    = $(flibs) $(mfusg) $(t2p2)
objects = $(addsuffix .obj, $(srcs))
object2 = $(addprefix $(TMPDIR)/, $(objects))
# display all targets
$(info $$objects = $(objects))

# program name
PROGRAM = texture2par

all: $(PROGRAM)

# $@ -- target name
# $< -- prerequisites
$(PROGRAM): $(objects)
	$(FC) $(MODSWITCH) -o $@ $(objects)

# % is the
%.obj: %
	$(FC) $(FCFLAGS) -c -o $@ $<

%.F90:

%.f90:

%.f:

%.F:

clean:
	rm -rf $(TMPDIR)\*
	rm -rf *.obj
	rm -rf $(PROGRAM)*