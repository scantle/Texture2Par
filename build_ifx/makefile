# Texture2Par command-line build for Intel IFX using GNU make.
#
# Build from this directory after the Intel oneAPI environment is loaded:
#
#   make -f makefile rebuild
#   make -f makefile CFG=debug rebuild
#
# On Windows, if GNU make is installed as mingw32-make:
#
#   mingw32-make -f makefile rebuild
#   mingw32-make -f makefile CFG=debug rebuild

FC = ifx
CFG ?= release

ifeq ($(OS),Windows_NT)
  RC ?= rc
  AR = lib
  EXEEXT := .exe
  LIBEXT := .lib
  MKDIR_P = if not exist "$(subst /,\,$1)" mkdir "$(subst /,\,$1)"
  RM_RF = if exist "$(subst /,\,$1)" rmdir /s /q "$(subst /,\,$1)"
  RM_F = if exist "$(subst /,\,$1)" del /q "$(subst /,\,$1)"
  T2P_COMMON = /nologo /fpp /libs:static /threads /winapp /Qm64
  FLIB_COMMON = /nologo /fpp /libs:static /threads /Qm64
  MFUSG_COMMON = /nologo /libs:static /threads /Qm64
  IFX_MOD = /module:$(MODDIR) /I$(MODDIR)
  IFX_REAL8 = /real-size:64
  IFX_FIXED = /extfor:f
  IFX_OBJ = /object:$@
  IFX_EXE = /exe:$@
  IFX_LINK = /nologo /libs:static /Qm64
  INC = /I
  ARFLAGS = /nologo /machine:x64 /out:$@
  RCFLAGS = /nologo /I$(SRCDIR)
else
  $(error This GNU makefile is currently verified only for Windows with Intel IFX)
endif

ifeq ($(CFG),debug)
  BUILDDIR := build/debug
  T2P_OPTFLAGS := /Od /debug:full /check:bounds /warn:interfaces /traceback /Zi /dbglibs
  LIB_OPTFLAGS := /Od /debug:full /check:bounds /warn:interfaces /traceback /Zi /dbglibs
  LINKDEBUG := /debug:full
else
  BUILDDIR := build/release
  T2P_OPTFLAGS := /O2 /traceback
  LIB_OPTFLAGS := /O2
  LINKDEBUG :=
endif

OBJDIR := $(BUILDDIR)/obj
MODDIR := $(BUILDDIR)/mod
FLIBMODDIR := $(MODDIR)/flibs
T2PMODDIR := $(MODDIR)/Texture2Par
LIBDIR := $(BUILDDIR)/lib
BINDIR := ../bin
PROGRAM := $(BINDIR)/Texture2Par$(EXEEXT)
MFUSGMODDIR := $(MODDIR)/mfUSGLib
FLIB_LIB := $(LIBDIR)/flibs$(LIBEXT)
MFUSG_LIB := $(LIBDIR)/mfUSGLib$(LIBEXT)
FLIB_DEP = $(objflibs)

LIBS := $(FLIB_LIB) $(MFUSG_LIB)

FLIBDIR := ../src/Libraries/flibs/src/strings
MFUSGDIR := ../src/Libraries/mfUSGLib
MFUSG2018DIR := ../src/Libraries/mfUSGLib/Jan_23_2018
SRCDIR := ../src
ROOTDIR := ..

VPATH := $(FLIBDIR) $(MFUSGDIR) $(MFUSG2018DIR) $(SRCDIR) $(ROOTDIR)

FLIB_MODFLAGS := $(subst $(MODDIR),$(FLIBMODDIR),$(IFX_MOD))
MFUSG_MODFLAGS := $(subst $(MODDIR),$(MFUSGMODDIR),$(IFX_MOD))
T2P_MODFLAGS := $(subst $(MODDIR),$(T2PMODDIR),$(IFX_MOD)) $(INC)$(FLIBMODDIR) $(INC)$(MFUSGMODDIR)

FLIB_FFLAGS := $(FLIB_COMMON) $(FLIB_MODFLAGS) $(INC)$(FLIBDIR) $(LIB_OPTFLAGS)
MFUSG_FFLAGS := $(MFUSG_COMMON) $(MFUSG_MODFLAGS) $(LIB_OPTFLAGS)
T2P_FFLAGS := $(T2P_COMMON) $(T2P_MODFLAGS) $(INC)$(FLIBDIR) $(T2P_OPTFLAGS) $(IFX_REAL8)
LDFLAGS := $(IFX_LINK) $(LINKDEBUG)

flibs := \
    m_vstring.f90 \
    m_vstringlist.F90

mfusg := \
    sparse.f \
    mf_fpath.f90 \
    pmmodule_usg.f90 \
    xmdlib_2.f \
    gmodules.f \
    gwf2fhb7u1.f \
    gwf2ghb7u1.f \
    pcgu7.f \
    gwf2mnw17.f \
    gwf2hfb7u1.f \
    gwf2sub7u1.f \
    gwf2str7u1.f \
    ppcgu1.f \
    gwf2wel7u1.f \
    gwf2chd7u1.f \
    gwf2QRT8u.f \
    gwf2mnw27.f \
    gwf2dpf1u1.f \
    xmd.f \
    gwf2drt8u.f \
    glo2p3du1.f \
    lak_gag_sfr_modules.f \
    gwf2riv7u1.f \
    disu2gncn1.f \
    gwf2drn7u1.f \
    glo2sgbu1.f \
    tvmu1.f \
    glo2btnu1.f \
    gwf2ets8u1.f \
    gwf2rch8u1.f \
    tvmu2.f \
    gwt2dptu1.f \
    gwf2evt8u1.f \
    density.f \
    gwt2bndsu1.f \
    utl7u1.f \
    cln2basu1.f \
    gwf2basu1.f \
    glo2sms-u1.f \
    gwf2sfr7u1.f \
    cln2props1.f \
    parutl7.f \
    gwf2gag7u1.f \
    disu2gncb1.f \
    gwf2bcf-lpf-u1.f \
    gwf2lak7u1.f \
    mf_t2p.f90 \
    glo2basu1.f

t2p2 := \
    sgemv.f \
    Tools.f90 \
    kdtree2.f90 \
    ColMap.f90 \
    version.f90 \
    fpath.f90 \
    Lists.f90 \
    Str_Module.f90 \
    time.f90 \
    cellblocks.f90 \
    FileIO.f90 \
    global.f90 \
    ErrorHandler.f90 \
    Options.f90 \
    variogram.f90 \
    weights.f90 \
    data_transformer.f90 \
    Grid.f90 \
    Sparse_Symmetric.f90 \
    Categories.f90 \
    DataSets.f90 \
    Interpolators.f90 \
    write_fmodel.f90 \
    kriging.f90 \
    PilotPoints.f90 \
    FlowModel.f90 \
    subsidence_funcs.f90 \
    sgsim.f90 \
    T2P.f90 \
    Read_Main_Input.f90 \
    t2p_log.f90 \
    Read_Datasets.f90 \
    Texture2Par.f90

objflibs := $(addprefix $(OBJDIR)/,$(addsuffix .obj,$(basename $(flibs))))
objmfusg := $(addprefix $(OBJDIR)/,$(addsuffix .obj,$(basename $(mfusg))))
objt2p := $(addprefix $(OBJDIR)/,$(addsuffix .obj,$(basename $(t2p2))))
objects := $(objt2p) $(LIBS)
$(objflibs): FFLAGS := $(FLIB_FFLAGS)
$(objmfusg): FFLAGS := $(MFUSG_FFLAGS)
$(objt2p): FFLAGS := $(T2P_FFLAGS)
resources := $(OBJDIR)/T2P.res

.PHONY: all release debug dirs clean rebuild

all: dirs $(PROGRAM)

release:
	$(MAKE) -f makefile CFG=release all

debug:
	$(MAKE) -f makefile CFG=debug all

rebuild:
	$(MAKE) -f makefile CFG=$(CFG) clean
	$(MAKE) -f makefile CFG=$(CFG) all

dirs:
	$(call MKDIR_P,$(BUILDDIR))
	$(call MKDIR_P,$(OBJDIR))
	$(call MKDIR_P,$(MODDIR))
	$(call MKDIR_P,$(FLIBMODDIR))
	$(call MKDIR_P,$(MFUSGMODDIR))
	$(call MKDIR_P,$(T2PMODDIR))
	$(call MKDIR_P,$(LIBDIR))
	$(call MKDIR_P,$(BINDIR))

$(PROGRAM): $(objects) $(resources) | dirs
	$(FC) $(LDFLAGS) $(IFX_EXE) $(objects) $(resources)

$(FLIB_LIB): $(objflibs) | dirs
	$(AR) $(ARFLAGS) $(objflibs)

$(MFUSG_LIB): $(objmfusg) | dirs
	$(AR) $(ARFLAGS) $(objmfusg)

$(OBJDIR)/T2P.res: $(SRCDIR)/T2P.rc $(SRCDIR)/resource.h $(SRCDIR)/T2P.ico | dirs
	$(RC) $(RCFLAGS) /fo$@ $<

$(OBJDIR)/m_vstring.obj: m_vstring.f90 | dirs
	$(FC) $(FLIB_FFLAGS) $(IFX_OBJ) -c $<

$(OBJDIR)/m_vstringlist.obj: m_vstringlist.F90 $(OBJDIR)/m_vstring.obj | dirs
	$(FC) $(FLIB_FFLAGS) $(IFX_OBJ) -c $<

$(OBJDIR)/%.obj: %.f90 | dirs
	$(FC) $(FFLAGS) $(IFX_OBJ) -c $<

$(OBJDIR)/%.obj: %.F90 | dirs
	$(FC) $(FFLAGS) $(IFX_OBJ) -c $<

$(OBJDIR)/%.obj: %.f | dirs
	$(FC) $(FFLAGS) $(IFX_FIXED) $(IFX_OBJ) -c $<

clean:
	$(call RM_RF,$(BUILDDIR))
	$(call RM_F,$(PROGRAM))
	$(call RM_F,$(BINDIR)/Texture2Par.pdb)
	$(call RM_F,$(BINDIR)/Texture2Par.ilk)

# Module-level prerequisites. The object list above is intentionally ordered
# for clean builds; these dependencies cover normal incremental GNU make builds.
$(OBJDIR)/xmd.obj: $(OBJDIR)/xmdlib_2.obj $(OBJDIR)/gmodules.obj
$(OBJDIR)/ppcgu1.obj: $(OBJDIR)/pmmodule_usg.obj
$(OBJDIR)/glo2p3du1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/gwf2QRT8u.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/gwf2str7u1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/gwf2drt8u.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/gwf2dpf1u1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/gwf2wel7u1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/tvmu2.obj: $(OBJDIR)/gmodules.obj $(OBJDIR)/lak_gag_sfr_modules.obj $(OBJDIR)/gwf2dpf1u1.obj
$(OBJDIR)/gwf2hfb7u1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/glo2sgbu1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/disu2gncn1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/tvmu1.obj: $(OBJDIR)/gmodules.obj $(OBJDIR)/lak_gag_sfr_modules.obj $(OBJDIR)/gwf2dpf1u1.obj
$(OBJDIR)/gwf2ghb7u1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/gwf2chd7u1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/gwf2fhb7u1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/gwf2sub7u1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/gwf2riv7u1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/gwf2drn7u1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/glo2btnu1.obj: $(OBJDIR)/gmodules.obj $(OBJDIR)/gwf2dpf1u1.obj
$(OBJDIR)/gwf2rch8u1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/gwf2ets8u1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/gwt2dptu1.obj: $(OBJDIR)/gmodules.obj $(OBJDIR)/gwf2dpf1u1.obj
$(OBJDIR)/gwf2evt8u1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/gwt2bndsu1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/gwf2basu1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/glo2sms-u1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/gwf2lak7u1.obj: $(OBJDIR)/gmodules.obj $(OBJDIR)/lak_gag_sfr_modules.obj
$(OBJDIR)/gwf2gag7u1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/gwf2bcf-lpf-u1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/glo2basu1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/gwf2sfr7u1.obj: $(OBJDIR)/gmodules.obj $(OBJDIR)/lak_gag_sfr_modules.obj
$(OBJDIR)/parutl7.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/utl7u1.obj: $(OBJDIR)/gmodules.obj
$(OBJDIR)/mf_t2p.obj: $(OBJDIR)/mf_fpath.obj $(OBJDIR)/gmodules.obj

$(OBJDIR)/FileIO.obj: $(OBJDIR)/Str_Module.obj $(FLIB_DEP) $(OBJDIR)/Lists.obj
$(OBJDIR)/global.obj: $(OBJDIR)/FileIO.obj
$(OBJDIR)/ErrorHandler.obj: $(OBJDIR)/global.obj $(FLIB_DEP)
$(OBJDIR)/variogram.obj: $(OBJDIR)/global.obj $(OBJDIR)/FileIO.obj $(FLIB_DEP) $(OBJDIR)/ErrorHandler.obj
$(OBJDIR)/data_transformer.obj: $(FLIB_DEP) $(OBJDIR)/ErrorHandler.obj
$(OBJDIR)/Sparse_Symmetric.obj: $(OBJDIR)/ErrorHandler.obj
$(OBJDIR)/Categories.obj: $(FLIB_DEP) $(OBJDIR)/FileIO.obj $(OBJDIR)/ErrorHandler.obj
$(OBJDIR)/Grid.obj: $(OBJDIR)/kdtree2.obj $(OBJDIR)/cellblocks.obj $(OBJDIR)/ErrorHandler.obj $(OBJDIR)/Tools.obj
$(OBJDIR)/DataSets.obj: $(OBJDIR)/global.obj $(OBJDIR)/Grid.obj $(OBJDIR)/Categories.obj
$(OBJDIR)/weights.obj: $(OBJDIR)/FileIO.obj $(OBJDIR)/ErrorHandler.obj
$(OBJDIR)/PilotPoints.obj: $(OBJDIR)/DataSets.obj $(OBJDIR)/global.obj $(OBJDIR)/Grid.obj $(FLIB_DEP) $(OBJDIR)/FileIO.obj $(OBJDIR)/ErrorHandler.obj
$(OBJDIR)/Options.obj: $(OBJDIR)/global.obj $(OBJDIR)/FileIO.obj $(FLIB_DEP) $(OBJDIR)/ErrorHandler.obj $(OBJDIR)/Tools.obj
$(OBJDIR)/write_fmodel.obj: $(OBJDIR)/ErrorHandler.obj $(OBJDIR)/FileIO.obj $(OBJDIR)/DataSets.obj $(MFUSG_LIB)
$(OBJDIR)/Interpolators.obj: $(OBJDIR)/global.obj $(OBJDIR)/weights.obj $(OBJDIR)/Grid.obj $(OBJDIR)/DataSets.obj $(OBJDIR)/ErrorHandler.obj $(OBJDIR)/Options.obj $(OBJDIR)/variogram.obj $(OBJDIR)/Categories.obj
$(OBJDIR)/kriging.obj: $(OBJDIR)/global.obj $(OBJDIR)/Interpolators.obj $(OBJDIR)/weights.obj $(OBJDIR)/Grid.obj $(OBJDIR)/variogram.obj $(OBJDIR)/DataSets.obj $(OBJDIR)/ErrorHandler.obj $(OBJDIR)/Options.obj $(FLIB_DEP)
$(OBJDIR)/sgsim.obj: $(OBJDIR)/global.obj $(OBJDIR)/Interpolators.obj $(OBJDIR)/weights.obj $(OBJDIR)/Grid.obj $(OBJDIR)/variogram.obj $(OBJDIR)/DataSets.obj $(OBJDIR)/ErrorHandler.obj $(OBJDIR)/Options.obj $(OBJDIR)/Sparse_Symmetric.obj $(OBJDIR)/kriging.obj $(OBJDIR)/kdtree2.obj
$(OBJDIR)/FlowModel.obj: $(OBJDIR)/FileIO.obj $(OBJDIR)/Grid.obj $(OBJDIR)/ErrorHandler.obj $(OBJDIR)/Categories.obj $(OBJDIR)/Options.obj $(OBJDIR)/cellblocks.obj $(OBJDIR)/global.obj $(OBJDIR)/DataSets.obj $(OBJDIR)/fpath.obj $(OBJDIR)/Interpolators.obj $(OBJDIR)/write_fmodel.obj $(MFUSG_LIB)
$(OBJDIR)/subsidence_funcs.obj: $(OBJDIR)/global.obj $(OBJDIR)/DataSets.obj $(OBJDIR)/Interpolators.obj $(OBJDIR)/kriging.obj $(OBJDIR)/variogram.obj $(OBJDIR)/Options.obj $(OBJDIR)/ErrorHandler.obj $(OBJDIR)/Grid.obj
$(OBJDIR)/T2P.obj: $(OBJDIR)/global.obj $(FLIB_DEP) $(OBJDIR)/ErrorHandler.obj $(OBJDIR)/Options.obj $(OBJDIR)/DataSets.obj $(OBJDIR)/FlowModel.obj $(OBJDIR)/PilotPoints.obj $(OBJDIR)/Interpolators.obj $(OBJDIR)/kriging.obj $(OBJDIR)/sgsim.obj $(OBJDIR)/Categories.obj $(OBJDIR)/variogram.obj $(OBJDIR)/Grid.obj $(OBJDIR)/subsidence_funcs.obj $(OBJDIR)/ColMap.obj
$(OBJDIR)/t2p_log.obj: $(OBJDIR)/FileIO.obj $(OBJDIR)/version.obj $(OBJDIR)/time.obj $(OBJDIR)/T2P.obj $(OBJDIR)/ColMap.obj $(OBJDIR)/ErrorHandler.obj $(FLIB_DEP)
$(OBJDIR)/Read_Datasets.obj: $(OBJDIR)/global.obj $(OBJDIR)/T2P.obj $(OBJDIR)/Str_Module.obj $(FLIB_DEP) $(OBJDIR)/FileIO.obj $(OBJDIR)/ErrorHandler.obj $(OBJDIR)/Grid.obj $(OBJDIR)/data_transformer.obj $(OBJDIR)/subsidence_funcs.obj $(OBJDIR)/t2p_log.obj $(OBJDIR)/ColMap.obj
$(OBJDIR)/Read_Main_Input.obj: $(OBJDIR)/T2P.obj $(FLIB_DEP) $(OBJDIR)/FileIO.obj $(OBJDIR)/ErrorHandler.obj $(OBJDIR)/FlowModel.obj $(OBJDIR)/ColMap.obj $(OBJDIR)/variogram.obj $(OBJDIR)/kriging.obj $(OBJDIR)/Tools.obj $(OBJDIR)/PilotPoints.obj $(OBJDIR)/global.obj $(OBJDIR)/subsidence_funcs.obj
$(OBJDIR)/Texture2Par.obj: $(OBJDIR)/global.obj $(OBJDIR)/version.obj $(OBJDIR)/T2P.obj $(OBJDIR)/FileIO.obj $(OBJDIR)/Read_Main_Input.obj $(OBJDIR)/t2p_log.obj $(OBJDIR)/ErrorHandler.obj $(OBJDIR)/Read_Datasets.obj
