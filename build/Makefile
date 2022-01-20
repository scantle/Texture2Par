# makefile created on 2021-07-30 14:34:49.615882
# by pymake (version 1.2.1) for the 'mf6' executable 
# using the 'gfortran' fortran compiler(s).

PROGRAM = Texture2Par

FC = gfortran

# set the optimization level (OPTLEVEL) if not defined
OPTLEVEL ?= -O2


MODSWITCH = -J mod

# determine OS
ifeq ($(OS), Windows_NT)
	detected_OS = Windows
	OS_macro = -D_WIN32
else
	detected_OS = $(shell sh -c 'uname 2>/dev/null || echo Unknown')
	ifeq ($(detected_OS), Darwin)
		OS_macro = -D__APPLE__
	else
		OS_macro = -D__LINUX__
	endif
endif

# set the fortran flags
ifeq ($(detected_OS), Windows)
	ifeq ($(FC), gfortran)
		FFLAGS ?= -static -fbacktrace $(OS_macro) -Wtabs -Wline-truncation -Wunused-label -Wunused-variable -pedantic -Wcharacter-truncation -cpp -freal-4-real-8
	endif
else
	ifeq ($(FC), gfortran)
		FFLAGS ?= -static -fbacktrace $(OS_macro) -Wtabs -Wline-truncation -Wunused-label -Wunused-variable -pedantic -Wcharacter-truncation -cpp -freal-4-real-8
	endif
	ifeq ($(FC), $(filter $(FC), ifort mpiifort))
		FFLAGS ?= -no-heap-arrays -fpe0 -traceback -fpp /real-size:64
		MODSWITCH = -module mod
	endif
endif

# set the ldflgs
ifeq ($(detected_OS), Windows)
	ifeq ($(FC), $(filter $(FC), gfortran))
		LDFLAGS ?= -static -static-libgfortran -static-libgcc -static-libstdc++ -lm
	endif
else
	ifeq ($(FC), gfortran)
		LDFLAGS ?= -lc
	endif
	ifeq ($(FC), $(filter $(FC), ifort mpiifort))
		LDFLAGS ?= -lc
	endif
endif

VPATH = \
../mf6/Model/GroundWaterFlow  \
../mf6/Model/ModelUtilities   \
../mf6/Timing                 \
../mf6/Utilities              \
../mf6/Utilities/Memory       \
../mf6/Utilities/TimeSeries   \
../mfUSGLib                   \
../mfUSGLib/Jan_23_2018       \
..

OBJMF6 = ./obj/mf6
OBJUSG = ./obj/usg
OBJT2P = ./obj/t2p

INCSWITCH = -I $(OBJMF6) -I $(OBJUSG) -I $(OBJT2P)


mf6 = \
$(OBJMF6)/kind.o                    \
$(OBJMF6)/Sparse.o                  \
$(OBJMF6)/Constants.o               \
$(OBJMF6)/compilerversion.o         \
$(OBJMF6)/SmoothingFunctions.o      \
$(OBJMF6)/OpenSpec.o                \
$(OBJMF6)/HashTable.o               \
$(OBJMF6)/defmacro.o                \
$(OBJMF6)/SimVariables.o            \
$(OBJMF6)/genericutils.o            \
$(OBJMF6)/List.o                    \
$(OBJMF6)/TimeSeriesRecord.o        \
$(OBJMF6)/version.o                 \
$(OBJMF6)/ArrayHandlers.o           \
$(OBJMF6)/Message.o                 \
$(OBJMF6)/Sim.o                     \
$(OBJMF6)/InputOutput.o             \
$(OBJMF6)/ArrayReaders.o            \
$(OBJMF6)/MemoryHelper.o            \
$(OBJMF6)/BlockParser.o             \
$(OBJMF6)/TimeSeries.o              \
$(OBJMF6)/DisvGeom.o                \
$(OBJMF6)/TimeSeriesLink.o          \
$(OBJMF6)/TableTerm.o               \
$(OBJMF6)/TimeSeriesFileList.o      \
$(OBJMF6)/Table.o                   \
$(OBJMF6)/Memory.o                  \
$(OBJMF6)/ListReader.o              \
$(OBJMF6)/MemoryList.o              \
$(OBJMF6)/MemoryManager.o           \
$(OBJMF6)/Connections.o             \
$(OBJMF6)/MemorySetHandler.o        \
$(OBJMF6)/ats.o                     \
$(OBJMF6)/tdis.o                    \
$(OBJMF6)/TimeSeriesManager.o       \
$(OBJMF6)/DiscretizationBase.o      \
$(OBJMF6)/gwf3disu8.o               \
$(OBJMF6)/gwf3dis8.o                \
$(OBJMF6)/gwf3disv8.o               \
$(OBJMF6)/TimeArray.o               \
$(OBJMF6)/TimeArraySeries.o         \
$(OBJMF6)/TimeArraySeriesLink.o     \
$(OBJMF6)/TimeArraySeriesManager.o

mfusg = \
$(OBJUSG)/gmodules.o \
$(OBJUSG)/lak_gag_sfr_modules.o \
$(OBJUSG)/pcgu7.o \
$(OBJUSG)/gwf2dpf1u1.o \
$(OBJUSG)/glo2p3du1.o \
$(OBJUSG)/xmdlib_2.o \
$(OBJUSG)/xmd.o \
$(OBJUSG)/tvmu1.o \
$(OBJUSG)/tvmu2.o \
$(OBJUSG)/glo2btnu1.o \
$(OBJUSG)/gwf2evt8u1.o \
$(OBJUSG)/gwf2riv7u1.o \
$(OBJUSG)/gwf2hfb7u1.o \
$(OBJUSG)/gwf2chd7u1.o \
$(OBJUSG)/gwf2ghb7u1.o \
$(OBJUSG)/gwf2rch8u1.o \
$(OBJUSG)/gwf2drn7u1.o \
$(OBJUSG)/gwf2drt8u.o \
$(OBJUSG)/gwf2QRT8u.o \
$(OBJUSG)/gwf2mnw17.o \
$(OBJUSG)/gwf2mnw27.o \
$(OBJUSG)/gwf2sub7u1.o \
$(OBJUSG)/gwf2wel7u1.o \
$(OBJUSG)/gwf2fhb7u1.o \
$(OBJUSG)/gwf2str7u1.o \
$(OBJUSG)/gwf2bcf-lpf-u1.o \
$(OBJUSG)/glo2sgbu1.o \
$(OBJUSG)/disu2gncn1.o \
$(OBJUSG)/cln2basu1.o \
$(OBJUSG)/cln2props1.o \
$(OBJUSG)/parutl7.o \
$(OBJUSG)/utl7u1.o \
$(OBJUSG)/gwf2basu1.o \
$(OBJUSG)/gwf2sfr7u1.o \
$(OBJUSG)/gwf2lak7u1.o \
$(OBJUSG)/glo2basu1.o \
$(OBJUSG)/gwf2gag7u1.o \
$(OBJUSG)/gwf2ets8u1.o \
$(OBJUSG)/gwt2bndsu1.o \
$(OBJUSG)/gwt2dptu1.o \
$(OBJUSG)/density.o \
$(OBJUSG)/glo2sms-u1.o

texture2par = \
$(OBJT2P)/fpath.o \
$(OBJT2P)/KDTree.o \
$(OBJT2P)/Tools.o \
$(OBJT2P)/utl.o \
$(OBJT2P)/spkrige.o \
$(OBJT2P)/MakePar.o \
$(OBJT2P)/mf_t2p.o \
$(OBJT2P)/mfsupport.o \
$(OBJT2P)/writeLPF.o \
$(OBJT2P)/readIWFM.o \
$(OBJT2P)/io.o \
$(OBJT2P)/Texture2Par.o


objects = $(mf6) $(mfusg) $(texture2par)



# Define the objects that make up the program
$(PROGRAM) : makedirs $(objects)
	$(FC) $(OPTLEVEL) -o $@ $(objects) $(LDFLAGS)

$(OBJMF6)/%.o : %.f90
	$(FC) $(OPTLEVEL) $(FFLAGS) -c $< -o $@ -I $(OBJMF6) $(MODSWITCH)

$(OBJUSG)/%.o : %.f
	$(FC) $(OPTLEVEL) $(FFLAGS) -c $< -o $@ -I $(OBJUSG) $(MODSWITCH)
	
$(OBJT2P)/%.o : %.f90
	$(FC) $(OPTLEVEL) $(FFLAGS) -c $< -o $@ -I $(OBJT2P) -I $(OBJMF6) -I $(OBJUSG) $(MODSWITCH)


# Make the bin directory for the executable
.PHONY: makedirs
makedirs:
	mkdir -p $(OBJUSG)
	mkdir -p $(OBJMF6)
	mkdir -p $(OBJT2P)
	mkdir -p mod


.PHONY: clean
clean:
	rm -rf $(OBJUSG)
	rm -rf $(OBJMF6)
	rm -rf $(OBJT2P)
	rm -rf mod