# directory structure
DIR_MAIN    = ./
DIR_BUILD   = $(DIR_MAIN)build/
DIR_H       = $(DIR_BUILD)include/
DIR_OBJ     = $(DIR_BUILD)obj/
DIR_CXX     = $(DIR_BUILD)src/
DIR_MAC     = $(DIR_MAIN)macros/

# search paths
vpath %.h   $(DIR_H)
vpath %.cxx $(DIR_CXX)
vpath %.o   $(DIR_OBJ)

# file lists
BIN_FEMTO  = therm2_femto
HSRC_FEMTO = Parser.cxx Configurator.cxx ParticleDB.cxx ParticleType.cxx DecayTable.cxx DecayChannel.cxx Messages.cxx
SRC_FEMTO  = $(HSRC_FEMTO:%=$(DIR_CXX)%) $(BIN_FEMTO:%=$(DIR_CXX)%.cxx)
OBJ_FEMTO  = $(SRC_FEMTO:$(DIR_CXX)%.cxx=$(DIR_OBJ)%.o)

# preprocessor
PREPROCESS  = -D_CXX_VER_="\"$(shell $(CXX) --version | grep $(CXX))\"" -D_ROOT_VER_="\"$(shell root-config --version)\""
ifdef DEBUG
  PREPROCESS  := $(PREPROCESS) -D_DEBUG_LEVEL_=$(DEBUG)
else
  PREPROCESS  := $(PREPROCESS) -D_DEBUG_LEVEL_=0
endif

# compilation
CXX         = g++
LD          = g++
CXXFLAGS    = -std=c++11 -O0 -g -Wall -Wextra -Wpedantic -I $(DIR_H) $(PREPROCESS) `root-config --cflags`
LFLAGS      = -lm -g `root-config --libs`

#################################################################################
# RULES                                                                         #
#################################################################################
 
all: $(BIN_FEMTO:%=$(DIR_OBJ)%)
	cp $^ $(DIR_MAIN)
	echo
	echo "Ready!"
	echo "Type \"./therm2_femto\" to generate two-particle corelation function"
	echo

$(DIR_OBJ)therm2_femto: $(OBJ_FEMTO)
	echo "Linking:   $@ ($(LD))"
	$(LD) $^ -o $@ $(LFLAGS)

$(DIR_OBJ)%.o: %.cxx
	@[ -d $(DIR_OBJ) ] || mkdir -p $(DIR_OBJ)
	echo "Compiling: $< ($(CXX))"
	$(CXX) -c $< -o $@ $(CXXFLAGS)

clean:
	rm -f $(DIR_OBJ)*.o
	rm -f $(DIR_MAC)*.d
	rm -f $(DIR_MAC)*.so
	rm -f $(DIR_MAC)*.pcm
	rm -f $(DIR_OBJ)$(BIN_FEMTO) $(DIR_MAIN)$(BIN_FEMTO)
	echo "*.o, *.so, *.d, *.pcm and binary files removed."

.SILENT :