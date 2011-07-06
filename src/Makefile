# =============================
# Spanner SpannerScan Makefile
# (c) 2011 Chip Stewart
# =============================
#export OBJ_DIR = $(abspath ../obj)
#export BIN_DIR = $(abspath ../bin)
#export INC_DIR = $(abspath ../)

INCLUDES = -I$(SRC_DIR) -I$(SAMTOOLS_DIR)

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= BamScan.cpp Histo.cpp BamUtil.cpp FastaFile.cpp Pars.cpp UtilityFunctions.cpp BamX.cpp
OBJECTS= $(SOURCES:.cpp=.o)
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean
