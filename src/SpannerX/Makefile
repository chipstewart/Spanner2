# =============================
# Spanner SpannerScan Makefile
# (c) 2011 Chip Stewart
# =============================

SOURCES= SpannerX.cpp
OBJECTS= $(SOURCES:.cpp=.o)
_EXT_OBJECTS=BamScan.o FastaFile.o Pars.o BamUtil.o Histo.o UtilityFunctions.o BamX.o FragmentPosFile.o

EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
INCLUDES = -I$(SRC_DIR) -I$(SAMTOOLS_DIR)

PROGRAM= SpannerX


all: $(PROGRAM)

.PHONY: all

$(PROGRAM): $(BUILT_OBJECTS) $(EXT_OBJECTS)
	@echo "  * linking $(PROGRAM)"
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/$@ $^ $(LIBS)

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

$(EXT_OBJECTS):
	$(MAKE) --no-print-directory -C $(SRC_DIR)
		
clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean
