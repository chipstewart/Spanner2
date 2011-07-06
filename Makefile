# ==========================
# BEDTools Makefile
# (c) 2009 Aaron Quinlan
# ==========================

# define our object and binary directories
export OBJ_DIR	= $(abspath ./obj)
export BIN_DIR	= $(abspath ./bin)
export SRC_DIR	= $(abspath ./src)
export CXX		= g++
export CXXFLAGS = -Wall -O2
# define samtools area for sam.h and libbam.a
export SAMTOOLS_DIR = /usr/local/samtools-0.1.16
export LIBS		= -L$(SAMTOOLS_DIR) -lbam -lz 


SUBDIRS = $(SRC_DIR) \
		  $(SRC_DIR)/SpannerScan \
		  $(SRC_DIR)/SpannerX \
		  
all:
	[ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR)
	[ -d $(BIN_DIR) ] || mkdir -p $(BIN_DIR)
	
	@echo "Building Spanner:"
	@echo "========================================================="
	
	@for dir in $(SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $$dir; \
		echo ""; \
	done


.PHONY: all

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean
