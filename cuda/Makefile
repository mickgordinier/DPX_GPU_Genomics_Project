NVCC = nvcc
NVCC_FLAGS =

CPP_SOURCES = ../c++/backtrack.cpp \
							../c++/parseInput.cpp \
							../c++/timing.cpp \
							../c++/printLock.cpp
				
LNW_CU_SOURCES := $(wildcard LNW/LinearNeedlemanWunsch*.cu)
LNW_EXECUTABLES := $(notdir $(basename $(LNW_CU_SOURCES)))

MANUAL_EXECUTABLES = BSW ANW LSW SW SWOPT BWS test

ALL_EXECUTABLES = $(LNW_EXECUTABLES) $(MANUAL_EXECUTABLES)

print_vars:
	@echo "LNW Source Files: $(LNW_CU_SOURCES)"
	@echo "LNW Executables:  $(LNW_EXECUTABLES)"
	@echo "Manual Targets:   $(MANUAL_EXECUTABLES)"

# Pattern rule for automatically detected LNW-style builds
$(LNW_EXECUTABLES): %: LNW/%.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -o $@


LSW: LinearSmithWaterman.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -o $@

SWOPT: LinearSmithWatermanOpt.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -o $@

BSW: LinearSmithWaterman.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -o $@

ANW: AffineNeedlemanWunsch.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -o $@

clean:
	rm -f $(ALL_EXECUTABLES) *.out
