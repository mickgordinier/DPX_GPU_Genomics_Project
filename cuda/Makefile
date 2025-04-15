NVCC = nvcc
NVCC_FLAGS = 

CPP_SOURCES = 	../c++/backtrack.cpp \
				../c++/parseInput.cpp \
				../c++/timing.cpp
				
CU_SOURCES = $(filter-out $(CU_FILTER), $(wildcard *.cu))
CU_HEADERS = $(wildcard *.h)
CU_EXE = $(CU_SOURCES:.cu=) LNW ANW LNW2 LNW3

print_vars:
	@echo "Source Files $(CU_SOURCES)"
	@echo "Header Files $(CU_HEADERS)"
	@echo "Object Files $(CU_EXE)"

test: test.cu
	$(NVCC) -o test test.cu

LNW: LinearNeedlemanWunsch.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -o $@

LNW2: LinearNeedlemanWunschV2.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -o $@

LNW3: LinearNeedlemanWunschV3.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -o $@

ANW: AffineNeedlemanWunsch.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -o $@

%: %.cu
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -o $@

%.out: %
	./$< > $@
	
clean:
	rm -f $(CU_EXE)

