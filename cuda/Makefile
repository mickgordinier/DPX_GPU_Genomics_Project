NVCC = nvcc
NVCC_FLAGS = 

CPP_SOURCES = 	../c++/backtrack.cpp \
				../c++/parseInput.cpp
				
CU_SOURCES = $(filter-out $(CU_FILTER), $(wildcard *.cu))
CU_HEADERS = $(wildcard *.h)
CU_EXE = $(CU_SOURCES:.cu=) LNW

print_vars:
	@echo "Source Files $(CU_SOURCES)"
	@echo "Header Files $(CU_HEADERS)"
	@echo "Object Files $(CU_EXE)"

test: test.cu
	$(NVCC) -o test test.cu

LNW: LinearNeedlemanWunschV4.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -o $@

%: %.cu
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -o $@

%.out: %
	./$< > $@
	
clean:
	rm -f $(CU_EXE)

