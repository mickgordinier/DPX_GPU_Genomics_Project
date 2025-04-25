NVCC = nvcc
NVCC_FLAGS =

CPP_SOURCES = 	../c++/backtrack.cpp \
				../c++/parseInput.cpp \
				../c++/timing.cpp \
				../c++/printLock.cpp
				
CU_SOURCES = $(filter-out $(CU_FILTER), $(wildcard *.cu))
CU_HEADERS = $(wildcard *.h)
CU_EXE = $(CU_SOURCES:.cu=) LNW ANW LNW2 LNW3 LNW8 LNW9 LNW10 LNW11 LNW12 LNW13 LNW14 LNW15 LNW17 LNW18 LSW


print_vars:
	@echo "Source Files $(CU_SOURCES)"
	@echo "Header Files $(CU_HEADERS)"
	@echo "Object Files $(CU_EXE)"

test: test.cu
	$(NVCC) -o test test.cu

LSW: LinearSmithWaterman.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -o $@

LNW: LinearNeedlemanWunsch.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -o $@

LNW2: LinearNeedlemanWunschV2.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -o $@

LNW3: LinearNeedlemanWunschV3.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -o $@

LNW4: LinearNeedlemanWunschV4.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -o $@

LNW8: LinearNeedlemanWunschV8.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -g -o $@

LNW9: LinearNeedlemanWunschV9.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -g -o $@

LNW10: LinearNeedlemanWunschV10.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -g -o $@

LNW12: LinearNeedlemanWunschV12.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -g -o $@

LNW13: LinearNeedlemanWunschV13.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -g -o $@

LNW14: LinearNeedlemanWunschV14.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -g -o $@

LNW15: LinearNeedlemanWunschV15.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -g -o $@

LNW17: LinearNeedlemanWunschV17.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -g -o $@

LNW18: LinearNeedlemanWunschV18.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -g -o $@

ANW: AffineNeedlemanWunsch.cu $(CPP_SOURCES)
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -o $@

%: %.cu
	$(NVCC) $(NVCC_FLAGS) $< $(CPP_SOURCES) -o $@

%.out: %
	./$< > $@
	
clean:
	rm -f $(CU_EXE)
