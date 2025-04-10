NVCC = nvcc
NVCC_FLAGS = 


CU_TARGET  = main.cu
CU_FILTER = $(CU_TARGET) testFakeDPX.cu
CU_SOURCES = $(filter-out $(CU_FILTER), $(wildcard *.cu))
CU_HEADERS = $(wildcard *.h)
CU_OBJECTS = $(CU_SOURCES:.cu=.o)

print_vars:
	@echo "Target $(CU_TARGET)"
	@echo "Source Files $(CU_SOURCES)"
	@echo "Header Files $(CU_HEADERS)"
	@echo "Object Files $(CU_OBJECTS)"

test: test.cu
	$(NVCC) -o test test.cu

%: %.cu
	$(NVCC) $(NVCC_FLAGS) $(CU_HEADERS) $% -o $@

%.out: %
	./$< > $@
	
clean:
	rm -f main $(CPP_OBJECTS)

