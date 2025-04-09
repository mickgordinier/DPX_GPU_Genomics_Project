GCC = g++
GCC_FLAGS = 


CPP_TARGET  = main.cpp
CPP_FILTER = $(CPP_TARGET) testFakeDPX.cpp
CPP_SOURCES = $(filter-out $(CPP_FILTER), $(wildcard *.cpp))
CPP_HEADERS = $(wildcard *.h)
CPP_OBJECTS = $(CPP_SOURCES:.cpp=.o)

print_vars:
	@echo "Target $(CPP_TARGET)"
	@echo "Source Files $(CPP_SOURCES)"
	@echo "Header Files $(CPP_HEADERS)"
	@echo "Object Files $(CPP_OBJECTS)"

main: $(CPP_TARGET) $(CPP_OBJECTS) $(CPP_HEADERS)
	$(GCC) $(GCC_FLAGS) main.cpp $(CPP_OBJECTS) -o $@

%.out: %
	./$< > $@
	
clean:
	rm -f main $(CPP_OBJECTS)

