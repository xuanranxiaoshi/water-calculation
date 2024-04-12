AVX512_SUPPORTED := $(shell (grep -q "avx512*" /proc/cpuinfo) && echo true || echo false)

ifeq ($(wildcard /opt/intel/oneapi/compiler/latest/bin/icpx),)
	CXX = g++
	COMPILER_MESSAGE = Using g++ as Compiler
	CXXFLAGS = -O3 -fopenmp -finline-functions -ftree-loop-optimize -ftree-vectorize -march=native
	ifeq ($(AVX512_SUPPORTED),true)
		CXXFLAGS += -mavx512f -mavx512dq -mavx512cd -mavx512bw -mavx512vl
	endif
else 
	CXX = /opt/intel/oneapi/compiler/latest/bin/icpx
	COMPILER_MESSAGE = Using Intel C++ Compiler as Compiler
	CXXFLAGS = -O3 -qopenmp -x HOST
	SETVAR_COMMAND = /opt/intel/oneapi/setvars.sh
endif

EXECUTABLE = calculate
SOURCE = calculate.cpp

all: $(EXECUTABLE)

$(EXECUTABLE): $(SOURCE)
	@echo $(COMPILER_MESSAGE)
	$(shell $(SETVAR_COMMAND))
	$(CXX) $(CXXFLAGS) $^ -o $@

profile: $(SOURCE)
	@echo "Compiling with profiling flags"
	$(shell $(SETVAR_COMMAND))
	$(CXX) -pg $(CXXFLAGS) $^ -o $(EXECUTABLE)_profiled

.PHONY: clean
clean:
	rm -f $(EXECUTABLE) $(EXECUTABLE)_profiled