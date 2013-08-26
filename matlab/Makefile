all: ScrambledSobol

mex_%.o : %.cpp
	$(CXX) -O0 -g -c -o $@ -I$(MATLAB_ROOT)/extern/include $^

mex_%.o : %.c
	$(CC) -O0 -g -c -o $@ -I$(MATLAB_ROOT)/extern/include $^

ScrambledSobol: mex_ScrambledSobol.o mex_MatousekAffineOwen.o
	$(CXX) -O0 -g -o $@ $^  -L$(MATLAB_ROOT)/bin/maci64 \
		-L$(MATLAB_ROOT)/runtime/maci64 \
		-Wl,-rpath -Wl,$(MATLAB_ROOT)/bin/maci64 \
		-Wl,-rpath -Wl,$(MATLAB_ROOT)/runtime/maci64 \
		-lmx -lmex -lmat -lmwmclmcrrt


