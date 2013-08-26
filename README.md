SCP3

Folder structure

# How to build.

Compilers and compiler options are specified in the `Makefile`. All you need to do to compile the code is type


	make

at the command line. It will build the default target, `SCP3`, using the
default compiler, `gcc`. The `Makefile` will also try to find a `BLAS/LAPACK`
library by searching for Intel MKL, ACML or the reference BLAS/LAPACK, in this
order.

By default, the optimized binary is built. To use debugging, type

	make DEBUG=1

If you want to use the Intel Compilers, just type

	make COMPILER=intel

Changing the `BLAS/LAPACK` library is also possible

	make BLAS=acml

to force usage of AMD's ACML library or

	make BLAS=xyz

if you know that the BLAS library is in some `libxyz.so` file. For example, if the Goto2 optimized BLAS is installed, just type

	make BLAS=goto2

By default, the WHBB / HBB2-pol code is not included. To change that, type

	make WHBB=1

To summarize, a typical command line could look like this

	make COMPILER=intel WHBB=1 BLAS=mkl_parallel SCP3

to build SCP3 using the Intel compiler, support for WHBB and HBB2-pol and the parallel `BLAS/LAPACK` from Intel.

This particular `Makefile` is also able to keep the source code in a directory
and compile the code in a different one, keeping executables and object files
separate from the source. It can very useful if at some point you need a debug
build to find a problem, but don't want to recompile all the files in the
current directory. Recompiling Burkhardt's Sobol code takes a loooong time. If
you want to use this feature

	mkdir ~/build/SCP3
	cd ~/build/SCP3
	ln -s ~/src/SCP3/Makefile .
	make

