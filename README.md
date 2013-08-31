Crash course

1. Switch to your home directory. This is not necessary at all, but make the
explanation easier.

2. If on HPC, prepare the environment by running

        module add igeorges/all

3. Clone the code from github, unless you already have it.

        git clone https://github.com/iong/SCP3.git

The source code is now in `SCP3`.
3. Switch to the directory where you want to keep the executables and create a
   symbolic link to the Makefile in the source directory

        cd /w2/igeorges/SCP3/build
        ln -s ~/SCP3/Makefile .

4. Run `make` to compile

For `Makefile` options and technical details see the [Wiki](https://github.com/iong/SCP3/wiki).
