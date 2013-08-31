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

4. Run `make` to compile. This will build `SCP3` and `SelectEW`. The latter reads the matrices stored `SCP3` in the `*.h5` files and computes the lowest eigenvalues and eigenvectors.

5. Switch to the folder where you want to keep the data and copy the queue submitting script from the `examples/` directory.

        cd /w2/igeorges/SCP3/TIP4P/hexamer/prism
        cp ~/SCP3/examples/scp3_prism.ge .
        cp ~/SCP3/examples/prism_T0K_1M_D50.dat .

6. Submit the job.

        qsub scp3_prism.ge

7. To obtain the lowest 10 eigenvalues of a Hamiltonian matrix, run

        /w2/igeorges/SCP3/build/SelectEW sM_2097152.h5 10

   The script `scp3_ew.ge` (respectively `scp3_ew.pbs` for PBS) in the
   `examples/` directory distributes diagonalization over the whole cluster. It
   submits individual jobs for each `*.h5` file in the current directory.

For `Makefile` options and technical details see the [Wiki](https://github.com/iong/SCP3/wiki).

TODO
====

- [ ] Actually create `scp3_ew.*` files.
- [ ] `SelectEW` outputs ASCII files.
- [ ] Describe `SelectEW` in Wiki
