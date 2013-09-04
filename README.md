Crash course

1. Switch to your home directory. This is not necessary at all, but it makes the
   explanation easier.

2. If on HPC, prepare the environment by running

        module add gcc/4.7.3 igeorges/all

   Vladimir already has this line in his `~/.bashrc`.
 
3. Clone the code from github, unless you already have it.

        git clone https://github.com/iong/SCP3.git

   The source code is now in `SCP3`.

4. Switch to the directory where you want to keep the executables and create a
   symbolic link to the Makefile in the source directory

        cd /w2/mandelsh/SCP3/build
        ln -s ~/SCP3/Makefile .

5. Run `make` to compile. This will build `SCP3` and `SelectEW`. The latter
   reads the matrices stored `SCP3` in the `*.h5` files and computes the lowest
   eigenvalues and eigenvectors.

5. Switch to the folder where you want to keep the data and copy the queue
   submitting script from the `examples/` directory.

        cd /w2/mandelsh/SCP3/TIP4P/hexamer/prism
        cp ~/SCP3/examples/scp3_prism.ge .
        cp ~/SCP3/examples/prism_T0K_1M_D50.dat .

6. Submit the job.

        qsub scp3_prism.ge

   `SCP3` will write the Hamiltonian matrix in binary HDF5 files in regular
   intervals, by default every 2^17 Sobol points.
   
   [HDF5](http://www.hdfgroup.org/) is a self describing file format meant for
   very large datasets that require very fast access. HDF5 supports compression
   and allows direct access to parts of the file without first parsing the
   entire contents.

7. To obtain the lowest 10 eigenvalues of a Hamiltonian matrix, run

        /w2/mandelsh/SCP3/build/eigensolver sM_2097152.h5 10

   The examples directory also includes a script to run `eigensolver` in the
   queuing system. I'd recommend to copy it at the top of the SCP3 hierarchy:

        cp ~/SCP3/examples/eigensolver.ge /w2/mandelsh/SCP3/

   The same script can be used for all `sM_*.h5` files:

        qsub -N sM_0131072.h5 /w2/mandelsh/SCP3/eigensolver.ge
        qsub -N sM_0262144.h5 /w2/mandelsh/SCP3/eigensolver.ge
        ...

   The `-N` option to `qsub` specifies the name of the job, i.e. the name which
   will be displayed by `qstat`. `eigensolver.ge` reads this name to find out
   which file to use. This way, a single script can be used for all Hamiltonian
   matrices for all structures for all potentials.

   Shortcut: with a single command one can submit a job for every `sM_*.h5`
   file in the current directory:

        for x in sM_*.h5; do qsub -N $x ../../../eigensolver.ge; done

8. I will always push my changes to <https://github.org/iong/SCP3>. To update
   your version of the code, type

        cd ~/SCP3
        git pull

The [Wiki](https://github.com/iong/SCP3/wiki) includes a detailed description of
* [How To Compile](https://github.com/iong/SCP3/wiki/How-to-Compile)
* [All `SCP3` and `eigensolver` options](https://github.com/iong/SCP3/wiki/Using-SCP3)
* [The `module` environment on HPC](https://github.com/iong/SCP3/wiki/Technicalities)
