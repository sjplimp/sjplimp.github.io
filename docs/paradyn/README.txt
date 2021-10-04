c ParaDyn - Parallel DYNAMO - molecular dynamics with EAM potentials
c 
c Authored by Steve Plimpton
c   (505) 845-7873, sjplimp@cs.sandia.gov
c   Dept 9221, MS 1111, Sandia National Labs, Albuquerque, NM  87185-1111
c
c Based on the serial DYNAMO code authored by
c   Stephen Foiles (foiles@ca.sandia.gov), Sandia National Labs, Livermore, CA
c   Murray Daw (daw@hubcap.clemson.edu), Clemson University
c
c See the README file for more information
c

ParaDyn 1998 - Release 1 (April 1998)

This distribution contains the ParaDyn code and related potential and
test files.  ParaDyn is a parallel molecular dynamics code which uses
embedded atom method (EAM) potentials to model metals and metal
alloys.  It is written in F77 and C and performs message-passing via
MPI calls.  Thus it should be portable to virtually any parallel or
single-processor machine.

All of ParaDyn's features are part of a newer parallel molecular
dynamics code called LAMMPS, which has a multitude of additional
features.  LAMMPS is generally faster than ParaDyn for reasonably
large problems.  ParaDyn is no longer under active development.  The
LAMMPS WWW Site at http://www.cs.sandia.gov/~sjplimp/lammps.html has
more details.

ParaDyn is an adaptation of the serial DYNAMO code written by Stephen
Foiles and Murray Daw.  My thanks to both Stephen and Murray for
assistance at various stages of this project.  Similar to DYNAMO,
ParaDyn's features include NPT dynamics, a choice of
temperature/pressure controls and boundary conditions, atom and region
constraints, and options for dynamics or energy minimization.  For
those who have used DYNAMO, one difference is that for portability
ParaDyn is run with a script file of input commands, rather than a
Fortran namelist deck.

The parallel techniques used in ParaDyn are discussed in the paper:

S. J. Plimpton and B. A. Hendrickson, "Parallel Molecular Dynamics
With the Embedded Atom Method", in Materials Theory and Modelling,
edited by J. Broughton, P. Bristowe, and J. Newsam, (MRS Proceedings
291, Pittsburgh, PA, 1993), p. 37.

The code and a postscript copy of the paper are available from me via
e-mail or on the Web at http://www.cs.sandia.gov/~sjplimp/main.html

Please send me e-mail if you want to be notified of any upgrades or
significant bug-fixes to ParaDyn.

Feel free to contact me regarding ParaDyn if you

(a) have specific questions about how to use it,
(b) encounter problems when running it or find any bugs,
(c) port it to a new machine,
(d) succeed in using or modifying it to solve some interesting problem.

Steve Plimpton

----------------------------------------------------------------------------

When you uncompress (gunzip paradyn.tar.gz) and untar (tar xvf
paradyn.tar) the distribution, you should have 2 files and 4
directories: src, potentials, tools, examples.

README          this file
OVERVIEW	explanation of all ParaDyn input commands

***** src directory:

*.f,*.c,*.h	ParaDyn source files
Crib		variable documentation for ParaDyn
Makefile	top-level generic Makefile
Makefile.xxx    low-level Makefiles for various machines

***** potentials directory:

cuu3,niu3,...	potential files for various elements in DYNAMO format

***** tools directory:

dyn2para.f	serial utility to convert DYNAMO coordinate 
			files to ParaDyn atom files - use for converting
			DYNAMO output files to ParaDyn input files -
			see header of file for usage syntax
randomize.f	serial utility to randomize atom order in an atom file
Makefile	Makefile to make tools

***** examples directory:

Contains input and output files for 4 benchmark calculations you can
use to test ParaDyn.  See the "Testing ParaDyn" section below.

----------------------------------------------------------------------------

Making ParaDyn:

The src directory contains the F77 and C source files for ParaDyn as
well as a top-level Makefile and examples of lower-level
machine-specific Makefile.xxx files.

Makefile is setup to create multiple versions of ParaDyn for various
machines.  Type "make" to see a list of supported targets.  Typing
"make target" will produce an executable named "pd_target".

The top-level Makefile works on most Unix platforms, but may be
incompatible with some, since the "make" command is not very standard.
If it is available on your system, an alternative is the Gnu make
command "gmake", which is standard.

When a particular target is made, it will create an Obj_target
sub-directory to store machine-specific *.o files.  This is to prevent
*.o files generated for different machines from getting mixed up.

If a Makefile.xxx file exists for your machine, you will likely need
to edit the section with system-specific paths, compiler options, MPI
libraries, etc.  To make ParaDyn for a new machine, simply copy one of
the Makefile.xxx files to use as a template, giving it a new suffix.
Then add the new suffix to the top-level Makefile in the appropriate
places.

ParaDyn is a parallel code that will run on any number of processors,
but also runs fine in serial on a single-processor - e.g. Cray vector
machines or any Unix workstation or (possibly) even a Wintel PC.  To
compile for a single processor, you have 2 choices.  If you have MPI
installed, you can compile and link as usual -- see the Makefile.sgi
file for example.  This will allow you to run on multiple (virtual)
processors of your workstation - e.g. using the mpirun command.  Or
you can compile without MPI, using the Makefile.serial file.  This
requires that you first make a library in the STUBS sub-directory that
provides hooks for the MPI routines called by ParaDyn.  Just type
"make" in the STUBS directory to create this library.

Note that ParaDyn is not optimized for the vector and RISC processing
single processors often support.  If you plan to run large
calculations on a vector machine such as the Cray Y-MP or C90 you may
be better off getting the original DYNAMO code from Stephen Foiles
which is optimized for vector processing.

There are 2 F77 compiler flags that include specific features in the
code.  If you add -DSTRESS to the F77FLAGS line, you will compute the
stresses on individual atoms which can be dumped to a file as desired.
This slows down the force computations by about 10%, so it is not
turned on by default.

If you have a Fortran compiler that supports the common "pointer"
extensions to F77, you can add -DDYNAMIC, which allows ParaDyn to set
array sizes at run-time.  The advantage is that you do not typically
need to recompile the code when you run different size problems on
different numbers of processors.  This flag also slows down the code
by about 10%, so you may not wish to use it.  Ideally you might run a
-DDYNAMIC version of the code while setting up or running a few
timesteps of your problem.  Then look at the "Array bounds" section at
the bottom of the log file to see the exact sizes needed, set the
parameters accordingly in paradyn.h, and recompile without the
-DDYNAMIC flag for the production runs of your problem.  If you change
the setting of either of these flags, make sure you recompile the
entire ParaDyn source, since the Makefile will not trigger this
recompilation automatically.  Also, if your compiler does not support
"pointer" statements you should delete memory.f from the list of *.f
files in the top-level Makefile.

There are other compiler flags you may want to experiment with for
optimal performance of the inter-processor communication routines in
ParaDyn.  See the discussion in the Running ParaDyn section below.

----------------------------------------------------------------------------

Testing ParaDyn:

In the examples directory, there are 4 *.in benchmark input files.
Two are small enough to run quickly on a single workstation processor
(cu_bulk.in and min.in).  All 4 can be run on a parallel machine.  You
can compare the results from running on your machine with the example
output files that are in the test directory for accuracy.  You should
get identical answers on any number of processors.  Files with a
example.xxx.P suffix are sample outputs from running the benchmarks on
P processors of machine xxx.  The test scripts will run all 4 examples
on different machines.

Here is how to run the benchmark tests:

(1) pd_sgi < cu_bulk.in

Runs 100 steps of a bulk Cu lattice of 500 atoms.
Creates cu_bulk.log as output.

(2) pd_sgi < gb_diff.in

Runs 3000 steps of a Cu grain boundary with 1720 atoms.
Needs gb_diff.atoms as additional input.
Creates gb_diff.log as output.
Creates gb_diff.diff as additional output
  so long as ParaDyn is compiled with correct diagnostic routine.

(3) pd_sgi < gb_random.in

Same as previous run except used randomized atom set with
  parallel method 2.
Runs 3000 steps of a Cu grain boundary with 1720 atoms.
Needs gb_diff.atoms_random as additional input
  (created from running gb_diff.atoms thru randomize).
Creates gb_diff.log_random as output.
Creates gb_diff.diff_random as additional output
  so long as ParaDyn is compiled with correct diagnostic routine.
Should get statistically similar answers to (2).

(4) pd_sgi < min.in

Minimizes the energy of a Cu impurity atom in a 500 atom Ni lattice.
Needs min.atoms as additional input.
Creates min.log and min.dump as output.

----------------------------------------------------------------------------

Running ParaDyn:

I have attempted to implement essentially all the features of DYANMO
in ParaDyn.  These features are accessed as commands in the input file
and are described in the OVERVIEW file.

The potential files used by ParaDyn are in the same format used by
DYNAMO for single elements or alloy systems.  A few of these are
provided in the potentials directory.  Stephen Foiles has a large
collection of these.

There are three input file commands which affect how fast a problem
runs on a given number of processors.  These are "parallel method",
"newton flag", and "neighbor method" and are described in the OVERVIEW
file.  When running a particular problem on a particular number of
processors, you may want to experiment with different values of these
three options to see what gives the fastest run time.  In principle
you should get identical answers using any combination of the three
options on any number of processors.  In practice, round-off errors
can cause slight differences and eventual divergence of dynamical
trajectories.

As discussed in the OVERVIEW file, when using parallel method = 2
(force-decomposition), for optimal performance you should run on P
processors where P = M*N and M is roughly equal to N.  If you run with
P a prime (or with widely differing factors), then your communication
costs will be closer to parallel method = 1.

Another issue affecting execution speed for parallel method = 2
(force-decomposition) is the ordering of atoms.  If the ordering is
regular either in the input data file ("read atoms") or as the lattice
is generated ("create atoms"), then force-decomposition will be more
likely to suffer from load-imbalance.  This can slow the code down in
parallel.  Thus you should first use the serial code "randomize" from
the tools directory to randomly permute the order of atoms in a "read
atoms" input file.  See the header of the randomize.f file for more
information.  Or you should use the "create atoms" command with a
non-zero 6th parameter - see the OVERVIEW file for details.  One
important note: two runs with different atom orderings should give
statistically similar results (e.g. thermodynamic averages), but will
typically not be identical.  This will happen if you use the "create
velocity" command, because it does not assign initial velocities the
same for both runs.  It will also happen if you use "read velocities"
on an input file of velocities that has not been permuted with exactly
the same re-ordering.

A final issue affecting execution speed is the communication routines
used in ParaDyn.  As discussed at the top of gs.c the default is to
use a recursive gather/scatter with blocking sends and receives.  This
is often the fastest option, but may cause some MPI implementations to
hang if you run big simulations.  This will typically happen before
the first thermodynamics print-out to the screen.  Sometimes this can
be controlled by setting appropriate environment variables that alter
MPI rules for when message buffering is done.  If the code hangs, you
can add a -DGS_IRECV switch to your CCFLAGS options in the
Makefile.xxx file and recompile gs.c.  This will use non-blocking
receives which require no buffer space and thus should not hang, but
is often a bit slower.  You can also experiement with -DGS_MPI which
will call the MPI_Allgatherv and MPI_Reduce_scatter routines directly.
In most MPI implementations these are slower than the recursive
gather/scatter routines I provide, but maybe you will be lucky!
Finally, if you have an optimized daxpy routine (e.g. BLAS routine) on
your machine, you can use the -DGS_AXPY switch to boost performance a
bit.  You must check that the syntax in gs.c for calling daxpy matches
your library routine.

There is also a -DSYNC switch you can compile the force.f file with,
which will do extra synchronization before calling communication
routines.  This will tend to slow down the code slightly, but will
give you a more accurate print-out of the time spent communicating.
Without this switch, load imbalance will typically be included in the
communication time; with this switch it should be included in the
force time.

If you get an error message when running ParaDyn about "boosting"
something, it means your arrays are not allocated large enough.  You
need to modify the appropriate parameter statement(s) at the top of
the param.h file and recompile.  Some of these errors are detected at
setup, others (like neighbor list overflow) may not occur until the
middle of a run.  When the latter happens the program will either
gracefully stop (if all processors incurred the same error) or hang.
Either way you should get an error message printed to the screen.  You
can also get an error message about running out of physical memory on
a processor.  If you can't set any parameters smaller and still run
your problem, then you need more processors!

A "boost" message should occur only rarely if you compiled with
-DDYNAMIC.  The most likely culprit will be the neighbor list arrays
For the -DDYNAMIC version you increase these sizes by setting the
"extra_neigh" parameter in paradyn.h to a larger value.  This is a
multiplier on the size for the neighbor list arrays that ParaDyn
estimates at run-time.  If you are running in static memory mode (no
-DDYNAMIC switch), you are more likely to get boost errors for
maxlocal, maxrow, or maxcol, all of which are functions of the number
of atoms N and the number of processors P.  Roughly speaking, you need
maxlocal >= N/P.  For atom-decomposition, you need maxrow >= N/P and
maxcol = N.  For force-decomposition, you need maxrow >= N/sqrt(P) and
maxcol >= N/sqrt(P).

I've tried to be pretty careful about detecting memory-overflow kinds
of errors.  If ParaDyn ever crashs or hangs without spitting out an
error message first due to algorithmic/parallelism problems (as
opposed to physics problems, like too big a timestep or putting 2
atoms on top of each other or the GS_IRECV problem discussed above),
it's probably a bug, so let me know about it.

To use ParaDyn to measure other properties of a system besides the
simple thermodynamic quantities (T,P,etc) it outputs to the screen and
log file you have 2 options.  You can simply dump atom positions
and/or velocities to disk and post-process the files to compute
desired quantities.  Or you write your own diagnostic routine to
compute the desired quantities on-the-fly as the simulation runs.  The
diagnostic.f file has an example routine that works in this fashion to
compute diffusion coefficients in a grain boundary system.  See the
top of the diagnostic.f file and the "diagnostic" command description
in the OVERVIEW file for more information.  If you don't want any
diagnostics you should compile with the diagnostic.hold file instead
of the diagnostic.f file.  Note that 2 of the example tests use the
diagnostic.f provided.  You can also add a user-defined force to
whatever atoms you want every time the force routine is called, by
adding a userforce routine in the file userforce.f.

There are two utilities in the tools directory.  For DYNAMO users
there is a dyn2para code that converts DYNAMO configuration files to
ParaDyn input files (atom and velocity dumps).  See the header of the
dyn2para.f file for the usage syntax.  The randomize utility was
discussed above.

----------------------------------------------------------------------------

Understanding ParaDyn:

If you wish to modify ParaDyn or understand it's inner workings you
may find the Crib file useful.  I'd like to think the source code is
so cleanly written as to be self-documenting (ha!).  The MRS paper is
the best high-level overview of what ParaDyn is doing as far as
parallelism.

If you can't understand something or want suggestions about modifying
some section of code, give me a call or send e-mail about it and I'll
explain (or obfuscate) things more fully.
