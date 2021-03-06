lifting
=======

This package contains the source code for the paper **The (not so) Trivial
Lifting in Two Dimensions**, by Ricardo Fukasawa, Laurent Poirrier and Álinson
S. Xavier.

Required Tools and Libraries
----------------------------

To produce the tables in the paper, the following tools and libraries were
used. Different versions may produce slightly different outputs.

- GNU Make, version 3.81
- CMake, version 3.7.2
- GCC, the GNU Compiler Collection, version 6.3.0
- Ruby, version 2.4.0
- IBM® ILOG® CPLEX®, version 12.6

Build instructions
------------------

1. Navigate to the folder `../build`
2. Run `cmake ..` followed by `make lifting-benchmark.run`
3. Two binaries (`lifting-benchmark.run` and `liblifting.a`) will be generated

Running the experiments
-----------------------

1. Build the project, following the instructions above.
2. Navigate to the folder `lifting/benchmark` and execute

        ./run_experiments.sh

3. Two CSV files will be generated inside the folder `lifting/benchmark/tables`,
	corresponding to the two tables that appear in the paper.

Modifying the instances
-----------------------

In order to run the experiments with a different set of instances,
the file `lifting/benchmark/instances/filtered/all.txt` should be modified.
Each line in this file describes the origin `f` and a lattice-free set `B`.
The set `B` is described by the coordinates of its vertices. Since the
benchmark code only deals with maximal lattice-free sets, it is also necessary
to specify the lattice-points that belong to each facet of `B`.
If `n` is the number of facets, `v` is an n-by-2 matrix of doubles corresponding
to the vertices and l is an n-by-2 matrix of doubles corresponding to the
lattice-points, then each line of the file should be written as

	f[0] f[1] n v[0][0] v[0][1] ... v[n-1][0] v[n-1][1] n l[0][0] l[0][1] ... l[n-1][0] l[n-1][1]

