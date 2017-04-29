# multirow
This repository contains the source code for the papers on multi-row intersection cuts by
**Ricardo Fukasawa**, **Laurent Poirrier** and **Álinson S. Xavier**. It also contains
the instances, the raw outputs and the scripts that were used to generate the tables 
that appear in the papers.

## Required tools and libraries

To produce the tables in the papers, the following tools and libraries were used.
Different versions may produce slightly different outputs.
Compiling the source code requires a modern C/C++ compiler, with support to C++11 features and OpenMP.

* GNU Make, version 3.81
* CMake, version 3.8
* GCC, the GNU Compiler Collection, version 6.3
* Ruby, version 1.9.3
* IBM® ILOG® CPLEX®, version 12.6

## Downloading and compiling

```
git clone https://github.com/iSoron/multirow.git
cd multirow
git submodule update
mkdir build
cd build
cmake ..
make
```

If CMake cannot find CPLEX, it may be necessary to modify the file `cmake/FindCPLEX.cmake`.

## Project structure

```
├── build             Directory that holds the compiled files
├── cmake             Custom CMake scripts
├── googletest        Google's C++ Test Framework
├── infinity          Source-code for the paper 'Intersection Cuts from the Infinity Norm'
│   ├── benchmark         Benchmark portion of the code, including scripts, instances and outputs
│   └── library           Actual code to compute infinity cuts
├── lifting           Source-code for the paper 'The (not so) Trivial Lifting in Two Dimensions'
│   ├── benchmark         Benchmark portion of the code, including scripts, instances and outputs
│   └── library           Actual code to perform the trivial lifting
├── multirow          Common code used by multiple papers
├── onerow            Source-code for the paper 'Intersection Cuts from Single-Row Relaxations'
│   ├── benchmark         Benchmark portion of the code, including scripts, instances and outputs
│   └── library           Actual code to generate wedge cuts
└── qxx               Auxiliary library by Laurent Poirrier
```
