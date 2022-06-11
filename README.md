# PPS
An Application to Generate Parametric Pseudo-Surfaces from Triangle Meshes

## INTRODUCTION

The code is organized as follows:

* bin            - subdirectory where the executables sample-pnt and sample-loop will be written to     
* data           - subdirectory where the example input files and empty output subdirectories are 
* doc            - subdirectory where doxygen documentation files are
* include        - subdirectory where include files for the libraries ppsfrompnt and ppsfromloop will be written to
* lib            - subdirectory where the lib file for the libraries ppsfrompnt and ppsfromloop will be written to
* CMakeLists.txt - cmake configuration file
* LICENSE.md     - license file
* README.md      - this file
* script         - subdirectory containing a shell script to execute the examples
* src            - subdirectory containing the source files of libraries and applications

In a nutshell, this code implements the surface representation described in the papers

Jean Gallier and Dianna Xu and Marcelo Siqueira  
Parametric pseudo-manifolds  
Differential Geometry and its Applications, 30(6), 2012, p. 702-736  
[DOI](https://doi.org/10.1016/j.difgeo.2012.09.002)

and

Marcelo Siqueira and Dianna Xu and Jean Gallier and Luis Gustavo Nonato and Dimas Martínez Morera and Luiz Velho  
A new construction of smooth surfaces from triangle meshes using parametric pseudo-manifolds  
Computers & Graphics, 33(3), 2009, p. 331-340  
[DOI](https://doi.org/10.1016/j.cag.2009.03.017)  

## INSTALLATION

You need cmake version 3.14 or higher and a C++ compiler that supports standard 17.  
To build the code, open a terminal, enter directory PPS, and run cmake as follows:  

```
cmake -S . -B build  
cmake --build build --config Release  
cmake --install build --prefix <full to the PPS directory or where else you want to have the executables>
```

If you have doxygen installed on your machine, the build will try to generate documentation too.

## EXAMPLES

You can run the executables 'sample-pnt' and 'sample-loop' on the data files in subdirectory

```
  <path to directory PPS>/examples
```

by executing the shell script 'run.sh' in subdirectory  

```
  <path to directory PPS>/script
```

##  LAST UPDATE

June 10, 2022

## CONTACT

If you run into trouble compiling or using the library, please email me at:

mfsiqueira@gmail.com

Have fun!

Marcelo Siqueira