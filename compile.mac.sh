#!/bin/csh
echo "This script compiles Garffield"
echo "  options: compile, clean"
echo $1
if ($1 == "clean") then
echo "Cleaning Garffield"
cd lib/lammps/lib/reax
make -f Makefile.gfortran clean
cd ../../src
make clean-all
cd ../../optlist
make clean
cd ../pgapack/lib/wolf
rm -f *.o *.a
cd ../../../../src
make clean
cd ../
echo "Finished cleaning"
endif

if ($1 == "compile" ) then
echo "Compiling Garffield"
cd lib/lammps/lib/reax
make -f Makefile.gfortran
cd ../../src
make mac_mpi
make makelib
make -f Makefile.lib mac_mpi
cd ../../optlist
make
cd ../pgapack/source
make
cd ../../../src
make
cd ../
echo "Finished compiling garffield"
endif

