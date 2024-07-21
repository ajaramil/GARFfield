#!/bin/csh 
echo "This script compiles Garffield"
echo "  options: <compile, clean> arch [garffield]"
echo $1 $2
# arch=ion,linux
# clean Garffield and all dependencies
if ($1 == "clean") then
  echo "Cleaning Garffield"
  # reaxlib
  cd lib/lammps/lib/reax
  make -f Makefile.gfortran clean
  cd ../../src
  # optlist
  make clean-all
  cd ../../optlist
  make clean
  # pgapack
  cd ../pgapack/lib/$2
  rm -f *.o *.a
  # garffield
  cd ../../../../src
  make clean
  cd ../
  echo "Finished cleaning"
endif

# make Garffield and all dependencies
if ($1 == "compile" ) then
  echo "Compiling Garffield"
  # reaxlib
  cd lib/lammps/lib/reax
  make -f Makefile.gfortran
  # lammps lib
  cd ../../src
  make -j 2 $2
  make makeshlib
  make -f Makefile.shlib $2
  # optlist
  cd ../../optlist
  make
  #pgapack
  cd ../pgapack/source
  make
  # garffield
  cd ../../../src
  make
  cd ../
  echo "Finished compiling garffield"
endif

# Only garffield code
if ($2 == "garffield") then
  cd src
  make clean
  make 
endif

