#!/bin/sh
### Shell script to create a ROOT loadable GCC shared lib out of .cxx source code
###
### NvE 23-may-2000 UU-SAP Utrecht
# 
### Name of the produced shared libraries
lib1=ralice.so
lib2=ralice.dylib
#
### Some MAC specific settings
export MACOSX_DEPLOYMENT_TARGET=10.3
unset LD_PREBIND
# 
### The option string for GCC shared lib compilation and linking ***
### For the GCC ROOT loadable shared lib the strict requirements are ***
### dropped to avoid many warnings from the rootcint generated code ***
gccroot="-c -g0 -ansi -pedantic -Wall -Wno-long-long -I$ROOTSYS/include"
#
echo "lib = " $lib
echo "gccroot = " $gccroot 
#
### Go to the directory with the source files
cd ..
#
### Create the dictionary files
rootcint -f zzzralicedict.cxx -c RALICEHeaders.h RALICELinkDef.h
# 
### Compile and create the ROOT loadable shared library
#
# Compilation phase
g++ $gccroot *.cxx   
#
# Creating ralice.so library 
g++ -v -bundle -undefined dynamic_lookup -o $lib1 *.o
#
# Creating ralice.dylib library 
g++ -v -dynamiclib -undefined dynamic_lookup -single_module -o $lib2 *.o
### On some systems the following extra "-read_only_relocs" flag might be needed
# g++ -v -dynamiclib -undefined dynamic_lookup -single_module -read_only_relocs -o $lib2 *.o
#
rm zzzralicedict.*
rm *.o
# 
### Move the created lib to the scripts directory and go there
mv $lib1 scripts
mv $lib2 scripts
cd scripts
#
echo ' ' 
echo '*** macgcclib done. Results in ' $lib1 ' and ' $lib2 
