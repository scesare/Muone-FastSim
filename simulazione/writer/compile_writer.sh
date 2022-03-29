#!/bin/sh

MAIN=write_MuE_MCevents

ROOTINCDIR=`$ROOTSYS/bin/root-config --incdir`
ROOTLIBS=`$ROOTSYS/bin/root-config --cflags --libs`

interface=../code

echo "Compiling ..."

rootcling -f MuEtreeWDict.C -I${interface} MuEtree.h MuEtreeWLinkDef.h

g++ -Wall -I${interface} -I${ROOTINCDIR} ${MAIN}.cc MuEtreeWDict.C ${ROOTLIBS} -lX11 -o ${MAIN}.exe
