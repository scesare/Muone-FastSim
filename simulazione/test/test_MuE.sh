#!/bin/sh

MAIN=MuE

ROOTINCDIR=`$ROOTSYS/bin/root-config --incdir`
ROOTLIBS=`$ROOTSYS/bin/root-config --cflags --libs`

code=../code


echo "Compiling ..."
rootcling -f MuEtreeDict.C -I${code} MuEtree.h MuEana.h MuEtreeLinkDef.h

FOPTIONS="-c -O2"
gfortran ${FOPTIONS} ${code}/summa.f
gfortran ${FOPTIONS} ${code}/hadr5n12.f
gfortran ${FOPTIONS} ${code}/muerad.f

OPTIONS="-O2 -Wall"
OPTIONS="-g -Wall"

g++ ${OPTIONS} -std=c++1y ${ROOTLIBS} -I${code} -I${ROOTINCDIR} ${code}/${MAIN}.cc MuEtreeDict.C ${code}/MuEtree.cc ${code}/Inputs.cc ${code}/Analysis.cc ${code}/FastSim.cc ${code}/Utils.cc ${code}/ElasticState.cc ${code}/ECAL.cc ${code}/EMShower.cc ${code}/RadialInterval.cc ${code}/GammaFunctionGenerator.cc ${code}/BaseNumericalRandomGenerator.cc ${code}/IncGamma.cc ${code}/dalpha.cc summa.o hadr5n12.o muerad.o -lgfortran -lX11 -o ${MAIN}.exe

#exit

OUTDIR=job_`date +"%y-%m-%d_%T"`

################################################################################
### FASTSIM input configuration ################################################
################################################################################
cat > fastSim.cfi <<!
<cfi>
0        # detector resolution model: 0=simplest-2par; 1=Antonio's 3 par
0        # 0:default; 1:MS only on X plane; 2: polar angle smearing
0.0425   # 0.042925 # total material thickness for model 0 (in X0) / 1.5cm Be = 0.0425 X0
0.02     # intr.ang.resol. for model 0 (in mrad) // optimal 10um/50cm=0.02 mrad
</cfi>
!
################################################################################

################################################################################
### ANALYSIS input configuration ################################################
################################################################################
cat > analysis.cfi <<!
<cfi>
1       # bool makeTree; 1(0) = do(not) produce the output Tree
1       # bool doTemplates; 1(0) = do (not) produce 2D template histos
35.     # max Theta for selection of events and template histograms
0       # 1(0): do(not) make angle correlation plots in fine bins
2       # parameterization of hadronic running: 0:pol2; 1:LL; 2:LLmod
0       # 0:nominal Lumi; 1:LowLumi; 
5       # range around expected average (as number of sigmas on each axis)
2       # number of divisions within one sigma interval
</cfi>
!
################################################################################

################################################################################
### MAIN input configuration ################################################
################################################################################
cat > input.cfi <<!
<cfi>
mysample.txt # string input_dirs_file; // file containing a list of directory paths where to look for NLO MC events
100                       # 10000000  long long n_events; // events to be processed (N>0:N; 0:all; <0: read histograms from existing results)10000000
results.root            # string histo_ifname; // path to input histo file (kinematical distributions when n_events<0)
${OUTDIR}               # string output_dir; // output directory name
fastSim.cfi             # string fastSim_ifname; // FastSim cfg file
analysis.cfi            # string analysis_ifname; // Analysis cfg file
</cfi>
!

echo "Running ..."

time ./${MAIN}.exe  input.cfi  > ${MAIN}.log 2>&1



mv input.cfi    ${OUTDIR}
mv fastSim.cfi  ${OUTDIR}
mv analysis.cfi ${OUTDIR}
mv ${MAIN}.log  ${OUTDIR} 

echo "Done."
