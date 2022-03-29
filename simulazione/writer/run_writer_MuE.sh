#!/bin/sh

MAIN=write_MuE_MCevents

INPUT=events-NLO.dat
OUTPUT=events-NLO.root

echo "Running ..."
./${MAIN}.exe  ${INPUT}  ${OUTPUT}  >  ${MAIN}.log 2>&1
echo "Done."
