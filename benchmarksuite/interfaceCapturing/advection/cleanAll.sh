#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

cleanCase () {
    echo "cleaning $1"
    cd $1
    echo "y " | ./rmCases
    cd ..
}

cleanCase discUniFlow
cleanCase discUniFlowTet
cleanCase vortexShearedDisc
cleanCase vortexShearedDiscTet
cleanCase Leveque
cleanCase LevequeTet
cleanCase deformationSphere
cleanCase deformationSphereTet