#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

cleanCase () {
    echo "cleaning $1"
    cd $1
    echo "y " | ./rmCases
    cd ..
}

cleanCase advectedCircle
cleanCase advectedCircleTet
cleanCase sinWave
cleanCase sinWaveTet
