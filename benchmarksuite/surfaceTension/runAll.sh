#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

./advectedCircle/runAll.sh
./advectedCircleTet/runAll.sh
./sinWave/runAll.sh
./sinWaveTet/runAll.sh
