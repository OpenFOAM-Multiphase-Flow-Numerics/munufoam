#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

./circle/runAll.sh
./circleTet/runAll.sh
./sphere/runAll.sh
./sphereTet/runAll.sh
