#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory

sphinx-build -b html ./doc ./docs/