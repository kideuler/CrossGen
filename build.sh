#!/bin/bash

# Parse arguments: default ALL=false, set ALL=true if --all is provided
ALL=false
if [ $# -gt 0 ]; then
    for arg in "$@"; do
        if [ "$arg" = "--all" ]; then
            ALL=true
        fi
    done
fi

# build the executables
export CROSSGEN_DIR="$(pwd)"
rm -rf build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_OPENGL_VIEWER=ON ../
make -j
cd ..


