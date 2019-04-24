#!/bin/bash
if cd build ; then
   rm -rf *
   cmake ..
   make
   ./clang-blueprint
else
    echo "build folder not found"
fi
