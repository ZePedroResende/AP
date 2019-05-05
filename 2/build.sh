#!/bin/bash
if cd build ; then
   rm -rf *
   cmake ..
   make
   ./rb
else
    echo "build folder not found"
fi
