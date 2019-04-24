#!/bin/bash
if cd build ; then
   rm -rf *
   cmake ..
   make
else
    echo "build folder not found"
fi
