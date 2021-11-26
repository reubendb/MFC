#!/bin/bash

./_install.sh

exit_code=$?

if [ $exit_code -ne 0 ]; then
    echo 
    echo
    echo 
    echo 
    echo "========================================================================"
    echo "| Fatal Error"
    echo "| ---> Exit Code ($exit_code)"
    echo "| ---> Some stages output logs to log/"
    echo "| ---> Please ensure your build environment is adequate (view README.md)"
    echo "========================================================================"
    echo 
    echo
    echo 
    echo
fi

exit $exit_code