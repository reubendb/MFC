#!/bin/bash

# Color ANSI Escape Sequences
FG_RED='\033[0;31m'
FG_NONE='\033[0m'

./_install.sh "$@"

exit_code=$?

if [ $exit_code -ne 0 ]; then
    echo -e "$FG_RED\n"
    echo "==========================================================================="
    echo "| Fatal Error:                                                            |"
    echo "| ---> Exit Code ($exit_code).                                             "
    echo "| ---> Take a look at the output above.                                   |"
    echo "| ---> Some stages output logs to log/.                                   |"
    echo "| ---> Please ensure your build environment is adequate (view README.md). |"
    echo "==========================================================================="
    echo -e "\n$FG_NONE"
fi

exit $exit_code
