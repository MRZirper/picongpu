#!/bin/bash                                                                                                                                                                                                                                                                    
#                                                                                                                                                                                                                                                                              
# This file is part of PIConGPU.                                                                                                                                                                                                                                               
# Copyright 2023 PIConGPU contributors                                                                                                                                                                                                                                         
# Authors: Hannes Wolf, Mika Soren Voss, Rene Widera                                                                                                                                                                                                                                       
# License: GPLv3+                                                                                                                                                                                                                                                              
#                                                                                                                                                                                                                                                                              
help()
{
  echo "Validate CurrentDeposition output data."
  echo "The test is evaluating the current density field with the corresponding analytic solution."
  echo ""
  echo "Usage:"
  echo "    validate.sh [-d dataPath] [inputSetPath]"
  echo ""
  echo "  -d | --data dataPath                 - path to simulation output data"
  echo "                                         Default: inputPath/simOutput/simData_%T.h5"
  echo "  -h | --help                          - show help"
  echo ""
  echo "  inputSetPath                         - path to directory of python test files"
  echo "                                         Default: ./lib/python/test/CurrentDeposition"
}

# options may be followed by
# - one colon to indicate they has a required argument
OPTS=`getopt -o d:h -l data:,help -- "$@"`
if [ $? != 0 ] ; then
    # something went wrong, getopt will put out an error message for us
    exit 1
fi

eval set -- "$OPTS"

# parser
while true ; do
    case "$1" in
        -d|--data)
            dataPath=$2
            shift
            ;;
        -h|--help)
            echo -e "$(help)"
            shift
            exit 0
            ;;
        --) shift; break;;
    esac
    shift
done

dataPath=$1
inputSetPath=$2

echo $dataPath

MAINTEST="./lib/python/test/CurrentDeposition"

if [ -z "$dataPath" ] ; then
    dataPath=$0/../simOutput/simData_%T.h5
else
    dataPath=$1
fi

echo $dataPath
python $MAINTEST/MainTest.py $dataPath
exit $?
