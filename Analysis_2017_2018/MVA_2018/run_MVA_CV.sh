#!/bin/bash
export ROOTSYS=/lustrehome/fsimone/root_build/
source $ROOTSYS/bin/thisroot.sh 

echo "set ROOT 6.20"

echo "root -l  MVA_code_2018_CV.cpp\(\"A\"\)"
root -l  MVA_code_2018_CV.cpp\(\"A\"\)

echo "root -l  MVA_code_2018_CV.cpp\(\"B\"\)"
root -l  MVA_code_2018_CV.cpp\(\"B\"\)

echo "root -l  MVA_code_2018_CV.cpp\(\"C\"\)"
root -l  MVA_code_2018_CV.cpp\(\"C\"\)
