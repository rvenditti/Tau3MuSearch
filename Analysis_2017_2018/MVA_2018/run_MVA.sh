#!/bin/bash
#export ROOTSYS=/lustrehome/fsimone/root/
#source $ROOTSYS/bin/thisroot.sh 

echo "run TMVA from ROOT > 6.14 (CMSSW_10_6_8)"

echo "root -l -b MVA_code_2018.cpp\(\"Alowpt\"\)"
root -l -b MVA_code_2018.cpp\(\"Alowpt\"\)

echo "root -l -b MVA_code_2018.cpp\(\"Ahighpt\"\)"
root -l -b MVA_code_2018.cpp\(\"Ahighpt\"\)

echo "root -l -b MVA_code_2018.cpp\(\"Blowpt\"\)"
root -l -b MVA_code_2018.cpp\(\"Blowpt\"\)

echo "root -l -b MVA_code_2018.cpp\(\"Bhighpt\"\)"
root -l -b MVA_code_2018.cpp\(\"Bhighpt\"\)

echo "root -l -b MVA_code_2018.cpp\(\"C\"\)"
root -l -b MVA_code_2018.cpp\(\"C\"\)

#echo "root -l -b MVA_code_2018.cpp\(\"ABC\"\)"
#root -l -b MVA_code_2018.cpp\(\"ABC\"\)
