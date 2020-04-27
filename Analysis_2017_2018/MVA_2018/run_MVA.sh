#!/bin/bash

echo "run TMVA from ROOT > 6.14 (CMSSW_10_6_8)"

echo "root -l MVA_code_2018.cpp\(\"A\"\)"
root -l MVA_code_2018.cpp\(\"A\"\)

echo "root -l MVA_code_2018.cpp\(\"B\"\)"
root -l MVA_code_2018.cpp\(\"B\"\)

echo "root -l MVA_code_2018.cpp\(\"C\"\)"
root -l MVA_code_2018.cpp\(\"C\"\)

echo "root -l MVA_code_2018.cpp\(\"ABC\"\)"
root -l MVA_code_2018.cpp\(\"ABC\"\)
