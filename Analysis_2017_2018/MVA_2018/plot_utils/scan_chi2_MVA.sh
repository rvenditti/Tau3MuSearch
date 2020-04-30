#!/bin/bash
#export ROOTSYS=/lustrehome/fsimone/root/
#source $ROOTSYS/bin/thisroot.sh 

echo "run TMVA from ROOT > 6.14 (CMSSW_10_6_8)"

for i in 10000 5000 2500 1000 500 100 50 40 30 20 15 10 8 6 4
do
    COMMAND='TString chi2cut = ''"fv_nC<'$i'";'
    echo $COMMAND
    sed -i -e "/TString chi2cut/s/.*/$COMMAND/" T3M_common.h

    #echo "root -l -b MVA_code_2018.cpp\(\"A\"\)"
    #root -l -b MVA_code_2018.cpp\(\"A\"\)

    #echo "root -l -b MVA_code_2018.cpp\(\"B\"\)"
    #root -l -b MVA_code_2018.cpp\(\"B\"\)

    #echo "root -l -b MVA_code_2018.cpp\(\"C\"\)"
    #root -l -b MVA_code_2018.cpp\(\"C\"\)

    ##echo "root -l -b MVA_code_2018.cpp\(\"ABC\"\)"
    ##root -l -b MVA_code_2018.cpp\(\"ABC\"\)

    root -l -g evaluate_fillbdt.C
    root -l -g BDT_optimal_cut.C
    root -l -g evaluate_fillmass.C
done
