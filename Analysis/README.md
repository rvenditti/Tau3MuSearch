The script createRunFile.py creates a folder in format YYYYMMDD_hhss and a bash script containing the commands for the sumbission of the analysis to the bath system.

usage: createRunFile.py [-h] [--run {2017B,2017C,2017D,2017E,2017F}] [--n N]
                        [--massregion {sgn,bkg}]
                        [--MCprocess {Ds,B0,Bp,DsPhiPi,MiniBias}]
                        {data,MC} {tau3mu,control}

Options to give to the script

positional arguments:
  {data,MC}             Specify if data or Monte Carlo
  {tau3mu,control}      Specify analysis type

optional arguments:
  -h, --help            show this help message and exit
  --run {2017B,2017C,2017D,2017E,2017F}
                        run in data
  --n N                 number of .root files per job
  --massregion {sgn,bkg}
                        Specify invariant mass range
  --MCprocess {Ds,B0,Bp,DsPhiPi,MiniBias}
                        process in Monte Carl0
                       
Depending on N, a certain number of files will be prepared in the YYYYMMDD_hhss folder:
- an executable file for the batch system (setting the environment, compile the .cpp, execute the code)
- a .cpp file for the given group of input files

Note: path to input ntuples are hardcoded in createRunFile.py.
