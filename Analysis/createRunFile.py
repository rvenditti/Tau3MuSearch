import sys
import os
import csv
import string
import datetime

# Define the parser
import argparse
parser = argparse.ArgumentParser(description="Options to give to the script")
# Positional arguments
parser.add_argument("dataset", type=str, choices=['data', 'MC'], help="Specify if data or Monte Carlo")
parser.add_argument("anatype", type=str, choices=['tau3mu', 'control'], help="Specify analysis type")
# Optional Arguments
parser.add_argument("--run", type=str, default='', choices=['2017B', '2017C', '2017D', '2017E', '2017F'], help="run in data")
parser.add_argument("--n", type=int, default=255, help="number of .root files per job")
parser.add_argument("--massregion", type=str, default='', choices=['sgn', 'bkg'], help="Specify invariant mass range")
parser.add_argument("--MCprocess", type=str, default='', choices=['Ds', 'B0', 'Bp', 'DsPhiPi', 'MiniBias'], help="process in Monte Carlo")
args = parser.parse_args()

#prepare output filename  and option string
if args.massregion:
   args.massregion = '_'+args.massregion
if args.dataset == 'data':
   out_filename = 'AnalysedTree_'+args.dataset+args.massregion+'_'+args.run+'_'+args.anatype
   temp = '_'+args.anatype
   option_string = ' "'+args.dataset+temp.replace("_tau3mu","")+args.massregion+'" "'+args.run+'"'
else:
   out_filename = 'AnalysedTree_'+args.dataset+args.massregion+'_'+args.MCprocess+'_'+args.anatype
   temp = '_'+args.anatype
   option_string = ' "'+args.dataset+temp.replace("_tau3mu","")+args.massregion+'" "'+args.MCprocess+'"'

startTime = datetime.datetime.now().strftime("%Y%m%d_%H%M")

# Create target Directory if don't exist
if not os.path.exists(startTime):
    os.mkdir(startTime)
    print('Directory '+startTime+' created\n')
else:    
    print('Directory '+startTime+' already exists\n')

if args.anatype == 'tau3mu':
   if args.dataset == 'data' and args.run == '2017B':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017B_AOD_v3/190807_095259'
   if args.dataset == 'data' and args.run == '2017C':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017C_AOD_v3/190807_214813'
   if args.dataset == 'data' and args.run == '2017D':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017D_AOD_v3_newSplit'
   if args.dataset == 'data' and args.run == '2017E':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017E_AOD_v3_newSplit'
   if args.dataset == 'data' and args.run == '2017F':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017F_AOD_v3_newSplit'
if args.anatype == 'control':
   if args.dataset == 'data' and args.run == '2017B':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimPhiPi_DoubleMuonLowMass_Run2017B_31July_TM/190807_091820'
   if args.dataset == 'data' and args.run == '2017C':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimPhiPi_DoubleMuonLowMass_Run2017C_31July_TM/190807_214658'
   if args.dataset == 'data' and args.run == '2017D':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimPhiPi_DoubleMuonLowMass_Run2017D_31July_TM_newSplit/190829_123749'
   if args.dataset == 'data' and args.run == '2017E':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimPhiPi_DoubleMuonLowMass_Run2017E_31July_TM_newSplit/190902_083458'
   if args.dataset == 'data' and args.run == '2017F':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimPhiPi_DoubleMuonLowMass_Run2017F_19Aug_newSplit/190819_222432'
if args.dataset == 'MC' and args.MCprocess == 'Ds':
   path = ''
if args.dataset == 'MC' and args.MCprocess == 'B0':
   path = '/lustre/cms/store/user/rosma/BdToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_BdSignalOfficial_2017_v3_newSplit/190905_085548'
if args.dataset == 'MC' and args.MCprocess == 'Bp':
   path = '/lustre/cms/store/user/rosma/BuToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_BuSignalOfficial_2017_v3_newSplit/'
if args.dataset == 'MC' and args.MCprocess == 'DsPhiPi':
   path = '/lustre/cms/store/user/fsimone/histoSkim_twoMuonTrack_MCDsPhiPi'

#generating the list of all .root files in given directory and subdirectories
fileList = []
for r, d, f in os.walk(path): # r=root, d=directories, f = files
    for file in f:
        if '.root' in file:
            fileList.append(os.path.join(r, file))

#prepare final script
final_script = open("submit_analysis_"+startTime+".sh", "w")
final_script.write("#!/bin/bash\n")
final_script.write("chmod 777 -R *\n")
final_script.write("cd "+startTime+"\n")

#loop to generate one .cpp+executable+batch system conf file for each group of "n" files
n_chunk = len(fileList)//args.n
print('Number of files is {0:2d}'.format(len(fileList)))
print('Number of jobs is {0:2d}'.format(n_chunk+1))
for file_index in range(n_chunk+1):
      chunk = '' 
      for idx, l in enumerate(fileList):
         if idx < args.n*(file_index+1) and idx >= args.n*file_index:
             l = l.rstrip()
             l = '        chain->AddFile("{}");\n'.format(l)
             chunk = chunk + l

      #analysis.cpp template
      with open("templates/Analysis_template.cpp", "r") as in_file:
          buf = in_file.readlines()

      cpp_filename = "Analysis_"+args.dataset+args.massregion+"_"+args.run+args.MCprocess+"_"+args.anatype+"_chunk"+str(file_index)+".cpp"
      with open(cpp_filename, "w") as out_file:
          for lb in buf:
              if lb == '        //AddFile_'+args.dataset+args.MCprocess+'_'+args.anatype+'\n':
                  #write group of files
                  out_file.write(chunk)
              elif lb == '            TString fileout = "AddOutput_'+args.dataset+args.MCprocess+'_'+args.anatype+'.root";\n':
                  #write output file name
                  out_file.write('        TString fileout = "'+out_filename+str(file_index)+'.root";\n')
              else: out_file.write(lb)

      #executable template
      with open("templates/launch_analysis_template.job", "r") as launch_infile:
          buf2 = launch_infile.readlines()

      launch_filename = "launch_analysis_"+args.dataset+"_"+args.run+args.MCprocess+"_"+args.anatype+"_"+str(file_index)+".job"
      with open(startTime+"/"+launch_filename, "w") as launch_outfile:
          for lb2 in buf2:
              if lb2 == "#compile\n":
                  launch_outfile.write("cd "+startTime+"\n")
                  launch_outfile.write("g++ -I $ROOTSYS/include ../"+cpp_filename+" `root-config --glibs` `root-config --libs` `root-config --cflags` -lTMVA -L $ROOTSYS/lib -o executable"+str(file_index)+"\n")
              elif lb2 == "#execute\n":
                  launch_outfile.write('./executable'+str(file_index)+option_string+'\n')
              else: launch_outfile.write(lb2)

      #myCondor template
      with open("templates/my_HTCondor_template.job", "r") as myCondor_infile:
          buf3 = myCondor_infile.readlines()

      condor_filename = "my_HTCondor_"+args.dataset+"_"+args.run+args.MCprocess+"_"+args.anatype+"_"+str(file_index)+".job"
      with open(startTime+"/"+condor_filename, "w") as myCondor_outfile:
          for lb3 in buf3:
              if lb3 == "Executable = launch_analysis_template.job\n":
                  myCondor_outfile.write("Executable = "+launch_filename+"\n")
              else: myCondor_outfile.write(lb3)

      #add lines to final script
      final_script.write("echo condor_submit "+condor_filename+" -name ettore\n")
      final_script.write("condor_submit "+condor_filename+" -name ettore\n")

final_script.close()
