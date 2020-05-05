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
parser.add_argument("--run", type=str, default='', choices=[         '2016B', '2016C', '2016D', '2016E', '2016F', '2016G', '2016Hv2', '2016Hv3', 
                                                                     '2017B', '2017C', '2017D', '2017E', '2017F', 
                                                            '2018A', '2018B', '2018C', '2018D'
                                                            ], help="run in data")
parser.add_argument("--n", type=int, default=255, help="number of .root files per job")
parser.add_argument("--massregion", type=str, default='', choices=['sgn', 'bkg'], help="Specify invariant mass range")
parser.add_argument("--MCprocess", type=str, default='', choices=['2016Ds', '2016B0', '2016Bp',
                                                                  '2017Ds', '2017B0', '2017Bp',
                                                                  '2018Ds', '2018B0', '2018Bp',
                                                                  'DsPhiPi', 'MiniBias'], help="process in Monte Carlo")
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
output_folder = '/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/'
code_folder = os.path.dirname(os.path.abspath(__file__))

if not os.path.exists(output_folder+startTime):
    os.mkdir(output_folder+startTime)
    print('Directory '+output_folder+startTime+' created\n')
else:    
    print('Directory '+output_folder+startTime+' already exists\n')

if args.anatype == 'tau3mu':
   if args.dataset == 'data' and args.run == '2016B':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run201B_AOD_v1/200217_093036'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run201B_AOD_v0/200104_124224'
   if args.dataset == 'data' and args.run == '2016C':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016C_AOD_v1/200214_182431'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016C_AOD_v0/200104_130806'
   if args.dataset == 'data' and args.run == '2016D':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016D_AOD_v1/200214_183128'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016D_AOD_v0/200104_130853'
   if args.dataset == 'data' and args.run == '2016E':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016E_AOD_v1/200214_183658'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016E_AOD_v0/200104_132747'
   if args.dataset == 'data' and args.run == '2016F':
      path = ''
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016F_AOD_v0/200104_133637'
   if args.dataset == 'data' and args.run == '2016G':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016G_AOD_v1bis/200214_184229'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016G_AOD_v4_IsoBS/200210_130222'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016G_AOD_v2/200204_160258'
   if args.dataset == 'data' and args.run == '2016Hv2':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016Hv2_v1/200214_184620'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016Hv2_v0/200104_145634'
   if args.dataset == 'data' and args.run == '2016Hv3':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016Hv3_v1/200214_184703'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2016Hv3_v0/200104_151137'

   if args.dataset == 'data' and args.run == '2017B':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017B_AOD_v8/200403_175023/'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017B_AOD_v7/191120_150348'
   if args.dataset == 'data' and args.run == '2017C':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017C_AOD_v8-bis/200406_135655/'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017C_AOD_v8/200215_184354/'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017C_AOD_v7/191120_151537'
   if args.dataset == 'data' and args.run == '2017D':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017D_AOD_v8/200406_133445/'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017D_AOD_v7/191120_160824'
   if args.dataset == 'data' and args.run == '2017E':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017E_AOD_v8/200406_134450/'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017E_AOD_v7/191120_155344'
   if args.dataset == 'data' and args.run == '2017F':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017F_AOD_v8/200406_135435/'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017F_AOD_v7/191120_155750'

   if args.dataset == 'data' and args.run == '2018A':
      path = ''
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018A_AOD_Mini_v1/200427_162535/' #MiniAOD
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018A_AOD_Mini/200421_184022'    #MiniAOD
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018A_AOD_v4/200217_120328'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018A_AOD_v3/200119_100145'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018A_AOD_v2/191122_153230/'
   if args.dataset == 'data' and args.run == '2018B':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018B_AOD_v5/200429_184137/' #new variables
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018B_AOD_v4/200217_120246'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018B_AOD_v2/191122_153359/'
   if args.dataset == 'data' and args.run == '2018C':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018C_AOD_v5/200429_184043/' #new variables
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018C_AOD_v4/200217_120028'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018C_AOD_v2/191122_153443/'
   if args.dataset == 'data' and args.run == '2018D':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018D_AOD_v5/200429_183916/' #new variables
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018D_AOD_v4/200217_112055'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018D_AOD_v2/191122_153513/'

if args.anatype == 'control':
   if args.dataset == 'data' and args.run == '2017B':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimPhiPi_DoubleMuonLowMass_Run2017B_v7/191120_143429'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimPhiPi_DoubleMuonLowMass_Run2017B_31July_TM/190807_091820'
   if args.dataset == 'data' and args.run == '2017C':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimPhiPi_DoubleMuonLowMass_Run2017C_v7/191120_144121'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimPhiPi_DoubleMuonLowMass_Run2017C_31July_TM/190807_214658'
   if args.dataset == 'data' and args.run == '2017D':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimPhiPi_DoubleMuonLowMass_Run2017D_v7/191120_144415'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimPhiPi_DoubleMuonLowMass_Run2017D_31July_TM_newSplit/190829_123749'
   if args.dataset == 'data' and args.run == '2017E':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimPhiPi_DoubleMuonLowMass_Run2017E_31July_v7/191120_144631'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimPhiPi_DoubleMuonLowMass_Run2017E_31July_TM_newSplit/190902_083458'
   if args.dataset == 'data' and args.run == '2017F':
      path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimPhiPi_DoubleMuonLowMass_Run2017F_v7/191120_144844'
      #path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimPhiPi_DoubleMuonLowMass_Run2017F_19Aug_newSplit/190819_222432'

if args.dataset == 'MC' and args.MCprocess == '2017Ds':
   path = '/lustre/cms/store/user/rosma/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_DsSignalOfficial_2017_v8/200403_140244/'
   #no HLT! #path = '/lustre/cms/store/user/fsimone/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_DsSignal_2017_noHLT_v2/191202_100003/'
   #path = '/lustre/cms/store/user/rosma/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_DsSignalOfficial_2017_v7/191120_145101'
   #path = '/lustre/cms/store/user/rosma/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_DsSignalOfficial_2017_v6/191106_163531'
   #path = '/lustre/cms/store/user/rosma/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_DsSignalOfficial_2017_v4_newSplit/190927_121422/0000/'
if args.dataset == 'MC' and args.MCprocess == '2017B0':
   path = '/lustre/cms/store/user/rosma/BdToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_BdSignalOfficial_2017_v8/200403_181258/'
   #path = '/lustre/cms/store/user/rosma/BdToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_BdSignalOfficial_2017_v7/191120_145848'
   #path = '/lustre/cms/store/user/rosma/BdToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_BdSignalOfficial_2017_v6/191106_163824'
   #path = '/lustre/cms/store/user/rosma/BdToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_BdSignalOfficial_2017_v3_newSplit/190905_085548'
if args.dataset == 'MC' and args.MCprocess == '2017Bp':
   path = '/lustre/cms/store/user/rosma/BuToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_BuSignalOfficial_2017_v8/200403_181414/'
   #path = '/lustre/cms/store/user/rosma/BuToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_BuSignalOfficial_2017_v7/191120_145654'
   #path = '/lustre/cms/store/user/rosma/BuToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_BuSignalOfficial_2017_v6/191107_075024'
   #path = '/lustre/cms/store/user/rosma/BuToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_BuSignalOfficial_2017_v3_newSplit/'
if args.dataset == 'MC' and args.MCprocess == 'DsPhiPi':
   path = '/lustre/cms/store/user/rosma/DsToPhiPi_ToMuMu_MuFilter_TuneCUEP8M1_13TeV-pythia8/DsToPhiPi_ToMuMu_MuFilter_SkimPhiPi_MC_2017_v7/191120_132753'
   #path = '/lustre/cms/store/user/rosma/DsToPhiPi_ToMuMu_MuFilter_TuneCUEP8M1_13TeV-pythia8/DsToPhiPi_ToMuMu_MuFilter_SkimPhiPi_MC_Official_2017_v1/190927_142524/0000/'
   #path = '/lustre/cms/store/user/fsimone/histoSkim_twoMuonTrack_MCDsPhiPi'

if args.dataset == 'MC' and args.MCprocess == '2016Ds':
      path = '/lustre/cms/store/user/rosma/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_DsSignal_2016_v2/200215_162405'
      #path = '/lustre/cms/store/user/rosma/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_DsSignal_2016_v1/200210_125927'
      #path = '/lustre/cms/store/user/rosma/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_DsSignal_2016_v0/200103_114559/'
if args.dataset == 'MC' and args.MCprocess == '2016B0':
      path = '/lustre/cms/store/user/rosma/BdToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_BdSignal_2016_v2/200215_162625/'
      #path = '/lustre/cms/store/user/rosma/BdToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_BdSignal_2016_v0/200112_114429/'
if args.dataset == 'MC' and args.MCprocess == '2016Bp':
      path = '/lustre/cms/store/user/rosma/BuToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_BuSignal_2016_v2/200215_162024/'
      #path = '/lustre/cms/store/user/rosma/BuToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_BuSignal_2016_v0/200112_114249/'

if args.dataset == 'MC' and args.MCprocess == '2018Ds':
      path = '/lustre/cms/store/user/rosma/DsToTau_TauTo3Mu_March2020/SkimTau3Mu_DsToTauTo3Mu_2018_CMSSW_10_2_1_March2020_v9/200414_092747/' #new variables
      #path = '/lustre/cms/store/user/rosma/DsToTau_TauTo3Mu_March2020/SkimTau3Mu_DsToTauTo3Mu_2018_CMSSW_10_2_1_March2020_v8/200330_080934'
      #path = '/lustre/cms/store/user/rosma/DsToTau_TauTo3Mu/SkimTau3Mu_DsToTauTo3Mu_2018_CMSSW_10_2_1_Jianv2_v8/200217_121547'
      #path = '/lustre/cms/store/user/rosma/DsToTau_TauTo3Mu/SkimTau3Mu_DsToTauTo3Mu_2018_CMSSW_10_2_1_v4/191228_143705'
      #path = '/lustre/cms/store/user/rosma/DsToTau_TauTo3Mu/SkimTau3Mu_DsToTauTo3Mu_2018_v1/191128_165636'
if args.dataset == 'MC' and args.MCprocess == '2018B0':
      path = '/lustre/cms/store/user/rosma/B0ToTau_TauTo3Mu/SkimTau3Mu_B0ToTauTo3Mu_2018_CMSSW_10_2_1_v9/200414_092736/' #new variables
      #path = '/lustre/cms/store/user/rosma/B0ToTau_TauTo3Mu/SkimTau3Mu_B0ToTauTo3Mu_2018_CMSSW_10_2_1_v8/200330_132303'
      #path = '/lustre/cms/store/user/rosma/BdTau3Mu/SkimTau3Mu_BdToTauTo3Mu_2018_CMSSW_10_2_1_v8/200217_122205'
      #path = '/lustre/cms/store/user/rosma/BdTau3Mu/SkimTau3Mu_BdToTauTo3Mu_2018_CMSSW_10_2_1_v6/200122_090425/0000/'
if args.dataset == 'MC' and args.MCprocess == '2018Bp':
      path = '/lustre/cms/store/user/rosma/BuTau3Mu/SkimTau3Mu_BuToTauTo3Mu_2018_CMSSW_10_2_1_v9/200430_082303/' #new variables
      #path = '/lustre/cms/store/user/rosma/BuTau3Mu/SkimTau3Mu_BuToTauTo3Mu_2018_CMSSW_10_2_1_v8/200217_122234'
      #path = '/lustre/cms/store/user/rosma/BuTau3Mu/SkimTau3Mu_BuToTauTo3Mu_2018_CMSSW_10_2_1_v6/200122_091603/0000/'

#generating the list of all .root files in given directory and subdirectories
fileList = []
for r, d, f in os.walk(path): # r=root, d=directories, f = files
    for file in f:
        if '.root' in file:
            fileList.append(os.path.join(r, file))

#prepare final script
final_script_path = output_folder+startTime+"/submit_analysis_"+startTime+".sh"
final_script = open(final_script_path, "w")
final_script.write("#!/bin/bash\n")
final_script.write("chmod 777 -R "+output_folder+startTime+"/*\n")
final_script.write("chmod 777 -R "+code_folder+"/*\n")
final_script.write("cd "+output_folder+startTime+"\n")

#loop to generate one .cpp+executable+batch system conf file for each group of "n" files
n_chunk = len(fileList)//args.n
print('using ntuples in '+path)
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
              elif lb == '        //OutFile_'+args.dataset+args.MCprocess+'_'+args.anatype+'\n':
                  #write output file name
                  out_file.write('        fileout = "'+out_filename+str(file_index)+'.root";\n')
              else: out_file.write(lb)

      #executable template
      with open("templates/launch_analysis_template.job", "r") as launch_infile:
          buf2 = launch_infile.readlines()

      launch_filename = "launch_analysis_"+args.dataset+"_"+args.run+args.MCprocess+"_"+args.anatype+"_"+str(file_index)+".job"
      with open(output_folder+startTime+"/"+launch_filename, "w") as launch_outfile:
          for lb2 in buf2:
              if lb2 == "#compile\n":
                  launch_outfile.write("cd "+output_folder+startTime+"\n")
                  launch_outfile.write("g++ -I $ROOTSYS/include "+code_folder+"/"+cpp_filename+" `root-config --glibs` `root-config --libs` `root-config --cflags` -lTMVA -L $ROOTSYS/lib -o executable"+str(file_index)+"\n")
              elif lb2 == "#execute\n":
                  launch_outfile.write('./executable'+str(file_index)+option_string+'\n')
              else: launch_outfile.write(lb2)

      #myCondor template
      with open("templates/my_HTCondor_template.job", "r") as myCondor_infile:
          buf3 = myCondor_infile.readlines()

      condor_filename = "my_HTCondor_"+args.dataset+"_"+args.run+args.MCprocess+"_"+args.anatype+"_"+str(file_index)+".job"
      with open(output_folder+startTime+"/"+condor_filename, "w") as myCondor_outfile:
          for lb3 in buf3:
              if lb3 == "Executable = launch_analysis_template.job\n":
                  myCondor_outfile.write("Executable = "+launch_filename+"\n")
              else: myCondor_outfile.write(lb3)

      #add lines to final script
      final_script.write("echo condor_submit "+condor_filename+" -name ettore\n")
      final_script.write("condor_submit "+condor_filename+" -name ettore\n")

final_script.write("cd "+code_folder+"\n")
final_script.close()
print("to submit analysis to batch system do:\nsource "+final_script_path+"\n")
