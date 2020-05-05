import sys
import os
import csv
import string
import datetime

# Define the parser
import argparse
parser = argparse.ArgumentParser(description="Options to give to the script")
# Positional arguments
parser.add_argument("dataset", type=str, choices=['data'], help="Specify 'data'")
parser.add_argument("anatype", type=str, choices=['tau3mu', 'control'], help="Specify analysis type")
# Optional Arguments
parser.add_argument("--run", type=str, default='', choices=['2017B', '2017C', '2017D', '2017E', '2017F', '2017Ds', '2017B0', '2017Bp', '2017DsPhiPi', '2017DsPhiPi_old',  '2018A', '2018B', '2018C', '2018D', '2018Ds', 'ParkingA1', 'ParkingA2', 'ParkingA3', 'ParkingA4', 'ParkingA5', 'ParkingA6', 'ParkingB1', 'ParkingB2', 'ParkingB3', 'ParkingB4', 'ParkingB5', 'ParkingB5_a', 'ParkingB5_b', 'ParkingB6', 'ParkingC1', 'ParkingC2', 'ParkingC3', 'ParkingC4', 'ParkingC5', 'ParkingD1', 'ParkingD2', 'ParkingD3', 'ParkingD4', 'ParkingD5_a', 'ParkingD5_b', 'ParkingDs', 'ParkingB0', 'ParkingBp', 'ParkingDsPhiPi', 'ParkingA1_mini', 'ParkingA2_mini', 'ParkingA3_mini', 'ParkingA4_mini', 'ParkingA5_mini', 'ParkingA6_mini'], help="run in data")
parser.add_argument("--n", type=int, default=255, help="number of .root files per job")
args = parser.parse_args()

#prepare output filename  and option string

if args.dataset == 'data':
   out_filename = 'AnalysedTree_'+args.dataset+'_'+args.run+'_'+args.anatype
   temp = '_'+args.anatype
   option_string = ' "'+args.dataset+temp.replace("_tau3mu","")+'" "'+args.run+'"'
else:
   out_filename = 'AnalysedTree_'+args.dataset+'_'+args.anatype
   temp = '_'+args.anatype
   option_string = ' "'+args.dataset+temp.replace("_tau3mu","")+'" "'

startTime = datetime.datetime.now().strftime("%Y%m%d_%H%M")

# Create target Directory if don't exist
if not os.path.exists(startTime):
    os.mkdir(startTime)
    print('Directory '+startTime+' created\n')
else:    
    print('Directory '+startTime+' already exists\n')

if args.anatype == 'tau3mu':
    if args.dataset == 'data' and args.run == '2017B':
       path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017B_AOD_v7/191120_150348'
    if args.dataset == 'data' and args.run == '2017C':
       path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017C_AOD_v7/191120_151537'
    if args.dataset == 'data' and args.run == '2017D':
       path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017D_AOD_v7/191120_160824'
    if args.dataset == 'data' and args.run == '2017E':
       path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017E_AOD_v7/191120_155344'
    if args.dataset == 'data' and args.run == '2017F':
       path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2017F_AOD_v7/191120_155750'
    if args.dataset == 'data' and args.run == '2017Ds':
        #path = '/lustre/cms/store/user/rosma/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_DsSignalOfficial_2017_v7/191120_145101/0000'
       path = '/lustre/cms/store/user/rosma/DsToTau_TauTo3Mu/SkimTau3Mu_DsToTauTo3Mu_2018_CMSSW_10_2_1_v7/200119_100659'
    if args.dataset == 'data' and args.run == '2017DsEtaMuNu':
       path = '/lustre/cms/store/user/caruta/DsEtaMuNu_EtaMuMuGamma/TreeMaker_DsEta_PHSP_v4/191121_123656/0000'
    if args.dataset == 'data' and args.run == '2017B0':
       path = '/lustre/cms/store/user/rosma/BdToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_BdSignalOfficial_2017_v7/191120_145848/0000'
    if args.dataset == 'data' and args.run == '2017Bp':
       path = '/lustre/cms/store/user/rosma/BuToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/TreeMaker_BuSignalOfficial_2017_v7/191120_145654/0000'
    if args.dataset == 'data' and args.run == '2018A':
       path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018A_AOD_v2/191122_153230'
    if args.dataset == 'data' and args.run == '2018B':
       path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018B_AOD_v2/191122_153359'
    if args.dataset == 'data' and args.run == '2018C':
       path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018C_AOD_v2/191122_153443'
    if args.dataset == 'data' and args.run == '2018D':
       path = '/lustre/cms/store/user/rosma/DoubleMuonLowMass/SkimTau3Mu_DoubleMuonLowMass_Run2018D_AOD_v2/191122_153513'
    if args.dataset == 'data' and args.run == '2018Ds':
       path = '/lustre/cms/store/user/rosma/DsToTau_TauTo3Mu/SkimTau3Mu_DsToTauTo3Mu_2018_v1/191128_165636/0000'
    if args.dataset == 'data' and args.run == 'ParkingA1':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH1/SkimTau3Mu_BParking_Run2018A_BPH1_v3/191202_141907'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingA1'
        #path = '/lustre/cms/store/user/caruta/TreeProva'
    if args.dataset == 'data' and args.run == 'ParkingA2':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH2/SkimTau3Mu_BParking_Run2018A_BPH2_v3/191202_141925'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingA2'
        #path = '/lustre/cms/store/user/caruta/ntuple_ParkingA2'
    if args.dataset == 'data' and args.run == 'ParkingA3':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH3/SkimTau3Mu_BParking_Run2018A_BPH3_v3/191202_141943'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingA3'
    if args.dataset == 'data' and args.run == 'ParkingA4':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH4/SkimTau3Mu_BParking_Run2018A_BPH4_v3/191202_142001'
        #path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingA4'
       path = '/lustre/cms/store/user/rosma/ParkingBPH4/SkimTau3Mu_BParking_Run2018A_BPH4_Mini_v4/'
    if args.dataset == 'data' and args.run == 'ParkingA5':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH5/SkimTau3Mu_BParking_Run2018A_BPH5_v3/191202_142022'
       #path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingA5'
       path = '/lustre/cms/store/user/rosma/ParkingBPH5/SkimTau3Mu_BParking_Run2018A_BPH5_Mini_v4-TightSkim'
    if args.dataset == 'data' and args.run == 'ParkingA6':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH6/SkimTau3Mu_BParking_Run2018A_BPH6_v3/191202_142045'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingA6'
    if args.dataset == 'data' and args.run == 'ParkingB1':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH1/SkimTau3Mu_BParking_Run2018B_BPH1_v3/191209_095556'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingB1'
    if args.dataset == 'data' and args.run == 'ParkingB2':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH2/SkimTau3Mu_BParking_Run2018B_BPH2_v3/191209_100209'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingB2'
    if args.dataset == 'data' and args.run == 'ParkingB3':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH3/SkimTau3Mu_BParking_Run2018B_BPH3_v3/191209_100230'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingB3'
    if args.dataset == 'data' and args.run == 'ParkingB4':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH4/SkimTau3Mu_BParking_Run2018B_BPH4_v3/191209_100249'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingB4'
    if args.dataset == 'data' and args.run == 'ParkingB5':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH5/SkimTau3Mu_BParking_Run2018B_BPH5_v3/191209_100309'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingB5'
    if args.dataset == 'data' and args.run == 'ParkingB6':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH6/SkimTau3Mu_BParking_Run2018B_BPH6_v3/191216_104254'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingB6'
    if args.dataset == 'data' and args.run == 'ParkingC1':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH1/SkimTau3Mu_BParking_Run2018C_BPH1/191128_175752'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingC1'
    if args.dataset == 'data' and args.run == 'ParkingC2':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH2/SkimTau3Mu_BParking_Run2018C_BPH2/191210_142930'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingC2'
    if args.dataset == 'data' and args.run == 'ParkingC3':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH3/SkimTau3Mu_BParking_Run2018C_BPH3/191210_142949'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingC3'
    if args.dataset == 'data' and args.run == 'ParkingC4':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH4/SkimTau3Mu_BParking_Run2018C_BPH4/191210_143009'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingC4'
    if args.dataset == 'data' and args.run == 'ParkingC5':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH5/SkimTau3Mu_BParking_Run2018C_BPH5/191210_143029'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingC5'
    if args.dataset == 'data' and args.run == 'ParkingD1':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH1/SkimTau3Mu_BParking_Run2018D_BPH1_v3-05May2019_v1/191218_113406'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingD1'
    if args.dataset == 'data' and args.run == 'ParkingD2':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH2/SkimTau3Mu_BParking_Run2018D_BPH2_v3-05May2019_v1/191223_101849'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingD2'
    if args.dataset == 'data' and args.run == 'ParkingD3':
       path = ''
    if args.dataset == 'data' and args.run == 'ParkingD4':
       path = ''
    if args.dataset == 'data' and args.run == 'ParkingD5_a':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH5/SkimTau3Mu_BParking_Run2018D_BPH5_v3-05May2019/191213_111841'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingD5_a'
    if args.dataset == 'data' and args.run == 'ParkingD5_b':
        #path = '/lustre/cms/store/user/rosma/ParkingBPH5/SkimTau3Mu_BParking_Run2018D_BPH5_v3-05May2019_v1_missLumi/200110_105015'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingD5_b'
    if args.dataset == 'data' and args.run == 'ParkingDs':
        #path = '/lustre/cms/store/user/rosma/DsToTau_TauTo3Mu/SkimTau3Mu_DsToTauTo3Mu_2018_CMSSW_10_6_1_v7/200218_121421/0000/'
        #path = '/lustre/cms/store/user/rosma/DsToTau_TauTo3Mu/SkimTau3Mu_DsToTauTo3Mu_2018_CMSSW_10_6_1_v7_noTrg/'
       path = '/lustre/cms/store/user/rosma/DsToTau_TauTo3Mu_March2020/SkimTau3Mu_DsToTauTo3Mu_2018_CMSSW_10_6_1_v7_noTrg_March2020/'
    if args.dataset == 'data' and args.run == 'ParkingB0':
        #path = '/lustre/cms/store/user/rosma/BdTau3Mu/SkimTau3Mu_BdToTauTo3Mu_2018_CMSSW_10_6_1_v7/200218_121505/0000/'
        #path = '/lustre/cms/store/user/rosma/BdTau3Mu/SkimTau3Mu_BdToTauTo3Mu_2018_CMSSW_10_6_1_v7_noTrg/200221_175010/0000/'
       path = '/lustre/cms/store/user/rosma/B0ToTau_TauTo3Mu/SkimTau3Mu_BdToTauTo3Mu_2018_CMSSW_10_6_1_v7_noTrg_March2020/200403_130342/0000/'
    if args.dataset == 'data' and args.run == 'ParkingBp':
        #path = '/lustre/cms/store/user/rosma/BuTau3Mu/SkimTau3Mu_BuToTauTo3Mu_2018_CMSSW_10_6_1_v7/200214_175607'
        #path = '/lustre/cms/store/user/rosma/BuTau3Mu/SkimTau3Mu_BuToTauTo3Mu_2018_CMSSW_10_6_1_v7_noTrg/200221_174938/0000/'
       path = ''
    if args.dataset == 'data' and args.run == 'ParkingA1_mini':
       path = '/lustre/cms/store/user/rosma/ParkingBPH1/SkimTau3Mu_BParking_Run2018A_BPH1_Mini_v2'
    if args.dataset == 'data' and args.run == 'ParkingA2_mini':
       path = '/lustre/cms/store/user/rosma/ParkingBPH2/SkimTau3Mu_BParking_Run2018A_BPH2_Mini_v2'
    if args.dataset == 'data' and args.run == 'ParkingA3_mini':
       path = '/lustre/cms/store/user/rosma/ParkingBPH3/SkimTau3Mu_BParking_Run2018A_BPH3_Mini_v3'
    if args.dataset == 'data' and args.run == 'ParkingA4_mini':
       path = '/lustre/cms/store/user/rosma/ParkingBPH4/SkimTau3Mu_BParking_Run2018A_BPH4_Mini_v2'
    if args.dataset == 'data' and args.run == 'ParkingA5_mini':
       path = '/lustre/cms/store/user/rosma/ParkingBPH5/SkimTau3Mu_BParking_Run2018A_BPH5_Mini_v2'
    if args.dataset == 'data' and args.run == 'ParkingA6_mini':
       path = '/lustre/cms/store/user/rosma/ParkingBPH6/SkimTau3Mu_BParking_Run2018A_BPH6_Mini_v2'

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
    if args.dataset == 'data' and args.run == '2018A':
       path = ''
    if args.dataset == 'data' and args.run == '2018B':
       path = ''
    if args.dataset == 'data' and args.run == '2018C':
       path = ''
    if args.dataset == 'data' and args.run == '2018D':
       path = ''
#   if args.dataset == 'data' and args.run == '2017ParkingDsPhiPi':
#      path = '/lustre/cms/store/user/rosma/DsToPhiMuMuPi_CMSSW_10_2_X_2018/SkimDsPhiPi_DsToPhiMuMuPi_2018_CMSSW_10_6_1_TrgSkim2017/200115_105419/0000/'
    if args.dataset == 'data' and args.run == 'ParkingA1':
       #path = '/lustre/cms/store/user/rosma/ParkingBPH1/SkimDsPhiPi_BParking_Run2018A_BPH1_v1/191202_144024'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingA1_Control'
    if args.dataset == 'data' and args.run == 'ParkingA2':
       #path = '/lustre/cms/store/user/rosma/ParkingBPH2/SkimDsPhiPi_BParking_Run2018A_BPH2_v1/191202_144042'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingA2_Control'
    if args.dataset == 'data' and args.run == 'ParkingA3':
       #path = '/lustre/cms/store/user/rosma/ParkingBPH3/SkimDsPhiPi_BParking_Run2018A_BPH3_v1/191202_144100'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingA3_Control'
    if args.dataset == 'data' and args.run == 'ParkingA4':
       #path = '/lustre/cms/store/user/rosma/ParkingBPH4/SkimDsPhiPi_BParking_Run2018A_BPH4_v1/191202_144119'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingA4_Control'
    if args.dataset == 'data' and args.run == 'ParkingA5':
       #path = '/lustre/cms/store/user/rosma/ParkingBPH5/SkimDsPhiPi_BParking_Run2018A_BPH5_v1/191202_144139'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingA5_Control'
    if args.dataset == 'data' and args.run == 'ParkingA6':
       #path = '/lustre/cms/store/user/rosma/ParkingBPH6/SkimDsPhiPi_BParking_Run2018A_BPH6_v1/191202_144157'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingA6_Control'
    if args.dataset == 'data' and args.run == 'ParkingB1':
       #path = '/lustre/cms/store/user/rosma/ParkingBPH1/SkimDsPhiPi_BParking_Run2018B_BPH1_v1/191209_104613'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingB1_Control'
    if args.dataset == 'data' and args.run == 'ParkingB2':
       #path = '/lustre/cms/store/user/rosma/ParkingBPH2/SkimDsPhiPi_BParking_Run2018B_BPH2_v1/191209_105824'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingB2_Control'
    if args.dataset == 'data' and args.run == 'ParkingB3':
       #path = '/lustre/cms/store/user/rosma/ParkingBPH3/SkimDsPhiPi_BParking_Run2018B_BPH3_v1/191209_121431'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingB3_Control'
    if args.dataset == 'data' and args.run == 'ParkingB4':
       #path = '/lustre/cms/store/user/rosma/ParkingBPH4/SkimDsPhiPi_BParking_Run2018B_BPH4_v1/191209_114258'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingB4_Control'
    if args.dataset == 'data' and args.run == 'ParkingB5_a':
       #path = '/lustre/cms/store/user/rosma/ParkingBPH5/SkimDsPhiPi_BParking_Run2018B_BPH5_v1/191209_115714'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingB5_a_Control'
    if args.dataset == 'data' and args.run == 'ParkingB5_b':
       #path = '/lustre/cms/store/user/rosma/ParkingBPH5/SkimDsPhiPi_BParking_Run2018B_BPH5_v1/200107_160004'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingB5_b_Control'
    if args.dataset == 'data' and args.run == 'ParkingB6':
       #path = '/lustre/cms/store/user/rosma/ParkingBPH6/SkimDsPhiPi_BParking_Run2018B_BPH6_v1/191209_093237'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingB6_Control'
    if args.dataset == 'data' and args.run == 'ParkingC1':
       path = ''
    if args.dataset == 'data' and args.run == 'ParkingC2':
       #path = '/lustre/cms/store/user/rosma/ParkingBPH2/SkimDsPhiPi_BParking_Run2018C_BPH2_v3/191210_145811'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingC2_Control'
    if args.dataset == 'data' and args.run == 'ParkingC3':
       #path = '/lustre/cms/store/user/rosma/ParkingBPH3/SkimDsPhiPi_BParking_Run2018C_BPH3_v3/191210_145828'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingC3_Control'
    if args.dataset == 'data' and args.run == 'ParkingC4':
       #path = '/lustre/cms/store/user/rosma/ParkingBPH4/SkimDsPhiPi_BParking_Run2018C_BPH4_v3/191210_145845'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingC4_Control'
    if args.dataset == 'data' and args.run == 'ParkingC5':
       #path = '/lustre/cms/store/user/rosma/ParkingBPH5/SkimDsPhiPi_BParking_Run2018C_BPH5_v3/191210_145903'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingC5_Control'
    if args.dataset == 'data' and args.run == 'ParkingD1':
       #path = '/lustre/cms/store/user/rosma/ParkingBPH1/SkimDsPhiPi_BParking_Run2018D_BPH1_v3-05May2019_v1/191223_102048'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingD1_Control'
    if args.dataset == 'data' and args.run == 'ParkingD2':
       #path = '/lustre/cms/store/user/rosma/ParkingBPH2/SkimDsPhiPi_BParking_Run2018D_BPH2_v3-05May2019_v1/200102_155630'
       path = '/lustre/cms/store/user/caruta/ntuple_Parking/ntuple_ParkingD2_Control'
    if args.dataset == 'data' and args.run == 'ParkingD3':
       path = ''
    if args.dataset == 'data' and args.run == 'ParkingD4':
       path = ''
    if args.dataset == 'data' and args.run == 'ParkingD5':
       path = ''
    if args.dataset == 'data' and args.run == 'ParkingDsPhiPi':
       path = '/lustre/cms/store/user/rosma/DsToPhiMuMuPi_CMSSW_10_2_X_2018/SkimDsPhiPi_DsToPhiMuMuPi_2018_CMSSW_10_2_1_v1_noTrg/200124_101254/0000'
       #path = '/lustre/cms/store/user/rosma/DsToPhiMuMuPi_CMSSW_10_2_X_2018/SkimDsPhiPi_DsToPhiMuMuPi_2018_CMSSW_10_6_1_v5/200114_101303/0000/'
    if args.dataset == 'data' and args.run == '2017DsPhiPi':
#      path = '/lustre/cms/store/user/rosma/DsToPhiMuMuPi_CMSSW_10_2_X_2018/SkimDsPhiPi_DsToPhiMuMuPi_2018_CMSSW_10_6_1_TrgSkim2017/200115_105419/0000/'
       path = '/lustre/cms/store/user/rosma/DsToPhiPi_ToMuMu_MuFilter_TuneCUEP8M1_13TeV-pythia8/DsToPhiPi_ToMuMu_MuFilter_SkimPhiPi_MC_2017_v7/191120_132753/0000/'

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

      cpp_filename = "Analysis_"+args.dataset+"_"+args.run+"_"+args.anatype+"_chunk"+str(file_index)+".cpp"
      with open(cpp_filename, "w") as out_file:
          for lb in buf:
              if lb == '        //AddFile_'+args.dataset+'_'+args.anatype+'\n':
                  #write group of files
                  out_file.write(chunk)
              elif lb == '            TString fileout = "AddOutput_'+args.dataset+'_'+args.anatype+'.root";\n':
                  #write output file name
                  out_file.write('        TString fileout = "'+out_filename+str(file_index)+'.root";\n')
              else: out_file.write(lb)

      #executable template
      with open("templates/launch_analysis_template.job", "r") as launch_infile:
          buf2 = launch_infile.readlines()

      launch_filename = "launch_analysis_"+args.dataset+"_"+args.run+"_"+args.anatype+"_"+str(file_index)+".job"
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

      condor_filename = "my_HTCondor_"+args.dataset+"_"+args.run+"_"+args.anatype+"_"+str(file_index)+".job"
      with open(startTime+"/"+condor_filename, "w") as myCondor_outfile:
          for lb3 in buf3:
              if lb3 == "Executable = launch_analysis_template.job\n":
                  myCondor_outfile.write("Executable = "+launch_filename+"\n")
              else: myCondor_outfile.write(lb3)

      #add lines to final script
      final_script.write("echo condor_submit "+condor_filename+" -name ettore\n")
      final_script.write("condor_submit "+condor_filename+" -name ettore\n")

final_script.close()
