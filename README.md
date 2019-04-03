# Tau3MuSearch
```bash
mkdir Tau3Mu
cd Tau3Mu
cmsrel CMSSW_9_4_4
cd CMSSW_9_4_4/src
cmsenv
git cms-init
git clone https://github.com/rvenditti/Tau3MuSearch
```
# Content

- MiniAna2017Tree module that makes trees storing variables related to the single reco muons and the triplet candidates

- Tau3MuSkim_cff.py is an example of configuration file selecting the triplet candidates (e.g. 3 loose muons with same vertex)

- run_Tau3MuSkim_cfg.py is the file that you have to "cmsRun" to generate the tree

- "Analysis" folder contains some root macros running on the trees
