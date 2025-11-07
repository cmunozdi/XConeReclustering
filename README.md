# XConeReclustering

## Setting up the CMSSW framework for generating PFNano files and reprocessing them with the XCone algorithm.
```bash
cmsrel CMSSW_15_0_3_patch1
cd CMSSW_15_0_3_patch1/src
cmsenv
git cms-addpkg PhysicsTools/NanoAOD
scram b -j 9 
```

Getting btv production (just to prepare and send to crab the jobs):
```bash
git clone https://github.com/cms-btv-pog/btvnano-prod.git
scram b -j 9 
source /cvmfs/cms.cern.ch/common/crab-setup.sh prod # note: this is new w.r.t. 106X instructions
```

Getting XCone reclustering processing:
```bash
git clone https://github.com/cmunozdi/XConeReclustering.git
```
Make sure that:
```bash
echo $FASTJET_CONTRIB_BASE
echo $CPLUS_INCLUDE_PATH
echo $LD_LIBRARY_PATH
```
shows the path correctly. If not, check with ``scram tool info fastjet-contrib``, you should see something like:
```bash
Tool info as configured in location /afs/cern.ch/user/c/cmunozdi/CMSSW_15_0_3_patch1
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Name : fastjet-contrib
Version : 1.051-6b2cbf7b2399385490e165663710d5e0
Revision : 1
++++++++++++++++++++

FASTJET_CONTRIB_BASE=/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/fastjet-contrib/1.051-6b2cbf7b2399385490e165663710d5e0
INCLUDE=/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/fastjet-contrib/1.051-6b2cbf7b2399385490e165663710d5e0/include
LIB=fastjetcontribfragile
LIBDIR=/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/fastjet-contrib/1.051-6b2cbf7b2399385490e165663710d5e0/lib
```
and run something like this with the info comming from ``scram``:
```bash
export FASTJET_CONTRIB_BASE=/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/fastjet-contrib/1.051-6b2cbf7b2399385490e165663710d5e0
export CPLUS_INCLUDE_PATH=$FASTJET_CONTRIB_BASE/include:$CPLUS_INCLUDE_PATH
export LD_LIBRARY_PATH=$FASTJET_CONTRIB_BASE/lib:$LD_LIBRARY_PATH
```
Once you have set this, you can compile everithing with:
```bash
scram b -j 10
```

## Run locally XCone:
```bash
cmsenv
source env_xcone.sh
python3 ProcessNanoToBoostedTopQuarkWithXCone.py --input <inputfile or inputdirectory> --output <output root file> --isMC
```
Replace `<inputfile or inputdirectory>` with the path to your input file or directory, and `<output root file>` with the desired name for the output ROOT file, and write `--isMC` only for MC.
 

# Generate dictionary for custom Structs on the header
1. Create the dictionary file ``selection_helpers_BoostedTopQuark_LinkDef.h``:
```cpp
#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

// Declara las estructuras y funciones que deseas incluir en el diccionario
#pragma link C++ struct Lepton+;
#pragma link C++ struct JetAntikTReclus+;
#pragma link C++ struct XConeReclusteredJets+;
#pragma link C++ struct JetReclus+;
#pragma link C++ struct TopJetReclus+;
#pragma link C++ function triggerLepton+;
#pragma link C++ function CombineLeptons+;
#pragma link C++ function buildXConeJets+;
#pragma link C++ function Get_pTmiss+;
#pragma link C++ function initPlugin+;


#endif
```

2. Generate dictionary with ``rootcling``:
```bash
rootcling -f selection_helpers_BoostedTopQuark_Dict.cpp \
    -c selection_helpers_BoostedTopQuark.h selection_helpers_BoostedTopQuark_LinkDef.h
```

3. Compile to get the ``.so`` (remember doing ``source env_xcone.sh`` first):
```bash
g++ -shared -fPIC -o selection_helpers_BoostedTopQuark.so \
    selection_helpers_BoostedTopQuark.cpp selection_helpers_BoostedTopQuark_Dict.cpp \
    `root-config --cflags --libs`
```

## Merging JSON files and getting the Luminosity (processed Lumi)
1. Generate the individual JSON files:
```bash
for dir in DATA_samples/data_2023_Muon/crab_Muon*; do crab report -d "$dir"; done
```

2. Merging them:
```bash
mergeJSON.py --output DATA_samples/data_2023_Muon/combinedMuon2023processLumis.json $(find DATA_samples/data_2023_Muon/crab_*/results/ -name "processedLumis.json")
```

3. Getting Luminosity:
```bash
brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i DATA_samples/data_2023_Muon/combinedMuon2023processLumis.json
```