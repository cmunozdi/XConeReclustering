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
 