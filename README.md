# XConeReclustering

## Setting up the CMSSW framework for generating PFNano files and reprocessing them with the XCone algorithm.
```bash
cmsrel CMSSW_13_03
cd CMSSW_13_0_3/src
cmsenv
git cms-addpkg PhysicsTools/NanoAOD
scram b -j 9 
```
Installing fastjet:
```bash
curl -O https://fastjet.fr/repo/fastjet-3.4.3.tar.gz
tar zxvf fastjet-3.4.3.tar.gz
cd fastjet-3.4.3/
./configure --prefix=$PWD/../fastjet-install
make
make check
make install
cd ..

wget http://fastjet.hepforge.org/contrib/downloads/fjcontrib-1.101.tar.gz
tar zxvf fjcontrib-1.101.tar.gz 
cd fjcontrib-1.101/ 
./configure --fastjet-config=$PWD/../fastjet-install/bin/fastjet-config
make
make check
make install
cd ..
```
Getting btv:
```bash
git clone https://github.com/cms-btv-pog/btvnano-prod.git
scram b -j 9 
source /cvmfs/cms.cern.ch/common/crab-setup.sh prod # note: this is new w.r.t. 106X instructions
```
Getting XCone processing:
```bash
git clone https://github.com/cmunozdi/XConeReclustering.git
```
Run locally XCone:
```bash
python3 ProcessNanoToBoostedTopQuarkWithXCone.py --input `inputfile or inputdirectory` --output `output root file` (--isMC if input is simualtion)
``` 