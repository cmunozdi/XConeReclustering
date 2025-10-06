#!/bin/bash

# Fijar el entorno
# export SCRAM_ARCH=el8_amd64_gcc11
# cd /afs/cern.ch/user/c/cmunozdi/CMSSW_15_0_3_patch1/src
# eval `scram runtime -sh`  # o simplemente: cmsenv

# Obtener la información de fastjet-contrib usando scram
FJC_INFO=$(scram tool info fastjet-contrib)
export FASTJET_CONTRIB_BASE=$(echo "$FJC_INFO" | grep "FASTJET_CONTRIB_BASE" | cut -d '=' -f 2)
export CPLUS_INCLUDE_PATH=$FASTJET_CONTRIB_BASE/include:$CPLUS_INCLUDE_PATH
export LD_LIBRARY_PATH=$FASTJET_CONTRIB_BASE/lib:$LD_LIBRARY_PATH

# Rivet
RIVET_INFO=$(scram tool info rivet)
export RIVET_BASE=$(echo "$RIVET_INFO" | grep "RIVET_BASE" | cut -d '=' -f 2)
export CPLUS_INCLUDE_PATH=$RIVET_BASE/include:$CPLUS_INCLUDE_PATH
export LD_LIBRARY_PATH=$RIVET_BASE/lib:$LD_LIBRARY_PATH
export PATH=$RIVET_BASE/bin:$PATH

# FastJet
FJ_INFO=$(scram tool info fastjet)
export FASTJET_BASE=$(echo "$FJ_INFO" | grep "FASTJET_BASE" | cut -d '=' -f 2)
export CPLUS_INCLUDE_PATH=$FASTJET_BASE/include:$CPLUS_INCLUDE_PATH
export LD_LIBRARY_PATH=$FASTJET_BASE/lib:$LD_LIBRARY_PATH
export PATH=$FASTJET_BASE/bin:$PATH

# correctionlib
corrlib_INFO=$(scram tool info correctionlib)
export CORRLIB_BASE=$(echo "$corrlib_INFO" | grep "CORRECTIONLIB_BASE" | cut -d '=' -f 2)
export CPLUS_INCLUDE_PATH=$(echo "$corrlib_INFO" | grep "INCLUDE" | cut -d '=' -f 2):$CPLUS_INCLUDE_PATH
export LD_LIBRARY_PATH=$(echo "$corrlib_INFO" | grep "LIBDIR" | cut -d '=' -f 2):$LD_LIBRARY_PATH
export PATH=$(echo "$corrlib_INFO" | grep "PATH" | cut -d '=' -f 2):$PATH

# CRAB y brilws después de cmsenv
source /cvmfs/cms.cern.ch/common/crab-setup.sh prod
source /cvmfs/cms-bril.cern.ch/cms-lumi-pog/brilws-docker/brilws-env
