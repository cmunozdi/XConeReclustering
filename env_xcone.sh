#!/bin/bash

# Obtener la información de fastjet-contrib usando scram
FJC_INFO=$(scram tool info fastjet-contrib)

# Extraer el valor de FASTJET_CONTRIB_BASE
export FASTJET_CONTRIB_BASE=$(echo "$FJC_INFO" | grep "FASTJET_CONTRIB_BASE" | cut -d '=' -f 2)

# Configurar las variables de entorno necesarias
export CPLUS_INCLUDE_PATH=$FASTJET_CONTRIB_BASE/include:$CPLUS_INCLUDE_PATH
export LD_LIBRARY_PATH=$FASTJET_CONTRIB_BASE/lib:$LD_LIBRARY_PATH

# Obtener la información de Rivet usando scram
RIVET_INFO=$(scram tool info rivet)

# Extraer el valor de RIVET_BASE
export RIVET_BASE=$(echo "$RIVET_INFO" | grep "RIVET_BASE" | cut -d '=' -f 2)

# Configurar las variables de entorno necesarias
export CPLUS_INCLUDE_PATH=$RIVET_BASE/include:$CPLUS_INCLUDE_PATH
export LD_LIBRARY_PATH=$RIVET_BASE/lib:$LD_LIBRARY_PATH
export PATH=$RIVET_BASE/bin:$PATH

# Obtener la información de fastjet usando scram
FJ_INFO=$(scram tool info fastjet)

# Extraer el valor de FASTJET_BASE
export FASTJET_BASE=$(echo "$FJ_INFO" | grep "FASTJET_BASE" | cut -d '=' -f 2)

# Configurar las variables de entorno necesarias
export CPLUS_INCLUDE_PATH=$FASTJET_BASE/include:$CPLUS_INCLUDE_PATH
export LD_LIBRARY_PATH=$FASTJET_BASE/lib:$LD_LIBRARY_PATH
export PATH=$FASTJET_BASE/bin:$PATH

source /cvmfs/cms.cern.ch/common/crab-setup.sh prod # note: this is new w.r.t. 106X instructions
source /cvmfs/cms-bril.cern.ch/cms-lumi-pog/brilws-docker/brilws-env
