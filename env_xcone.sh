#!/bin/bash

# Obtener la información de fastjet-contrib usando scram
FJC_INFO=$(scram tool info fastjet-contrib)

# Extraer el valor de FASTJET_CONTRIB_BASE
export FASTJET_CONTRIB_BASE=$(echo "$FJC_INFO" | grep "FASTJET_CONTRIB_BASE" | cut -d '=' -f 2)

# Configurar las variables de entorno necesarias
export CPLUS_INCLUDE_PATH=$FASTJET_CONTRIB_BASE/include:$CPLUS_INCLUDE_PATH
export LD_LIBRARY_PATH=$FASTJET_CONTRIB_BASE/lib:$LD_LIBRARY_PATH

source /cvmfs/cms.cern.ch/common/crab-setup.sh prod # note: this is new w.r.t. 106X instructions