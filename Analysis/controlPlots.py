# . /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc11-opt/setup.sh
import ROOT
import numpy as np
import pandas as pd
import os
from glob import glob

#Bool variables to know over which sample should we run the code
is_ttbar_semi = True
is_ttbar_bck = False
is_singletop = False
is_wjets = False
is_qcdmu = False
is_dy = False
is_vv = False
is_data0 = False
is_data1 = False
OnlyGeneratePKLFiles = False
# Directory to save output files
output_dir = "/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2022postEE/DifferentBranchesMuonEraG" #"/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/DifferentBranchesMuonv1to4_JetVetoMap"
sufix_plotdir = "_G"  # Suffix for the plot directory

save_roots_boosted_selection = False


# Charging the header (adjust the path if necessary)
ROOT.gInterpreter.Declare('#include "./btagSFs/header_btag_scaleFactors.h"')
# Initializing the CorrectionSet ONLY ONCE
ROOT.gInterpreter.ProcessLine("initializeCorrectionSet();")
ROOT.gInterpreter.ProcessLine("initializeBTagCorrectionSet();")
ROOT.gInterpreter.ProcessLine("initializeMUOCorrectionSet();")
ROOT.gInterpreter.ProcessLine("initializePUReweightingCorrectionSet();")
ROOT.gInterpreter.ProcessLine("initializeJetVetoMap();")

ROOT.ROOT.EnableImplicitMT()

def get_pattern(input_dir):
    if any(keyword in input_dir for keyword in ["wjets", "ttbar"]):
        return "**/250818*/**/*.root"  # Specific pattern for wjets and ttbar
    else:
        return "**/250*/**/*.root"  # General pattern for other directories

def create_rdf_with_weights(input_dir, cross_section, luminosity, max_files=None, isTTbar=False):
    # Choose the pattern based on whether it's QCDMu or not
#     pattern = "**/250622*/**/*.root"  if "QCDMu" in input_dir else  "**/250625*/**/*.root" if "TTbar" in input_dir else "**/250625*/**/*.root"
    pattern = "**/25*/**/*.root"
    all_root_files = glob(os.path.join(input_dir, pattern), recursive=True)

    # Filter files that do not end with coor.root, lumi.root, or runs.root
    candidate_files = [
        f for f in all_root_files
        if not (f.endswith('events.root') or f.endswith('lumi.root') or f.endswith('runs.root'))
    ]

    valid_files = []
    total_gen_events_sumW = 0

    for file in candidate_files:
        f = ROOT.TFile.Open(file)
        if not f or f.IsZombie():
            f.Close()
            continue

        # Check if it has TTree "Events"
        if f.Get("Events"):
            valid_files.append(file)

        # Count generated events if there is TTree "Runs"
        runs_tree = f.Get("Runs")
        if runs_tree:
            for entry in runs_tree:
                total_gen_events_sumW += entry.genEventSumw

        f.Close()

        # If max_files is defined and we already have enough valid files, we exit the loop
        if max_files is not None and len(valid_files) >= max_files:
            break

    if not valid_files:
        print(f"Warning: No valid files with TTree 'Events' found in directory {input_dir}.")
        return None  # Or return an empty RDataFrame if needed

    if total_gen_events_sumW == 0:
        print(f"Warning: No generated events found in {input_dir}.")
        return None  # Or return an empty RDataFrame if needed

    event_weight = (luminosity * cross_section) / total_gen_events_sumW

    # Create the RDataFrame only with the valid files
    rdf = ROOT.RDataFrame("Events", valid_files)
    rdf = rdf.Define("eventWeight", f"{event_weight}*genWeight") # > 0 ? genWeight : abs(genWeight)
    return rdf


luminosity = 3.082753036*1000 #Muon era G correct JEC for data

#18.062658998*1000 #Muon0&1 v1to4 correct JEC for data 
#7.220180367*1000 #Muon0&1 v1to3 correct JEC for data 
#10.842478631*1000 #Muon0&1 v4 correct JEC for data
#1.606429085*1000 #Muon0&1 v3 correct JEC for data 
#1.270749295*1000 #Muon0&1 v2 correct JEC for data 
#4.343001987*1000 #Muon0&1 v1 correct JEC for data
#18.062658998*1000 #This new value is using /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json instead of /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json 18.084440726*1000 pb^-1 

cross_sections_singletop = [87.200000, 4.534000, 4.663095, 19.970880, 19.302840, 145.000000, 7.244000, 4.663095, 19.970880, 19.302840] #[87.200000, 1.477177, 19.302840, 145.000000, 2.360095, 19.302840] # #
cross_sections_wjets = [368.200000, 421.900000, 25.600000, 54.770000, 0.878500, 3.119000, 4427.000000, 1598.000000, 0.105300, 0.526100]
cross_sections_ttbar = [923.6*0.4392] #July 2022 value: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO#Updated_reference_cross_sections
cross_sections_ttbar_bck = [923.6*0.1061, 923.6*0.4544]
cross_sections_qcdmu = [402900, 96200, 22980, 7763, 699.1, 68.24, 21.37, 3.913, 1.323] #[1.323000, 23240.000000, 7763.000000, 699.100000, 68.240000, 404400.000000, 21.370000, 3.913000, 95900.000000]   #
cross_sections_dy = [20950, 6688] #[141500, 20950, 6688] #[20950, 141500, 6688] #
cross_sections_vv = [116.8, 54.3, 16.7] #[80.23, 29.1, 12.75] #Values from TOP-22-012 https://cms-results.web.cern.ch/cms-results/public-results/publications/TOP-22-012/index.html

if is_singletop:
    input_dirs_singletop = [
        #Without JetTightID
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/SingleTop/TbarBQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/SingleTop/TbarBtoLminusNuB-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/SingleTop/TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/SingleTop/TBbarQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/SingleTop/TBbartoLplusNuBbar-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/SingleTop/TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",

        # #With JetTightIDLepVeto inside lepton trigger function
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop/TbarBQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop/TbarBtoLminusNuB-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop/TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop/TbarWplusto4Q_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop/TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop/TBbarQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop/TBbartoLplusNuBbar-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop/TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop/TWminusto4Q_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop/TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",

        #JetTightID during matching lep-jet + LepVeto outside
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/SingleTop/TbarBQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/SingleTop/TbarBtoLminusNuB-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/SingleTop/TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/SingleTop/TbarWplusto4Q_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/SingleTop/TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/SingleTop/TBbarQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/SingleTop/TBbartoLplusNuBbar-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/SingleTop/TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/SingleTop/TWminusto4Q_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/SingleTop/TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",

        #Full production samples 2023 pre-BPix
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/SingleTop/TbarBQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/SingleTop/TbarBtoLminusNuB-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/SingleTop/TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/SingleTop/TbarWplusto4Q_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/SingleTop/TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/SingleTop/TBbarQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/SingleTop/TBbartoLplusNuBbar-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/SingleTop/TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/SingleTop/TWminusto4Q_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/SingleTop/TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",

        #Full production samples 2022 post-EE
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/SingleTop/TbarBQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/SingleTop/TbarBtoLminusNuB-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/SingleTop/TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/SingleTop/TbarWplusto4Q_TuneCP5_13p6TeV_powheg-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/SingleTop/TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/SingleTop/TBbarQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/SingleTop/TBbartoLplusNuBbar-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/SingleTop/TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/SingleTop/TWminusto4Q_TuneCP5_13p6TeV_powheg-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/SingleTop/TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",        


    ]

if is_wjets:
    input_dirs_wjets = [
        #Without JetTightID
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-200to400_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-400to600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-400to600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-40to100_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-40to100_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",

        # #With JetTightIDLepVeto inside lepton trigger function
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-200to400_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-400to600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-400to600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-40to100_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-40to100_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",

        #JetTightID during matching lep-jet + LepVeto outside
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-200to400_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-400to600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-400to600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-40to100_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-40to100_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",

        #Full production samples 2023 pre-BPix
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-200to400_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-400to600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-400to600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-40to100_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-40to100_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",

        #Full production samples 2022 post-EE
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/WJets/WtoLNu-2Jets_PTLNu-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/WJets/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/WJets/WtoLNu-2Jets_PTLNu-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/WJets/WtoLNu-2Jets_PTLNu-200to400_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/WJets/WtoLNu-2Jets_PTLNu-400to600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/WJets/WtoLNu-2Jets_PTLNu-400to600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/WJets/WtoLNu-2Jets_PTLNu-40to100_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/WJets/WtoLNu-2Jets_PTLNu-40to100_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/WJets/WtoLNu-2Jets_PTLNu-600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/WJets/WtoLNu-2Jets_PTLNu-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",        


    ]

if is_ttbar_semi:
    input_dirs_ttbar = [
        #Without JetTightID
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/TTbar/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8_moreStats",

        # #With JetTightIDLepVeto inside lepton trigger function
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/TTbar/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",

        #JetTightID during matching lep-jet + LepVeto outside
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/TTbar/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",

        #With PT reweighting
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/TTbar_TopPtReweighting/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8"

        #Full production samples 2023 pre-BPix
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/TTbar/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",

        #Full production samples 2022 post-EE
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/TTbar/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",        

    ]

if is_ttbar_bck:
    input_dirs_ttbar_bck = [
        #Without JetTightID
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/TTbar/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/TTbar/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8",

        # #With JetTightIDLepVeto inside lepton trigger function
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/TTbar/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/TTbar/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8",

        #JetTightID during matching lep-jet + LepVeto outside
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/TTbar/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/TTbar/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8",

        #With PT reweighting
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/TTbar_TopPtReweighting/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/TTbar_TopPtReweighting/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8",

        #Full production samples 2023 pre-BPix
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/TTbar/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/TTbar/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8",

        #Full production samples 2022 post-EE
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/TTbar/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/TTbar/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8",        


    ]

if is_qcdmu:
    input_dirs_qcdmu = [
        #Without JetTightID
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/QCDMu/QCD_PT-170to300_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/QCDMu/QCD_PT-300to470_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/QCDMu/QCD_PT-470to600_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/QCDMu/QCD_PT-600to800_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/QCDMu/QCD_PT-800to1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/QCDMu/QCD_PT-1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",

        # #With JetTightIDLepVeto inside lepton trigger function
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-120to170_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-15to20_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-170to300_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-20to30_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-300to470_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-30to50_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-470to600_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-50to80_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-600to800_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-800to1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-80to120_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",

        #JetTightID during matching lep-jet + LepVeto outside
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/QCDMu/QCD_PT-170to300_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/QCDMu/QCD_PT-300to470_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/QCDMu/QCD_PT-470to600_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/QCDMu/QCD_PT-600to800_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/QCDMu/QCD_PT-800to1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/QCDMu/QCD_PT-1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/QCDMu/QCD_PT-120to170_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/QCDMu/QCD_PT-50to80_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/QCDMu/QCD_PT-80to120_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",

        #Full production samples 2023 pre-BPix
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/QCDMu/QCD_PT-170to300_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/QCDMu/QCD_PT-300to470_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/QCDMu/QCD_PT-470to600_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/QCDMu/QCD_PT-600to800_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/QCDMu/QCD_PT-800to1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/QCDMu/QCD_PT-1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/QCDMu/QCD_PT-120to170_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/QCDMu/QCD_PT-80to120_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",

        #Full production samples 2022 post-EE
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/QCDMu/QCD_PT-50to80_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/QCDMu/QCD_PT-80to120_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/QCDMu/QCD_PT-120to170_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/QCDMu/QCD_PT-170to300_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/QCDMu/QCD_PT-300to470_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/QCDMu/QCD_PT-470to600_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/QCDMu/QCD_PT-600to800_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/QCDMu/QCD_PT-800to1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/QCDMu/QCD_PT-1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/QCDMu/QCD_PT-120to170_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/QCDMu/QCD_PT-80to120_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",        



    ]

if is_dy:
    input_dirs_dy = [
        #Without JetTightID
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/DY/DYto2L-2Jets_MLL-4to10_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/DY/DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/DY/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu-2Jets_MLL-105To160_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu_MLL-10to50_TuneCP5_13p6TeV_powheg-pythia8_moreStats",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu_MLL-120to200_TuneCP5_13p6TeV_powheg-pythia8_moreStats",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu_MLL-1500to2500_TuneCP5_13p6TeV_powheg-pythia8_moreStats",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu_MLL-200to400_TuneCP5_13p6TeV_powheg-pythia8_moreStats",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu_MLL-2500to4000_TuneCP5_13p6TeV_powheg-pythia8_moreStats",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu_MLL-4000to6000_TuneCP5_13p6TeV_powheg-pythia8_moreStats",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu_MLL-400to800_TuneCP5_13p6TeV_powheg-pythia8_moreStats",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu_MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8_moreStats",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu_MLL-6000_TuneCP5_13p6TeV_powheg-pythia8_moreStats",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu_MLL-800to1500_TuneCP5_13p6TeV_powheg-pythia8_moreStats",

        # #With JetTightIDLepVeto inside lepton trigger function
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/DY/DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/DY/DYto2L-2Jets_MLL-4to10_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/DY/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",

        #JetTightID during matching lep-jet + LepVeto outside
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/DY/DYto2L-2Jets_MLL-4to10_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/DY/DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/DY/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",

        #Full production samples 2023 pre-BPix
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/DY/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",

        #Full production samples 2022 post-EE
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/DY/DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/DY/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",        


    ]

if is_vv:
    input_dirs_vv = [
        #Without JetTightID
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/VV/WW_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/VV/WZ_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/VV/ZZ_TuneCP5_13p6TeV_pythia8",

        # #With JetTightIDLepVeto inside lepton trigger function
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/VV/WW_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/VV/WZ_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/MC/2023preBPix/VV/ZZ_TuneCP5_13p6TeV_pythia8",

        #JetTightID during matching lep-jet + LepVeto outside
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/VV/WW_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/VV/WZ_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto/MC/2023preBPix/VV/ZZ_TuneCP5_13p6TeV_pythia8",

        #Full production samples 2023 pre-BPix
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/VV/WW_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/VV/WZ_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/VV/ZZ_TuneCP5_13p6TeV_pythia8",

        #Full production samples 2022 post-EE
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/VV/WW_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/VV/WZ_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2022postEE/VV/ZZ_TuneCP5_13p6TeV_pythia8",        
    ]

print("Initiating the conversion of RDataFrames to DataFrames.")



ROOT.gInterpreter.Declare("""
int match_lepton_to_muon(float lep_pt, float lep_eta, float lep_phi,
                         const ROOT::RVec<float>& mu_pt,
                         const ROOT::RVec<float>& mu_eta,
                         const ROOT::RVec<float>& mu_phi) {
    /*std::cout << "Trigger lpton pt: " << lep_pt << " eta: " << lep_eta << " phi: " << lep_phi << std::endl;*/
    for (size_t i = 0; i < mu_pt.size(); ++i) {
        /*std::cout << "\tMuon index: " << i << " muon pt: " << mu_pt[i] << " eta: " << mu_eta[i] << " phi: " << mu_phi[i] << std::endl;*/
        if (std::abs(mu_pt[i] - lep_pt) < 1e-3 &&
            std::abs(mu_eta[i] - lep_eta) < 1e-3 &&
            std::abs(mu_phi[i] - lep_phi) < 1e-3) {

            return i;
        }
    }
    std::cout << "Warning: No matching muon found for lepton with pt: " << lep_pt << std::endl;
    return -1;
}
""")

ROOT.gInterpreter.Declare("""
bool pass_JetIdTight(int jet_idx, const ROOT::RVec<float>& Jet_eta,
                     const ROOT::RVec<float>& Jet_neHEF, const ROOT::RVec<float>& Jet_neEmEF,
                     const ROOT::RVec<int>& Jet_chMultiplicity, const ROOT::RVec<int>& Jet_neMultiplicity,
                     const ROOT::RVec<float>& Jet_chHEF, const ROOT::RVec<float>& Jet_muEF,
                     const ROOT::RVec<float>& Jet_chEmEF){
    bool Jet_passJetIdTight = false;
    if(abs(Jet_eta[jet_idx]) <= 2.6){
        Jet_passJetIdTight = (Jet_neHEF[jet_idx] < 0.99) && (Jet_neEmEF[jet_idx] < 0.9) && (Jet_chMultiplicity[jet_idx]+Jet_neMultiplicity[jet_idx] > 1) && (Jet_chHEF[jet_idx] > 0.01) && (Jet_chMultiplicity[jet_idx] > 0);
    }else if((abs(Jet_eta[jet_idx]) > 2.6) && (abs(Jet_eta[jet_idx]) <= 2.7)){
        Jet_passJetIdTight = (Jet_neHEF[jet_idx] < 0.90) && (Jet_neEmEF[jet_idx] < 0.99);
    }else if((abs(Jet_eta[jet_idx]) > 2.7) && (abs(Jet_eta[jet_idx]) <= 3.0)){
        Jet_passJetIdTight = (Jet_neHEF[jet_idx] < 0.90);
    }else if((abs(Jet_eta[jet_idx]) > 3.0) ){
        Jet_passJetIdTight = (Jet_neMultiplicity[jet_idx] >=2 ) && (Jet_neEmEF[jet_idx] < 0.4);
    }

    bool Jet_passJetIdTightLepVeto = false;
    if(abs(Jet_eta[jet_idx]) <= 2.7){
        Jet_passJetIdTightLepVeto = Jet_passJetIdTight && (Jet_muEF[jet_idx] < 0.8) && (Jet_chEmEF[jet_idx] < 0.8);
    }else Jet_passJetIdTightLepVeto = Jet_passJetIdTight;

    return Jet_passJetIdTightLepVeto;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::RVec<bool> pass_JetIdTightLepVeto(const ROOT::RVec<float>& Jet_eta,
                                   const ROOT::RVec<float>& Jet_neHEF, const ROOT::RVec<float>& Jet_neEmEF,
                                   const ROOT::RVec<int>& Jet_chMultiplicity, const ROOT::RVec<int>& Jet_neMultiplicity,
                                   const ROOT::RVec<float>& Jet_chHEF, const ROOT::RVec<float>& Jet_muEF,
                                   const ROOT::RVec<float>& Jet_chEmEF){
    ROOT::RVec<bool> Jet_passJetIdTight = ROOT::RVec<bool>(Jet_eta.size(), false);
    ROOT::RVec<bool> Jet_passJetIdTightLepVeto = ROOT::RVec<bool>(Jet_eta.size(), false);
    for(size_t jet_idx = 0; jet_idx < Jet_eta.size(); jet_idx++) {
        if(abs(Jet_eta[jet_idx]) <= 2.6){
            Jet_passJetIdTight[jet_idx] = (Jet_neHEF[jet_idx] < 0.99) && (Jet_neEmEF[jet_idx] < 0.9) && (Jet_chMultiplicity[jet_idx]+Jet_neMultiplicity[jet_idx] > 1) && (Jet_chHEF[jet_idx] > 0.01) && (Jet_chMultiplicity[jet_idx] > 0);
        }else if((abs(Jet_eta[jet_idx]) > 2.6) && (abs(Jet_eta[jet_idx]) <= 2.7)){
            Jet_passJetIdTight[jet_idx] = (Jet_neHEF[jet_idx] < 0.90) && (Jet_neEmEF[jet_idx] < 0.99);
        }else if((abs(Jet_eta[jet_idx]) > 2.7) && (abs(Jet_eta[jet_idx]) <= 3.0)){
            Jet_passJetIdTight[jet_idx] = (Jet_neHEF[jet_idx] < 0.90);
        }else if((abs(Jet_eta[jet_idx]) > 3.0) ){
            Jet_passJetIdTight[jet_idx] = (Jet_neMultiplicity[jet_idx] >=2 ) && (Jet_neEmEF[jet_idx] < 0.4);
        }

        if(abs(Jet_eta[jet_idx]) <= 2.7){
            Jet_passJetIdTightLepVeto[jet_idx] = Jet_passJetIdTight[jet_idx] && (Jet_muEF[jet_idx] < 0.8) && (Jet_chEmEF[jet_idx] < 0.8);
        }else Jet_passJetIdTightLepVeto[jet_idx] = Jet_passJetIdTight[jet_idx];
    }
    return Jet_passJetIdTightLepVeto;
}
""")


def apply_event_selection(rdfs, isMC=True, TWPbtag2023=0.73):
    filtered_rdfs = []

    for rdf in rdfs:
        branches = set(rdf.GetColumnNames())
        if "Num_lep_above_15" in branches:
            rdf = rdf.Redefine("Num_lep_above_15", "Sum(Muon_pt > 15 && abs(Muon_eta) < 2.4) + Sum(Electron_pt > 15 && abs(Electron_eta) < 2.4)")
        else:
            rdf = rdf.Define("Num_lep_above_15", "Sum(Muon_pt > 15 && abs(Muon_eta) < 2.4) + Sum(Electron_pt > 15 && abs(Electron_eta) < 2.4)")

        # rdf = rdf.Redefine('muon', 'triggerLepton(Muon_pt, Muon_eta, Muon_phi, Muon_pdgId, Muon_jetIdx, Muon_tightId, PFCands_pt, PFCands_eta, PFCands_phi, \
        #      PFCands_pdgId, Jet_pt, Jet_eta, Jet_phi, Jet_mass, Jet_passJetIdTight, Jet_rawFactor, Jet_area, Rho_fixedGridRhoFastjetAll, run, 25, true, false)') \
        #          .Redefine('Electron_tightId', 'Electron_cutBased >= 4') \
        #          .Redefine('electron', 'triggerLepton(Electron_pt, Electron_eta, Electron_phi, Electron_pdgId, Electron_jetIdx, Electron_tightId, \
        #             PFCands_pt, PFCands_eta, PFCands_phi, PFCands_pdgId, Jet_pt, Jet_eta, Jet_phi, Jet_mass, Jet_passJetIdTight, Jet_rawFactor, Jet_area, Rho_fixedGridRhoFastjetAll, run, 25, false, false)') \
        #          .Redefine('lepton', 'CombineLeptons(muon, electron)') \
        #          .Redefine('n_leptons', 'lepton.n_lep') \
                #  .Define('selected_lepton1_idx', 'ROOT::VecOps::ArgMax(lepton1.pt > 55 && abs(lepton1.eta) < 2.4)') \
                #  .Define('is_muon1_trg', 'abs(lepton1.pdgId[selected_lepton1_idx]) == 13') \
                #  .Define('printLepInfo', 'printValuesOfLepton(lepton)')


        rdf_filtered = (
            rdf
            # 1. Detector selection
            .Filter("pass_detector_selection", "Detector selection") #n_leptons == 1 && n_jets > 0 && pt_miss > 37.5 && n_fatjets > 0
            .Define("Jet_vetoMap", "is_event_not_vetoed_by_jetVetoMap(Jet_eta, Jet_phi)")

#             .Filter("lepton.n_lep==1")

            # 2. Lepton selection
            # .Redefine('muon', 'triggerLepton(Muon_pt, Muon_eta, Muon_phi, Muon_pdgId, Muon_jetIdx, Muon_tightId, PFCands_pt, PFCands_eta, PFCands_phi, \
            #      PFCands_pdgId, Jet_pt, Jet_eta, Jet_phi, Jet_mass, Jet_passJetIdTight, Jet_rawFactor, Jet_area, Rho_fixedGridRhoFastjetAll,, 25, true, false)') \
            # .Redefine('Electron_tightId', 'Electron_cutBased >= 4') \
            # .Redefine('electron', 'triggerLepton(Electron_pt, Electron_eta, Electron_phi, Electron_pdgId, Electron_jetIdx, Electron_tightId, \
            #      PFCands_pt, PFCands_eta, PFCands_phi, PFCands_pdgId, Jet_pt, Jet_eta, Jet_phi, Jet_mass, Jet_passJetIdTight, Jet_rawFactor, Jet_area, Rho_fixedGridRhoFastjetAll, run, 25, false, false)') \
            # .Redefine('lepton', 'CombineLeptons(muon, electron)') \
            # .Redefine('n_leptons', 'lepton.n_lep') \
#             .Define('lepton_trg', 'triggerLepton(Muon_pt, Muon_eta, Muon_phi, Muon_pdgId, Muon_jetIdx, Muon_tightId, PFCands_pt, \
#                      PFCands_eta, PFCands_phi, PFCands_pdgId, Jet_pt, Jet_eta, Jet_phi, true, true)')
#             .Define('muon_trg', 'triggerLepton(Muon_pt, Muon_eta, Muon_phi, Muon_pdgId, Muon_jetIdx, Muon_tightId, PFCands_pt, \
#                       PFCands_eta, PFCands_phi, PFCands_pdgId, Jet_pt, Jet_eta, Jet_phi, Jet_mass, Jet_rawFactor, Jet_area, \
#                       Rho_fixedGridRhoFastjetAll, 25, true, true)')
#             .Define('electron_trg', 'triggerLepton(Electron_pt, Electron_eta, Electron_phi, Electron_pdgId, Electron_jetIdx, \
#                     Electron_tightId, PFCands_pt, PFCands_eta, PFCands_phi, PFCands_pdgId, Jet_pt, Jet_eta, Jet_phi, Jet_mass, \
#                     Jet_rawFactor, Jet_area, Rho_fixedGridRhoFastjetAll, 25, false, true)')
#             .Define('lepton_trg', 'CombineLeptons(muon_trg, electron_trg)')
# # #             This is for when lepton branch have already included dR and pt_rel
# # #             .Define('lepton_trg', 'triggerLepton(lepton.pt, lepton.eta, lepton.phi, lepton.pdgId, \
# # #                      PFCands_pt, PFCands_eta, PFCands_phi, PFCands_pdgId, Jet_eta, Jet_phi, true, true)')
#             .Define('n_leptons_trg', 'lepton_trg.n_lep')
#             .Filter('n_leptons_trg == 1')
#             .Define('lepton_trg_pt', 'lepton_trg.pt')
#             .Define('lepton_trg_eta', 'lepton_trg.eta')
#             .Define('lepton_trg_phi', 'lepton_trg.phi')
#             .Define('lepton_trg_pdgId', 'lepton_trg.pdgId')
#             .Define('lepton_trg_dR', 'lepton_trg.dR_to_jet')
#             .Define('lepton_trg_ptrel', 'lepton_trg.pt_rel_to_jet')
#             .Define('lepton_trg_n_lep', 'lepton_trg.n_lep')
            # .Filter('Sum(lepton.pt>55 && abs(lepton.eta)<2.4) == 1')
# #             .Filter("lepton.n_lep>=1")
#             .Define('lepton_trg_pt', 'lepton.pt[0]')
#             .Define('lepton_trg_eta', 'lepton.eta[0]')
#             .Define('lepton_trg_phi', 'lepton.phi[0]')
#             .Define('lepton_trg_pdgId', 'lepton.pdgId[0]')
#             .Define('lepton_trg_dR', 'lepton.dR_to_jet[0]')
#             .Define('lepton_trg_ptrel', 'lepton.pt_rel_to_jet[0]')
#             .Define('lepton_trg_n_lep', 'lepton.n_lep')
            .Define('selected_lepton_idx', 'ROOT::VecOps::ArgMax(lepton.pt > 55 && abs(lepton.eta) < 2.4)') \
            .Define('is_muon_trg', 'abs(lepton.pdgId[selected_lepton_idx]) == 13') \
            # .Filter('is_muon_trg') \
            # .Filter('(lepton.dR_to_jet[selected_lepton_idx] > 3.)&&(lepton.dR_to_jet[selected_lepton_idx] < 3.2)') \
            .Define('lepton_trg_pt', 'lepton.pt[selected_lepton_idx]') \
            .Define('lepton_trg_eta', 'lepton.eta[selected_lepton_idx]') \
            .Define('lepton_trg_phi', 'lepton.phi[selected_lepton_idx]') \
            .Define('lepton_trg_pdgId', 'lepton.pdgId[selected_lepton_idx]') \
            .Define('lepton_trg_dR', 'lepton.dR_to_jet[selected_lepton_idx]') \
            .Define('lepton_trg_ptrel', 'lepton.pt_rel_to_jet[selected_lepton_idx]') \
            .Define('lepton_trg_jet_pt', 'lepton.closest_jet_pt[selected_lepton_idx]') \
            .Define('num_muon_selected', 'Sum((Muon_pt == lepton_trg_pt) && (Muon_eta == lepton_trg_eta) && (Muon_phi == lepton_trg_phi) && (Muon_highPtId == 2) && (Muon_tkIsoId >= 1))') \
            # .Filter('num_muon_selected == 1', "Just 1 selected muon") \
            .Define(
               "lepton_tuneP_pt",
               "Muon_tunepRelPt[ROOT::VecOps::ArgMax((Muon_pt == lepton_trg_pt) && (Muon_eta == lepton_trg_eta) && (Muon_phi == lepton_trg_phi) && (Muon_highPtId == 2) && (Muon_tkIsoId >= 1))] * lepton_trg_pt"
            ) \
            # .Redefine('lepton_trg_pt', 'lepton_tuneP_pt') \
            .Define('lepton_trg_n_lep', 'Sum((lepton.pt > 55) && (abs(lepton.eta) < 2.4))') \
            .Define('lepton_jet_corrFactor', 'lepton.closest_jet_pt[selected_lepton_idx] / Jet_pt[lepton.closest_jet_idx[selected_lepton_idx]]') \
            #Si ya existe la rama Num_lep_above_15, redefinela, si no, créala
            # .Redefine('Num_lep_above_15', "Sum(Muon_pt > 15 && abs(Muon_eta) < 2.4) + Sum(Electron_pt > 15 && abs(Electron_eta) < 2.4)")
            # .Filter('Num_lep_above_15 == 1', "Only one lepton above 15 GeV") \
            # .Define('lepton_trg_jet_idx2', 'lepton.closest_jet_idx[selected_lepton_idx]') \
            # .Define('Jet_passJetIdTightLepVeto2', 'pass_JetIdTight(lepton_trg_jet_idx2, Jet_eta, Jet_neHEF, Jet_neEmEF, Jet_chMultiplicity, Jet_neMultiplicity, Jet_chHEF, Jet_muEF, Jet_chEmEF)') \
            # .Filter('Jet_passJetIdTightLepVeto2', 'Jet passes tight ID with lepton veto') \
            # .Define('lepton_muon_idx',
            #    'match_lepton_to_muon(lepton_trg_pt, lepton_trg_eta, lepton_trg_phi, Muon_pt, Muon_eta, Muon_phi)') \
            # .Filter('lepton_muon_idx > -1 && Muon_tkIsoId[lepton_muon_idx] >= 1',
            #    'Muon_tkIsoId == 1 if lepton matches a muon')
            # .Filter('lepton_trg_n_lep == 1')

            # 2. Al menos un jet con las condiciones
            # .Filter(
            #     f"Sum(Jet_pt > 30 && abs(Jet_eta) < 2.4 && Jet_btagDeepFlavB > {TWPbtag2023}) >= 1",
            #     "Jet selection"
            # )  Jet_passJetIdTightLepVeto
            # .Define("Jet_passAllJetIdTightLepVeto", "pass_JetIdTightLepVeto(Jet_eta, Jet_neHEF, Jet_neEmEF, Jet_chMultiplicity, Jet_neMultiplicity, Jet_chHEF, Jet_muEF, Jet_chEmEF)") \
            # .Redefine("Jet_passAllJetIdTightLepVeto", "std::vector<bool>(Jet_passAllJetIdTightLepVeto.begin(), Jet_passAllJetIdTightLepVeto.end())") \
            .Define(
                "n_bjets",
                f"Sum((Jet_pt > 30) &&( abs(Jet_eta) < 2.4) && (Jet_btagDeepFlavB > {TWPbtag2023}) && Jet_passJetIdTight)"
            )
            .Define(
                "bjet_pt",
                f"Jet_pt[(Jet_pt > 30) &&( abs(Jet_eta) < 2.4) && (Jet_btagDeepFlavB > {TWPbtag2023}) && Jet_passJetIdTight]"
            )
            .Define(
                "bjet_eta",
                f"Jet_eta[(Jet_pt > 30) &&( abs(Jet_eta) < 2.4) && (Jet_btagDeepFlavB > {TWPbtag2023}) && Jet_passJetIdTight]"
            )
            .Define(
                "bjet_btag",
                f"Jet_btagDeepFlavB[(Jet_pt > 30) &&( abs(Jet_eta) < 2.4) && (Jet_btagDeepFlavB > {TWPbtag2023}) && Jet_passJetIdTight]"
            )
            # .Filter("n_bjets >= 1", "At least one b-tagged jet with pt > 30 GeV and |eta| < 2.4 passing tight ID with lepton veto")
            # .Filter(
            #     f"Sum(Jet_pt > 15)>=4",
            #     "4 jets at least"
            # )

            .Define("n_jets", "Sum(Jet_pt > 30 && abs(Jet_eta) < 2.4 && Jet_passJetIdTight)")
            # 3. pt_miss > 50
#             .Filter("pt_miss > 50", "Missing pt selection")
            # .Filter("PuppiMET_pt > 50", "Missing pt selection")

            # 4. Create the XCone clustering
#             .Define("jet_XConeFromPFCands", "buildXConeJets(PFCands_pt, PFCands_eta, PFCands_phi, PFCands_mass, PFCands_pdgId)")\
#             .Define("n_fatjets", "jet_XConeFromPFCands.fatjets.n_jets")
# #             .Filter("n_fatjets>0")
# # #             .Filter("Sum(jet_XConeFromPFCands.fatjets.pt>200)>0")
# #             # 4. Fatjet sin el lepton con pt > 400
#             .Filter(
#                 f"Sum(fatjets.pt[fatjets.findLepton == 0] > 400) >=1",
#                 "Fatjet pt without lepton"
#             )


#             # 5. Subjets dentro del topjet con pt>30 y eta<2.5
            .Define("num_subjets_in_topjet_pt30_eta2p5", f"Sum(subjets.pt[topjets.subjets_in_topjet] > 30 && abs(subjets.eta[topjets.subjets_in_topjet]) < 2.5)")
            # .Filter("num_subjets_in_topjet_pt30_eta2p5==3", "Subjet pt and eta inside topjet")

            # 6. There must be 3 subjets inside the topjet
            # .Filter("topjets.n_subjets == 3", "Three subjets inside the topjet")
            # 7. The XCone jet from hadronic side must have pt > 350 GeV
            # .Filter("topjets.pt > 350")
            # .Filter("fatjets.n_jets > 0", "At least one fatjet from XCone")

#             .Filter("HLT_Mu50", "triggerMuons")
            .Define("pass_trigger",
                    f"(abs(lepton.pdgId[selected_lepton_idx]) == 13 && (HLT_Mu50))") # || \
#                      (abs(lepton.pdgId[0]) == 11 && lepton.pt[0] > 120 && (HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Photon175)) || \
#                      (abs(lepton.pdgId[0]) == 11 && lepton.pt[0] < 120 && HLT_Ele30_WPTight_Gsf)")
            # .Filter("pass_trigger", "Trigger selection")

            # MET filters
            # .Define("pass_MET_filters_2023", "Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_BadPFMuonDzFilter && Flag_hfNoisyHitsFilter && Flag_eeBadScFilter")
            # .Filter("pass_MET_filters_2023", "MET filters selection")

        )
        if isMC:
            
            # MET filters
            # rdf_filtered = rdf_filtered.Define("pass_MET_filters_2023", "Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_BadPFMuonDzFilter && Flag_hfNoisyHitsFilter && Flag_eeBadScFilter")
            # Apply b-tagging correction if it is MC
            branches = set(rdf_filtered.GetColumnNames())
            if "btagWeight" in branches:
                rdf_filtered = rdf_filtered.Redefine("btagWeight", "compute_btagWeight(Jet_pt, Jet_eta, Jet_hadronFlavour, Jet_btagDeepFlavB)") \
                                        #    .Redefine("eventWeight", "eventWeight * btagWeight")
            else:
                rdf_filtered = rdf_filtered.Define("btagWeight", "compute_btagWeight(Jet_pt, Jet_eta, Jet_hadronFlavour, Jet_btagDeepFlavB)") \
                                        #    .Redefine("eventWeight", "eventWeight * btagWeight")

            if "TopPtWeight_NNLOpNLOEW" not in branches:
            #     rdf_filtered = rdf_filtered.Redefine("eventWeight", "eventWeight * TopPtWeight_NNLOpNLOEW")
            # else:
                rdf_filtered = rdf_filtered.Define("TopPtWeight_NNLOpNLOEW", "1.")
                

            rdf_filtered = rdf_filtered.Define("muoWeight", "compute_MUOWeight(lepton_trg_pt, lepton_trg_eta)") \
                                    #    .Redefine("eventWeight", "eventWeight * muoWeight")

            rdf_filtered = rdf_filtered.Define("PUReweighting", "compute_PUReweight(Pileup_nTrueInt)") \
                                    #    .Redefine("eventWeight", "eventWeight * PUReweighting")

        # else:
        #     #ONLY FOR DATA SINCE SEEMS THAT IT IS NOT WORKING DURING NTUPLE PRODUCTION
        #     rdf_filtered = rdf_filtered.Filter("topjets.n_subjets ==3 && topjets.pt > 200", "At least one topjet with 3 subjets and pt > 200 GeV")

        filtered_rdfs.append(rdf_filtered)

    return filtered_rdfs


columns_mc = ["eventWeight", "lepton_trg_pt", "lepton_tuneP_pt", "lepton_trg_eta", "lepton_trg_jet_pt", "pass_trigger", "PUReweighting", "pass_MET_filters_2022", "Num_lep_above_15", "num_muon_selected", "is_muon_trg", "PuppiMET_pt", "n_bjets", "bjet_pt", "bjet_eta", "n_jets", "TopPtWeight_NNLOpNLOEW", "muoWeight", "btagWeight", "topjets.pt", "topjets.eta", "topjets.mass", "num_subjets_in_topjet_pt30_eta2p5", "topjets.n_subjets", "subjets.pt", "subjets.eta", "Jet_pt", "Jet_eta", "Jet_phi", "Jet_vetoMap", "Jet_hadronFlavour", "Jet_btagDeepFlavB", "Pileup_nTrueInt", "HLT_HighPtTkMu100", "HLT_CascadeMu100", "HLT_IsoMu24", "HLT_Mu50"]

# columns_mc = ["eventWeight", "lepton", "lepton_trg_pt", "lepton_tuneP_pt", "lepton_trg_eta", "lepton_trg_jet_pt", "pass_trigger", "PUReweighting", "pass_MET_filters_2023", "Num_lep_above_15", "num_muon_selected", "is_muon_trg", "PuppiMET_pt", "n_bjets", "TopPtWeight_NNLOpNLOEW", "muoWeight", "btagWeight", "topjets.pt", "num_subjets_in_topjet_pt30_eta2p5", "topjets.n_subjets", "Jet_pt", "Jet_eta", "Jet_hadronFlavour", "Jet_btagDeepFlavB", "Pileup_nTrueInt"]

columns_data = [col for col in columns_mc if col not in ["eventWeight", "btagWeight", "muoWeight", "PUReweighting", "Jet_hadronFlavour", "TopPtWeight_NNLOpNLOEW", "Pileup_nTrueInt"]]
# columns_data = columns_data + ["Jet_passJetIdTightLepVeto"]

# DATA: MUON 0
if is_data0:
    print("Creating RDataFrames for Data 0 samples...")
    input_dirs_data0 = "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/Data/2022postEE/Muon/Run2022G-22Sep2023-v1_BTV_Run3_2022_Comm_MINIAODv4" #"/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/Data/2023preBPix/Muon0_CorrectedJEC/Run2023C-22Sep2023_v3*" #-v1_BTV_Run3_2023_Comm_MINIAODv4/" #"/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/Data/2023preBPix/Muon0/**"  #"/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/Data/2023preBPix/Muon0/**" #
    root_files_data0 = glob(os.path.join(input_dirs_data0, "**/**/*.root"))#"250818*/**/*.root")) #"**/**/*.root"))
    rdf_data0 = ROOT.RDataFrame("Events", root_files_data0)

    print("Applying event selection to Data 0...")
    rdf_data0 = apply_event_selection([rdf_data0], isMC=False)[0]
    if save_roots_boosted_selection:
        columns = list(rdf_data0.GetColumnNames())
        rdf_clean = rdf_data0
        for col in columns:
            try:
                rdf_clean = rdf_clean.FilterAvailable(col)
                print(f"Column {col} filtered.")
            except:
                print(f"Column {col} could not be filtered.")
        rdf_clean.Snapshot('Events', "/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/OnlyPreProductionSelectionMuonV2/rootFilesBoostedSelection/data0_boosted_selection.root", columns)


    print("Initiating the conversion of RDataFrames to DataFrames for Data 0...")
    df_data0 = pd.DataFrame(rdf_data0.AsNumpy(columns=columns_data))
    print("Conversion completed for Data 0.")

    print("Saving DataFrame of Data 0 in pickle format...")
    # Using general variable output path    
    df_data0.to_pickle(f"{output_dir}/dataE.pkl")
    print("File saved: dataE.pkl")
    print("Number of events: ", len(df_data0))
    del rdf_data0  # Release memory
    if OnlyGeneratePKLFiles:
        del df_data0
else:
    if(not OnlyGeneratePKLFiles):
        print("Loading DataFrame of Data 0 from existing pickle file...")
        # Load several pkl and concatenate them
        pkl_files_data0 = [
                            # f"{output_dir}/data0_v1.pkl",
                            # f"{output_dir}/data0_v2.pkl",
                            # f"{output_dir}/data0_v3.pkl",
                            # f"{output_dir}/data0_v4.pkl",
                            f"{output_dir}/dataE.pkl"
                          ]

        dfs = [pd.read_pickle(pkl) for pkl in pkl_files_data0]
        df_data0 = pd.concat(dfs, ignore_index=True)

        del dfs  # Release memory from the list of DataFrames

        # df_data0.to_pickle(f"{output_dir}/data0_v1to3.pkl")
        print("Completion of loading for Data 0.")

# # DATA: MUON 1
# if is_data1:
#     print("Creating RDataFrames for Data 1 samples...")
#     input_dirs_data1 = "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/Data/2023preBPix/Muon1_CorrectedJEC/Run2023C-22Sep2023_v3*" #4-v2_BTV_Run3_2023_Comm_MINIAODv4/" #"/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples/Data/2023preBPix/Muon1/**" #"/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightID/Data/2023preBPix/Muon1/**" #
#     root_files_data1 = glob(os.path.join(input_dirs_data1, "**/**/*.root"))#"250818*/**/*.root")) #"**/**/*.root"))
#     rdf_data1 = ROOT.RDataFrame("Events", root_files_data1)

#     print("Applying event selection to Data 1...")
#     rdf_data1 = apply_event_selection([rdf_data1], isMC=False)[0]
#     if save_roots_boosted_selection:
#         columns = list(rdf_data0.GetColumnNames())
#         rdf_data1.Snapshot("Events", "/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/OnlyPreProductionSelectionMuonV2/rootFilesBoostedSelection/data1_boosted_selection.root", columns)

#     print("Initiating the conversion of RDataFrames to DataFrames for Data 1...")
#     df_data1 = pd.DataFrame(rdf_data1.AsNumpy(columns=columns_data))
#     print("Conversion completed for Data 1.")

#     print("Saving DataFrame of Data 1 in pickle format...")
#     df_data1.to_pickle(f"{output_dir}/data1_v3.pkl")
#     print("File saved: data1_v3.pkl")
#     print("Number of events: ", len(df_data1))
#     del rdf_data1  # Release memory
#     if OnlyGeneratePKLFiles:
#         del df_data1
# else:
#     if(not OnlyGeneratePKLFiles):
#         print("Loading DataFrame of Data 1 from existing pickle file...")
#         # Load several pkl and concatenate them
#         pkl_files_data1 = [
#                             f"{output_dir}/data1_v1.pkl",
#                             f"{output_dir}/data1_v2.pkl",
#                             f"{output_dir}/data1_v3.pkl",
#                             f"{output_dir}/data1_v4.pkl"
#                           ]
#         dfs = [pd.read_pickle(pkl) for pkl in pkl_files_data1]
#         df_data1 = pd.concat(dfs, ignore_index=True)
#         del dfs  # Release memory from the list of DataFrames

#         # df_data1.to_pickle(f"{output_dir}/data1_v1to3.pkl")
#         print("Completion of loading for Data 1.")
if not OnlyGeneratePKLFiles:
    df_data = df_data0 #pd.concat([df_data0, df_data1], ignore_index=True)
    del df_data0
    # del df_data1

factor = 1. #luminosity / (18.062658998 * 1000)


# MC: SINGLE TOP
if is_singletop:
    print("Creating RDataFrames for Single Top samples...")
    rdfs_singletop = [
        create_rdf_with_weights(path, xs, luminosity)#, max_files=5)
        for path, xs in zip(input_dirs_singletop, cross_sections_singletop)
    ]

    print("Applying event selection to Single Top samples...")
    rdfs_singletop = apply_event_selection(rdfs_singletop)
    if save_roots_boosted_selection:
        columns = list(rdf_data0.GetColumnNames())
        rdfs_singletop.Snapshot("Events", "/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/OnlyPreProductionSelection/rootFilesBoostedSelection/singletop_boosted_selection.root", columns)

    print("Initiating the conversion of RDataFrames to DataFrames for SingleTop...")
    df_singletop = pd.concat([
        pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        for rdf in rdfs_singletop
    ], ignore_index=True)
    print("Conversion completed for SingleTop.")
    print("Saving DataFrame of SingleTop in pickle format...")
    df_singletop.to_pickle(f"{output_dir}/singletop.pkl")
    print("File saved: singletop.pkl")
    print("Number of events: ", len(df_singletop))
    del rdfs_singletop  # Release memory
    
    if OnlyGeneratePKLFiles:
        del df_singletop
    else: df_singletop['eventWeight'] = df_singletop['eventWeight'] * factor
else:
    if(not OnlyGeneratePKLFiles):
        print("Loading DataFrame of Single Top from existing pickle file...")
        df_singletop = pd.read_pickle(f"{output_dir}/singletop.pkl")
        df_singletop['eventWeight'] = df_singletop['eventWeight'] * factor
        # df_singletop.to_pickle(f"{output_dir}/singletop.pkl")
        print("Completion of loading for Single Top.")

# MC: W+JETS
if is_wjets:
    for i, (path, xs) in enumerate(zip(input_dirs_wjets, cross_sections_wjets)):
        print(f"Processing sample {i+1}/{len(input_dirs_wjets)}: {path}")
        rdf = create_rdf_with_weights(path, xs, luminosity)#, max_files=5)
        rdf = [rdf]  # put in list to use apply_event_selection
        rdf = apply_event_selection(rdf)
        rdf = rdf[0]
        df_part = pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        pickle_path = f"{output_dir}/wjets_part{i}.pkl"
        df_part.to_pickle(pickle_path)
        print("File saved:", pickle_path)

        del df_part, rdf
        ROOT.gSystem.ProcessEvents()  # Clean ROOT buffers

    if not OnlyGeneratePKLFiles:
        # Define the pattern to find all partial pickle files
        pkl_files = sorted(glob(f"{output_dir}/wjets_part*.pkl"))

        # Read all pickle files into a list of DataFrames
        dfs = [pd.read_pickle(pkl) for pkl in pkl_files]

        for i, df in enumerate(dfs):
            df['eventWeight'] = df['eventWeight'] * factor

        # Concatenate all DataFrames into one
        df_wjets = pd.concat(dfs, ignore_index=True)

        del dfs  # Release memory
    # print("Creating RDataFrames for W+jets samples...")
    # rdfs_wjets = [
    #     create_rdf_with_weights(path, xs, luminosity)#, max_files=100)
    #     for path, xs in zip(input_dirs_wjets, cross_sections_wjets)
    # ]

    # print("Applying event selection to W+jets samples...")
    # rdfs_wjets = apply_event_selection(rdfs_wjets)
    # if save_roots_boosted_selection:
    #     columns = list(rdf_data0.GetColumnNames())
    #     rdfs_wjets.Snapshot("Events", "/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/OnlyPreProductionSelection/rootFilesBoostedSelection/wjets_boosted_selection.root", columns)

    # print("Iniciando la conversión de RDataFrames a DataFrames para WJets...")
    # df_wjets = pd.concat([
    #     pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
    #     for rdf in rdfs_wjets
    # ], ignore_index=True)
    # print("Conversión completada para WJets.")
    # # df_wjets = df_wjets[
    # #     (abs(df_wjets['lepton_trg_pdgId']) == 13) &
    # # #     (
    # # #         df_wjets['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
    # # #         df_wjets['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
    # # #         == 1
    # # #     )
    # #     (
    # #         df_wjets.apply( lambda row: sum(pt > 30 and abs(eta) < 2.4 and btag > 0.6553 and passjetid
    # #                                             for pt, eta, btag, passjetid in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'], row['Jet_passJetIdTightLepVeto'])) >= 1, axis=1)
    # #     )
    # # ]
    # print("Guardando DataFrame de WJets en formato pickle...")
    # df_wjets.to_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/OnlyPreProductionSelection/wjets.pkl")
    # print("Archivo guardado: wjets.pkl")
    # print("Number of events: ", len(df_wjets))
    # del rdfs_wjets  # Liberar memoria
else:
    if(not OnlyGeneratePKLFiles):
        # Define the pattern to find all partial pickle files
        pkl_files = sorted(glob(f"{output_dir}/wjets_part*.pkl"))

        # Read all pickle files into a list of DataFrames
        print("Loading DataFrame of W+jets from existing pickle files...")
        dfs = [pd.read_pickle(pkl) for pkl in pkl_files]
        for i, df in enumerate(dfs):
            df['eventWeight'] = df['eventWeight'] * factor
            # df.to_pickle(f"{output_dir}/wjets_part{i}.pkl")

        # Concatenate all DataFrames into one
        df_wjets = pd.concat(dfs, ignore_index=True)
        print("Loading completed for W+jets.")
        del dfs  # Release memory

    # print("Cargando DataFrame de W+jets desde archivo pickle existente...")
    # df_wjets = pd.read_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/OnlyPreProductionSelection/wjets.pkl")
    # print("Carga completada para W+jets.")

# MC: TTBAR SEMILEPTONIC
if is_ttbar_semi:
    print("Creating RDataFrames for ttbar semileptonic samples...")
    rdfs_ttbar = [
        create_rdf_with_weights(path, xs, luminosity, isTTbar=True)#, max_files=5)
        for path, xs in zip(input_dirs_ttbar, cross_sections_ttbar)
    ]

    print("Applying event selection to ttbar semileptonic samples...")
    rdfs_ttbar= apply_event_selection(rdfs_ttbar)
    if save_roots_boosted_selection:
        columns = list(rdf_data0.GetColumnNames())
        rdfs_ttbar.Snapshot("Events", "/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/OnlyPreProductionSelection/rootFilesBoostedSelection/ttbar_semi_boosted_selection.root", columns)

    print("Initiating the conversion of RDataFrames to DataFrames for TTBar...")
    df_ttbar = pd.concat([
        pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        for rdf in rdfs_ttbar
    ], ignore_index=True)
    print("Conversion completed for TTBar.")
    print("Saving DataFrame of TTBar in pickle format...")
    df_ttbar.to_pickle(f"{output_dir}/ttbar_semi.pkl")
    print("File saved: ttbar_semi.pkl")
    print("Number of events: ", len(df_ttbar))
    del rdfs_ttbar  # Release memory
    if OnlyGeneratePKLFiles:
        del df_ttbar
    else: df_ttbar['eventWeight'] = df_ttbar['eventWeight'] * factor
else:
    if (not OnlyGeneratePKLFiles):
        print("Loading DataFrame of TTBar semileptonic from existing pickle file...")
        df_ttbar = pd.read_pickle(f"{output_dir}/ttbar_semi.pkl")
        df_ttbar['eventWeight'] = df_ttbar['eventWeight'] * factor
        # df_ttbar.to_pickle(f"{output_dir}/ttbar_semi.pkl")
        print("Loading completed for TTBar semileptonic.")

    # # Suponiendo que el peso total se llama "weight" (ajusta si es distinto)
    # df_ttbar["eventWeight"] = df_ttbar["eventWeight"] / 0.73
    # print("Guardando ttbar signal sin extra weight")
    # # Guarda un nuevo pickle sin el factor
    # df_ttbar.to_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/UsingLoosIsoAndTunepID_with0p7ttbarWeight_topPtRew_btagCorr_MuoCorr/ttbar_semi_beforeBTagMoreStats.pkl")

# MC: TTBAR BACKGROUND
if is_ttbar_bck:
    print("Creating RDataFrames for ttbar background samples...")
    rdfs_ttbar_bck = [
        create_rdf_with_weights(path, xs, luminosity, isTTbar=True)#, max_files=5)
        for path, xs in zip(input_dirs_ttbar_bck, cross_sections_ttbar_bck)
    ]

    print("Applying event selection to ttbar background samples...")
    rdfs_ttbar_bck= apply_event_selection(rdfs_ttbar_bck)
    if save_roots_boosted_selection:
        columns = list(rdf_data0.GetColumnNames())
        rdfs_ttbar_bck.Snapshot("Events", "/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/OnlyPreProductionSelection/rootFilesBoostedSelection/ttbar_bck_boosted_selection.root", columns)

    print("Initiating the conversion of RDataFrames to DataFrames for TTBar Background...")
    df_ttbar_bck = pd.concat([
        pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        for rdf in rdfs_ttbar_bck
    ], ignore_index=True)
    print("Conversion completed for TTBar Background.")
    print("Saving DataFrame of TTBar Background in pickle format...")
    df_ttbar_bck.to_pickle(f"{output_dir}/ttbar_bck.pkl")
    print("File saved: ttbar_bck.pkl")
    print("Number of events: ", len(df_ttbar_bck))
    del rdfs_ttbar_bck  # Release memory
    if not OnlyGeneratePKLFiles:
        df_ttbar_bck['eventWeight'] = df_ttbar_bck['eventWeight'] * factor
    else:
        del df_ttbar_bck
else:
    if (not OnlyGeneratePKLFiles):
        print("Loading DataFrame of TTBar Background from existing pickle file...")
        df_ttbar_bck = pd.read_pickle(f"{output_dir}/ttbar_bck.pkl")
        df_ttbar_bck['eventWeight'] = df_ttbar_bck['eventWeight'] * factor
        # df_ttbar_bck.to_pickle(f"{output_dir}/ttbar_bck.pkl")
        print("Loading completed for TTBar Background.")

    # df_ttbar_bck["eventWeight"] = df_ttbar_bck["eventWeight"] / 0.73
    # print("Guardando ttbar bck sin extra weight")
    # #Guarda un nuevo pickle sin el factor
    # df_ttbar_bck.to_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/UsingLoosIsoAndTunepID_with0p7ttbarWeight_topPtRew_btagCorr_MuoCorr/ttbar_bck_beforeBTag.pkl")

# MC: QCD MUON
if is_qcdmu:
    print("Creating RDataFrames for QCD muon samples...")
    rdfs_qcdmu = [
        create_rdf_with_weights(path, xs, luminosity)#, max_files=5)
        for path, xs in zip(input_dirs_qcdmu, cross_sections_qcdmu)
    ]

    print("Applying event selection to QCD muon samples...")
    rdfs_qcdmu = apply_event_selection(rdfs_qcdmu)
    if save_roots_boosted_selection:
        columns = list(rdf_data0.GetColumnNames())
        rdfs_qcdmu.Snapshot("Events", "/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/OnlyPreProductionSelection/rootFilesBoostedSelection/qcdmu_boosted_selection.root", columns)

    print("Initiating the conversion of RDataFrames to DataFrames for QCDMu...")
    df_qcdmu = pd.concat([
        pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        for rdf in rdfs_qcdmu
    ], ignore_index=True)
    print("Conversion completed for QCDMu.")
    print("Saving DataFrame of QCDMu in pickle format...")
    df_qcdmu.to_pickle(f"{output_dir}/qcdmu.pkl")
    print("File saved: qcdmu.pkl")
    print("Number of events: ", len(df_qcdmu))
    del rdfs_qcdmu  # Release memory
    if not OnlyGeneratePKLFiles:
        df_qcdmu['eventWeight'] = df_qcdmu['eventWeight'] * factor
    else:
        del df_qcdmu
else:
    if (not OnlyGeneratePKLFiles):
        print("Loading DataFrame of QCDMu from existing pickle file...")
        df_qcdmu = pd.read_pickle(f"{output_dir}/qcdmu.pkl")
        df_qcdmu['eventWeight'] = df_qcdmu['eventWeight'] * factor
        # df_qcdmu.to_pickle(f"{output_dir}/qcdmu.pkl")
        print("Loading completed for QCDMu.")

# MC: DY
if is_dy:
    print("Creating RDataFrames for DY samples...")
    rdfs_dy = [
        create_rdf_with_weights(path, xs, luminosity)#, max_files=5)
        for path, xs in zip(input_dirs_dy, cross_sections_dy)
    ]

    print("Applying event selection to DY samples...")
    rdfs_dy = apply_event_selection(rdfs_dy)
    if save_roots_boosted_selection:
        columns = list(rdf_data0.GetColumnNames())
        rdfs_dy.Snapshot("Events", "/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/OnlyPreProductionSelection/rootFilesBoostedSelection/dy_boosted_selection.root", columns)

    print("Initiating the conversion of RDataFrames to DataFrames for DY...")
    df_dy = pd.concat([
        pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        for rdf in rdfs_dy
    ], ignore_index=True)
    print("Conversion completed for DY.")
    print("Saving DataFrame of DY in pickle format...")
    df_dy.to_pickle(f"{output_dir}/dy.pkl")
    print("File saved: dy.pkl")
    print("Number of events dy: ", len(df_dy))
    del rdfs_dy  # Release memory
    if not OnlyGeneratePKLFiles:
        df_dy['eventWeight'] = df_dy['eventWeight'] * factor
    else:
        del df_dy
else:
    if (not OnlyGeneratePKLFiles):
        print("Loading DataFrame of DY from existing pickle file    ...")
        df_dy = pd.read_pickle(f"{output_dir}/dy.pkl")
        df_dy['eventWeight'] = df_dy['eventWeight'] * factor
        # df_dy.to_pickle(f"{output_dir}/dy.pkl")
        print("Loading completed for DY.")

# MC: VV
if is_vv:
    print("Creating RDataFrames for VV samples...")
    rdfs_vv = [
        create_rdf_with_weights(path, xs, luminosity)#, max_files=5)
        for path, xs in zip(input_dirs_vv, cross_sections_vv)
    ]

    print("Applying event selection to VV samples...")
    rdfs_vv = apply_event_selection(rdfs_vv)
    if save_roots_boosted_selection:
        columns = list(rdf_data0.GetColumnNames())
        rdfs_vv.Snapshot("Events", "/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/OnlyPreProductionSelection/rootFilesBoostedSelection/vv_boosted_selection.root", columns)

    print("Initiating the conversion of RDataFrames to DataFrames for VV...")
    df_vv = pd.concat([
        pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        for rdf in rdfs_vv
    ], ignore_index=True)
    print("Conversion completed for VV.")
    print("Saving DataFrame of VV in pickle format...")
    df_vv.to_pickle(f"{output_dir}/vv.pkl")
    print("File saved: vv.pkl")
    print("Number of events vv: ", len(df_vv))
    del rdfs_vv  # Release memory
    if not OnlyGeneratePKLFiles:
        df_vv['eventWeight'] = df_vv['eventWeight'] * factor
    else:
        del df_vv
else:
    if (not OnlyGeneratePKLFiles):
        print("Loading DataFrame of VV from existing pickle file...")
        df_vv = pd.read_pickle(f"{output_dir}/vv.pkl")
        df_vv['eventWeight'] = df_vv['eventWeight'] * factor
        # df_vv.to_pickle(f"{output_dir}/vv.pkl")
        print("Loading completed for VV.")

if OnlyGeneratePKLFiles:
    print("Only PKL files were generated. Exiting the program.")
    exit(0)


def convert_cpp_vectors_to_list(df):
    def to_list_safe(x):
        try:
            if hasattr(x, 'size') and hasattr(x, '__getitem__'):
                return [x[i] for i in range(x.size())]
            elif hasattr(x, "__iter__") and not isinstance(x, (str, bytes)):
                return list(x)
        except Exception:
            return x
        return x

    for col in df.columns:
        if df[col].dtype == 'object':
            df[col] = df[col].apply(to_list_safe)

    return df




df_data = convert_cpp_vectors_to_list(df_data)
print("DataFrame of Data converted to lists.")
df_ttbar = convert_cpp_vectors_to_list(df_ttbar)
print("DataFrame of TTBar converted to lists.")
df_ttbar_bck = convert_cpp_vectors_to_list(df_ttbar_bck)
print("DataFrame of TTBar Background converted to lists.")
df_singletop = convert_cpp_vectors_to_list(df_singletop)
print("DataFrame of SingleTop converted to lists.")
df_wjets = convert_cpp_vectors_to_list(df_wjets)
print("DataFrame of WJets converted to lists.")
df_qcdmu = convert_cpp_vectors_to_list(df_qcdmu)
print("DataFrame of QCDMu converted to lists.")
df_dy = convert_cpp_vectors_to_list(df_dy)
print("DataFrame of DY converted to lists.")
df_vv = convert_cpp_vectors_to_list(df_vv)
print("DataFrame of VV converted to lists.")

import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.ticker import ScalarFormatter
from matplotlib.gridspec import GridSpec

hep.style.use("CMS")

def save_histograms_to_root(hist_data, hist_mc, bin_edges, labels_mc, output_file):
    """
    Save histograms to a ROOT file.

    Args:
        hist_data (array): Data histogram.
        hist_mc (list of arrays): Monte Carlo histograms.
        bin_edges (array): Bin edges.
        labels_mc (list): Labels for the MC channels.
        output_file (str): Name of the ROOT file.
    """
    # Create the ROOT file
    root_file = ROOT.TFile(output_file, "RECREATE")

    # Save the data histogram
    h_data = ROOT.TH1F("Data", "Data", len(bin_edges) - 1, bin_edges)
    for i, value in enumerate(hist_data):
        h_data.SetBinContent(i + 1, value)
    h_data.Write()

    # Save the MC histograms
    for i, (hist, label) in enumerate(zip(hist_mc, labels_mc)):
        h_mc = ROOT.TH1F(f"MC_{label}", label, len(bin_edges) - 1, bin_edges)
        for j, value in enumerate(hist):
            h_mc.SetBinContent(j + 1, value)
        h_mc.Write()

    # Close the ROOT file
    root_file.Close()
    print(f"Histograms saved to: {output_file}")


def plot_variable(
    branch,
    dfs_mc,
    df_data,
    labels_mc,
    colors_mc,
    bins=30,
    xlim=None,
    ylim=None,
    logy=False,
    xlabel=None,
    ylabel="Events",
    title=None,
    print_percentages=False,
    output_dir=f"{output_dir}/plots{sufix_plotdir}",
    nth_element=None,  # For plotting the nth element of a list-like branch
    cut_id = None, #Adding cut id to distinguish different plots
):
    if title is None:
        title = branch
    def df_to_numpy(df, branch):
        if branch not in df.columns:
            return np.array([])

        data = df[branch].dropna()
        flat_data = []

        for item in data:
            try:
                # If it's iterable, we assume it may contain jets
                for val in item:
                    if isinstance(val, (int, float, np.integer, np.floating)):
                        flat_data.append(val)
            except TypeError:
                # Not iterable, but may be a single number
                if isinstance(item, (int, float, np.integer, np.floating)):
                    flat_data.append(item)
                # Otherwise, it's ignored

        return np.array(flat_data)
    def df_to_numpy_with_weights(df, branch, weight_branch="eventWeight"):
        if branch not in df.columns:
            return np.array([]), np.array([])

        data = df[[branch, weight_branch]].dropna()
        flat_data = []
        flat_weights = []

        for _, row in data.iterrows():
            jets = row[branch]
            weight = row[weight_branch]

            try:
                for jet in jets:
                    if isinstance(jet, (int, float, np.integer, np.floating)):
                        flat_data.append(jet)
                        flat_weights.append(weight)
            except TypeError:
                # jets not iterable → we ignore
                continue

        return np.array(flat_data), np.array(flat_weights)

    # Check that all bins are >= 0
    def verify_bins_non_negative(hist, label):
        if np.any(hist < 0):
            print(f"Warning: The histogram '{label}' contains bins with negative values.")
            print(f"Negative bins: {hist[hist < 0]}")
            # Optional: Adjust negative values to 0
            hist[hist < 0] = 0
        return hist
    # Create the output directory if it doesn't exist
    print("Creating output directories if they do not exist...")
    os.makedirs(output_dir, exist_ok=True)
    output_dir_hist = output_dir.replace("plots", "Histograms")
    os.makedirs(output_dir_hist, exist_ok=True)

    # Create figure and subplots
    print("Creating figure and subplots...")
    fig = plt.figure(figsize=(8, 12))
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.1)
    gs = GridSpec(3, 1, height_ratios=[3, 1, 1], hspace=0.1)
    ax_main = fig.add_subplot(gs[0])
    ax_fraction = fig.add_subplot(gs[1], sharex=ax_main)
    ax_ratio = fig.add_subplot(gs[2], sharex=ax_main)

    # Labels and luminosity
    print("Adding labels and luminosity...")
    ax_main.text(
        0, 1.05,
        f"Private work (CMS Data) {'cut=' + str(cut_id) if cut_id is not None else ''}",
        fontsize=24,
        verticalalignment='top',
        horizontalalignment='left',
        fontproperties="Tex Gyre Heros:italic",
        transform=ax_main.transAxes
    )
    ax_main.text(
        1, 1.05,
        f"{luminosity/1000:.2f} fb$^{{-1}}$",
        fontsize=22,
        verticalalignment='top',
        horizontalalignment='right',
        transform=ax_main.transAxes
    )

    # Prepare the data
    print("Preparing data...")
    bin_edges = np.linspace(xlim[0], xlim[1], bins + 1) if xlim else bins
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Filter the i-th element of the branch in the data
    print("Filtering data for the specified nth element...")
    if (nth_element is not None) and (nth_element > -1):
        df_data_filtered = df_data[branch].apply(lambda x: x[nth_element] if len(x) > nth_element else np.nan).dropna()
    elif (nth_element == -1):
        print(f"Using all elements in branch '{branch}'")
        df_data_filtered = df_to_numpy(df_data, branch)
    else:
        # print("Estoy aquiiiiiiii")
        df_data_filtered = df_data[branch]

    # Histogram of data
    print("Creating data histogram...")
    hist_data, _ = np.histogram(df_data_filtered, bins=bin_edges)
    errors_data = np.sqrt(hist_data)

    # Histograms of MC
    print("Creating MC histograms...")
    hist_mc = []
    mc_distributions = []
    mc_weights = []
    for df, label, color in zip(dfs_mc, labels_mc, colors_mc):
        weights = df["eventWeight"] if "eventWeight" in df.columns else np.ones(len(df))
        # print(f"Length of {label} DataFrame: {len(df)}")
        # print(f"Number of weights for {label}: {len(weights)}")
        # print(f"Shape of {label} DataFrame: {df.shape}")
        # print(f"Shape of weights for {label}: {weights.shape}")
        if (nth_element is not None) and (nth_element > -1):
            mc_filtered = df[branch].apply(lambda x: x[nth_element] if len(x) > nth_element else np.nan).dropna()
            weights_filtered = weights[df[branch].apply(lambda x: len(x) > nth_element)]
        elif nth_element == -1:
            mc_filtered, weights_filtered = df_to_numpy_with_weights(df, branch, weight_branch="eventWeight")
        else:
            mc_filtered = df[branch]
            weights_filtered = weights

        hist, _ = np.histogram(mc_filtered, bins=bin_edges, weights=weights_filtered)
        # Check that histogram values are >= 0
        # hist = verify_bins_non_negative(hist, label)
        hist_mc.append(hist)
        mc_distributions.append(mc_filtered)
        mc_weights.append(weights_filtered)

    # Save histograms to a ROOT file
    print("Saving histograms to ROOT file...")
    root_output_file = os.path.join(output_dir_hist, f"{title}.root")
    save_histograms_to_root(hist_data, hist_mc, bin_edges, labels_mc, root_output_file)
    # print(f"number of elements in ttbar: {len(dfs_mc[0])}")
    # print(f"number of elements in weights: {len(mc_weights[0])}")
    # print(f"shape of ttbar DataFrame: {dfs_mc[0].shape}")
    # print(f"shape of weights for ttbar: {mc_weights[0].shape}")
    # print("\n")

    # Draw stacked MC histograms
    print("Drawing stacked MC histograms...")
    ax_main.hist(
        [bin_edges[:-1]] * len(hist_mc),  # Use bin edges for each histogram
        bins=bin_edges,
        weights=hist_mc,  # Use corrected histograms
        stacked=True,
        label=labels_mc,
        color=colors_mc,
        alpha=0.7,
    )

    # Configure the main plot
    print("Configuring main plot...")
    ax_main.errorbar(
        bin_centers,
        hist_data,
        yerr=errors_data,
        fmt="o",
        label="Data",
        color="black",
    )
    ax_main.set_ylabel(ylabel)
    if xlim:
        ax_main.set_xlim(*xlim)
    if ylim:
        ax_main.set_ylim(*ylim)
    if logy:
        ax_main.set_yscale("log")
    else:
        ax_main.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax_main.get_yaxis().get_offset_text().set_x(-0.13)
        ax_main.get_yaxis().get_offset_text().set_y(1.25)
        ax_main.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    #set always ymin = 0
    ax_main.set_ylim(bottom=0)
    handles, labels = ax_main.get_legend_handles_labels()
    ax_main.legend(handles[::-1], labels[::-1], frameon=False, fontsize=16) #loc="upper right", bbox_to_anchor=(1,1), borderaxespad=0,
    ax_main.grid(True, linestyle="--", alpha=0.5)
    ax_main.tick_params(axis='x', labelbottom=False)
    # if title:
    #     ax_main.set_title(title)

    # Calculate fractions for the fraction subplot
    print("Calculating fractions for the fraction subplot...")
    total_mc = np.sum(hist_mc, axis=0)
    total_mc[total_mc == 0] = 1  # Avoid division by zero
    fractions = [hist / total_mc for hist in hist_mc]

    # Draw fractions in the fraction subplot (stacked)
    print("Drawing fraction subplot...")
    ax_fraction.hist(
        [bin_centers] * len(fractions),
        bins=bin_edges,
        weights=fractions,
        stacked=True,
        color=colors_mc,
        alpha=0.7,
    )
    ax_fraction.set_ylabel("Fraction")
    ax_fraction.set_ylim(0, 1)
    ax_fraction.grid(True, linestyle="--", alpha=0.5)
    ax_fraction.tick_params(axis='x', labelbottom=False)

    # Calculate Data/MC ratio for the ratio subplot
    print("Calculating Data/MC ratio for the ratio subplot...")
    total_mc_sum = np.sum(hist_mc, axis=0)
    # total_mc_sum[total_mc_sum == 0] = -999  # Avoid division by zero
    ratio = np.where(total_mc_sum > 0, (hist_data) / total_mc_sum, np.nan)
    ratio_err = np.where(total_mc_sum > 0, abs(errors_data / total_mc_sum), np.nan)

    # Draw ratio in the ratio subplot
    print("Drawing ratio subplot...")
    ax_ratio.errorbar(
        bin_centers,
        ratio,
        yerr=ratio_err,
        fmt="o",
        color="black",
    )
    ax_ratio.axhline(1, color="red", linestyle="--", linewidth=1)
    ax_ratio.set_ylabel(r'$\frac{\mathrm{Data}}{\mathrm{MC}}$')
    ax_ratio.set_xlabel(xlabel or branch)
    ax_ratio.set_ylim(0.5, 1.5)
    ax_ratio.grid(True, linestyle="--", alpha=0.5)

    # Save the plot as a file
    print("Saving the plot...")
    output_file = os.path.join(output_dir, f"{title}.png")
    # hep.cms.label("Private", loc=0, rlabel=f"{luminosity/1000:.2f} fb$^{{-1}}$", ax=ax_main, fontsize=22)
    # plt.tight_layout()
    plt.savefig(output_file)
    print(f"Plot saved to: {output_file}")

    # Close the plot to free memory
    print("Closing the plot to free memory...")
    plt.close()

    if print_percentages:
        total_events_per_mc = [np.sum(hist) for hist in hist_mc]
        total_events = sum(total_events_per_mc)
        percentages = [100 * val / total_events for val in total_events_per_mc]

        print("\nContribución de cada canal MC:")
        for label, events, percent in zip(labels_mc, total_events_per_mc, percentages):
            print(f"{label}: {events:.1f} eventos ({percent:.2f}%)")
        print(f"Total MC (ponderado): {total_events:.1f} eventos\n")

df_data_muon = df_data
print("Number of events in df_data after conversion: ", len(df_data_muon))
df_ttbar_muon = df_ttbar
print("Number of events in df_ttbar after conversion: ", len(df_ttbar_muon))
df_ttbar_bck_muon = df_ttbar_bck
print("Number of events in df_ttbar_bck after conversion: ", len(df_ttbar_bck_muon))
df_singletop_muon = df_singletop
print("Number of events in df_singletop after conversion: ", len(df_singletop_muon))
df_wjets_muon = df_wjets
print("Number of events in df_wjets after conversion: ", len(df_wjets_muon))
df_qcdmu_muon = df_qcdmu
print("Number of events in df_qcdmu after conversion: ", len(df_qcdmu_muon))
df_dy_muon = df_dy
print("Number of events in df_dy after conversion: ", len(df_dy_muon))
df_vv_muon = df_vv
print("Number of events in df_vv after conversion: ", len(df_vv_muon))


# print("\n\nApplying tight b-tagging ...")
# # Filtrar los DataFrames para muones
# df_data_muon = df_data[
#     (abs(df_data['lepton_trg_pdgId']) == 13) &
#     (df_data['lepton_trg_pt'] > 55) &
#     (abs(df_data['lepton_trg_eta']) < 2.4) &
# #     (
# #         df_data['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
# #         df_data['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
# #         == 1
# #     )
#     (
#         df_data.apply( lambda row: sum(pt > 30 and abs(eta) < 2.4 and btag > 0.6553 and passjetid
#                                               for pt, eta, btag, passjetid in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'], row['Jet_passJetIdTightLepVeto'])) >= 1, axis=1)
#     )
#  ]
# print("Number of events in data muon: ", len(df_data_muon))

# df_ttbar_muon = df_ttbar[
#     (abs(df_ttbar['lepton_trg_pdgId']) == 13) &
#     (df_ttbar['lepton_trg_pt'] > 55) &
#     (abs(df_ttbar['lepton_trg_eta']) < 2.4) &
# #     (
# #         df_ttbar['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
# #         df_ttbar['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
# #         == 1
# #     )
#     (
#         df_ttbar.apply( lambda row: sum(pt > 30 and abs(eta) < 2.4 and btag > 0.6553 and passjetid
#                                               for pt, eta, btag, passjetid in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'], row['Jet_passJetIdTightLepVeto'])) >= 1, axis=1)
#     )
#  ]
# print("Number of events in ttbar muon: ", len(df_ttbar_muon))

# df_ttbar_bck_muon = df_ttbar_bck[
#     (abs(df_ttbar_bck['lepton_trg_pdgId']) == 13) &
#     (df_ttbar_bck['lepton_trg_pt'] > 55) &
#     (abs(df_ttbar_bck['lepton_trg_eta']) < 2.4) &
# #     (
# #         df_ttbar_bck['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
# #         df_ttbar_bck['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
# #         == 1
# #     )
#     (
#         df_ttbar_bck.apply( lambda row: sum(pt > 30 and abs(eta) < 2.4 and btag > 0.6553 and passjetid
#                                               for pt, eta, btag, passjetid in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'], row['Jet_passJetIdTightLepVeto'])) >= 1, axis=1)
#     )
#  ]
# print("Number of events in ttbar background muon: ", len(df_ttbar_bck_muon))

# df_singletop_muon = df_singletop[
#     (abs(df_singletop['lepton_trg_pdgId']) == 13) &
#     (df_singletop['lepton_trg_pt'] > 55) &
#     (abs(df_singletop['lepton_trg_eta']) < 2.4) &
# #     (
# #         df_singletop['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
# #         df_singletop['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
# #         == 1
# #     )
#     (
#         df_singletop.apply( lambda row: sum(pt > 30 and abs(eta) < 2.4 and btag > 0.6553 and passjetid
#                                               for pt, eta, btag, passjetid in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'], row['Jet_passJetIdTightLepVeto'])) >= 1, axis=1)
#     )
#  ]
# print("Number of events in single top muon: ", len(df_singletop_muon))

# df_wjets_muon = df_wjets[
#     (abs(df_wjets['lepton_trg_pdgId']) == 13) &
#     (df_wjets['lepton_trg_pt'] > 55) &
#     (abs(df_wjets['lepton_trg_eta']) < 2.4) &
# #     (
# #         df_wjets['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
# #         df_wjets['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
# #         == 1
# #     )
#     (
#         df_wjets.apply( lambda row: sum(pt > 30 and abs(eta) < 2.4 and btag > 0.6553 and passjetid
#                                               for pt, eta, btag, passjetid in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'], row['Jet_passJetIdTightLepVeto'])) >= 1, axis=1)
#     )
#  ]
# print("Number of events in W+jets muon: ", len(df_wjets_muon))

# df_qcdmu_muon = df_qcdmu[
#     (abs(df_qcdmu['lepton_trg_pdgId']) == 13) &
#     (df_qcdmu['lepton_trg_pt'] > 55) &
#     (abs(df_qcdmu['lepton_trg_eta']) < 2.4) &
# #     (
# #         df_qcdmu['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
# #         df_qcdmu['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
# #         == 1
# #     )
#     (
#         df_qcdmu.apply( lambda row: sum(pt > 30 and abs(eta) < 2.4 and btag > 0.6553 and passjetid
#                                               for pt, eta, btag, passjetid in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'], row['Jet_passJetIdTightLepVeto'])) >= 1, axis=1)
#     )
#  ]
# print("Number of events in QCD muon: ", len(df_qcdmu_muon))

# df_dy_muon = df_dy[
#     (abs(df_dy['lepton_trg_pdgId']) == 13) &
#     (df_dy['lepton_trg_pt'] > 55) &
#     (abs(df_dy['lepton_trg_eta']) < 2.4) &
# #     (
# #         df_others['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
# #         df_others['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
# #         == 1
# #     )
#     (
#         df_dy.apply( lambda row: sum(pt > 30 and abs(eta) < 2.4 and btag > 0.6553 and passjetid
#                                               for pt, eta, btag, passjetid in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'], row['Jet_passJetIdTightLepVeto'])) >= 1, axis=1)
#     )
#  ]
# print("Number of events in DY processes muon: ", len(df_dy_muon))

# df_vv_muon = df_vv[
#     (abs(df_vv['lepton_trg_pdgId']) == 13) &
#     (df_vv['lepton_trg_pt'] > 55) &
#     (abs(df_vv['lepton_trg_eta']) < 2.4) &
# #     (
# #         df_others['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
# #         df_others['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
# #         == 1
# #     )
#     (
#         df_vv.apply( lambda row: sum(pt > 30 and abs(eta) < 2.4 and btag > 0.6553 and passjetid
#                                               for pt, eta, btag, passjetid in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'], row['Jet_passJetIdTightLepVeto'])) >= 1, axis=1)
#     )
#  ]
# print("Number of events in VV processes muon: ", len(df_vv_muon))

print("Plots")

#LEPTONS
# #Lepton pt
# plot_variable( branch="lepton_trg_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
#               bins=70, logy=False, xlim=(0, 400), #ylim=(0.75, 1e3),
#               xlabel=r'$p_T$ muon [GeV]', title='muon_pt'
#              )

# #Lepton TuneP pt              
# plot_variable( branch="lepton_tuneP_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
#               bins=70, logy=False, xlim=(0, 400), #ylim=(0.75, 1e3),
#               xlabel=r'$p_T$ muon [GeV]', title='muon_tuneP_pt'
#              )

# #Lepton eta
# plot_variable(branch="lepton_trg_eta", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
#               bins=15, xlim=(-3, 3), #ylim=(0, 0.54*1e3), #logy=True,
#               xlabel=r'$\eta$ muon', title='muon_eta'
#              )

# Now we create a function to filter one by one the different cuts:
# The cuts ids are:
# 0: PreProductionSelection: 
# 1: Trigger
# 2: PU reweight
# 3: MET filter
# 4: Exactly one lepton
# 5: MET > 50
# 6: n_bjets >= 1
# 7: TopJet.pt > 350
# 8: N_subjets == 3
# 9: Top pt reweight
# 10: Muon corrections
# 11: btagging corrections
# By indicating the array position we filter until that cut (including it)


def filter_by_cut(df, cut_id):
    if cut_id < 0 or cut_id > 11:
        raise ValueError("cut_id must be between 0 and 11")
    if cut_id == 0:
        return df
    if cut_id > 0:
        df = df[
                (df['is_muon_trg'] == True) & 
                ((df['HLT_Mu50']) & (df['Jet_vetoMap'])) #"HLT_HighPtTkMu100", "HLT_CascadeMu100", "HLT_IsoMu24", "HLT_Mu50"
            ]
    if cut_id > 1:
        if 'PUReweighting' in df.columns:
            mask = np.isfinite(df['PUReweighting'])
            df.loc[mask, 'eventWeight'] = df.loc[mask, 'PUReweighting'] * df.loc[mask, 'eventWeight']
    if cut_id > 2:
        df = df[df['pass_MET_filters_2022'] == True]
    if cut_id > 3:
        df = df[
                (df['num_muon_selected'] == 1) &
                (df['Num_lep_above_15'] == 1) &
                (df['is_muon_trg'] == True)
            ]
    if cut_id > 4:
        df = df[df['PuppiMET_pt'] > 50]
    if cut_id > 5:
        df = df[df['n_bjets'] >= 1]
    if cut_id > 6:
        df = df[df['topjets.pt'] > 350]
    if cut_id > 7:
        df = df[
            (df['num_subjets_in_topjet_pt30_eta2p5'] >= 3) &
            (df['topjets.n_subjets'] == 3)
        ]
    if cut_id > 8:
        if 'TopPtWeight_NNLOpNLOEW' in df.columns:
            mask = np.isfinite(df['TopPtWeight_NNLOpNLOEW'])
            df.loc[mask, 'eventWeight'] = df.loc[mask, 'TopPtWeight_NNLOpNLOEW'] * df.loc[mask, 'eventWeight']
    if cut_id > 9:
        if 'muoWeight' in df.columns:
            mask = np.isfinite(df['muoWeight'])
            df.loc[mask, 'eventWeight'] = df.loc[mask, 'muoWeight'] * df.loc[mask, 'eventWeight']
    if cut_id > 10:
        if 'btagWeight' in df.columns:
            mask = np.isfinite(df['btagWeight'])
            df.loc[mask, 'eventWeight'] = df.loc[mask, 'btagWeight'] * df.loc[mask, 'eventWeight']
    return df


# Dictionary to map cut ids to descriptive names
cut_names = {
    0: "PreProductionSelection",
    1: "Trigger",
    2: "PUreweight",
    3: "METfilters",
    4: "Exactly1Muon",
    5: "PuppiMET_gt_50",
    6: "n_bjets_ge_1",
    7: "TopJet_pt_gt_350",
    8: "N_subjets_eq_3",
    9: "TopPtWeight_NNLOpNLOEW",
    10: "muoWeight",
    11: "btagWeight"
}

for cut_id in range(12):
    # if cut_id == 8: continue  # Skip cut 8 as per instruction
    df_data_muon = filter_by_cut(df_data, cut_id)
    df_ttbar_muon = filter_by_cut(df_ttbar, cut_id)
    df_ttbar_bck_muon = filter_by_cut(df_ttbar_bck, cut_id)
    df_singletop_muon = filter_by_cut(df_singletop, cut_id)
    df_wjets_muon = filter_by_cut(df_wjets, cut_id)
    df_qcdmu_muon = filter_by_cut(df_qcdmu, cut_id)
    df_dy_muon = filter_by_cut(df_dy, cut_id)
    df_vv_muon = filter_by_cut(df_vv, cut_id)

    #Lepton pt
    plot_variable( branch="lepton_trg_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
                labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], cut_id=cut_id,
                bins=70, logy=False, xlim=(0, 400), #ylim=(0.75, 1e3),
                xlabel=r'$p_T$ muon [GeV]', title=f'muon_pt_{cut_id}_{cut_names[cut_id]}'
                )

    #Lepton TuneP pt              
    plot_variable( branch="lepton_tuneP_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
                labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], cut_id=cut_id,
                bins=70, logy=False, xlim=(0, 400), #ylim=(0.75, 1e3),
                xlabel=r'$p_T$ muon TuneP [GeV]', title=f'muon_tuneP_pt_{cut_id}_{cut_names[cut_id]}'
                )

    #Lepton eta
    plot_variable(branch="lepton_trg_eta", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
                labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], cut_id=cut_id,
                bins=15, xlim=(-3, 3), #ylim=(0, 0.54*1e3), #logy=True,
                xlabel=r'$\eta$ muon', title=f'muon_eta_{cut_id}_{cut_names[cut_id]}'
                )

    #Number of leptons above 15 GeV
    plot_variable(branch="Num_lep_above_15", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
                labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], cut_id=cut_id,
                bins=np.arange(-0.5, 5.5, 1), #ylim=(0, 3.5*1e3), #logy=True,
                xlabel='Number of leptons above 15 GeV', title=f'n_lep_{cut_id}_{cut_names[cut_id]}'
                )


    #Number of bjets
    plot_variable(branch="n_bjets", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
                labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], cut_id=cut_id,
                bins=np.arange(-0.5, 5.5, 1), #ylim=(0, 3.5*1e3), #logy=True,
                xlabel='Number of b-jets', title=f'n_bJets_{cut_id}_{cut_names[cut_id]}'
                )

    #Leading bjet pt
    plot_variable(branch="bjet_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
                labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], cut_id=cut_id,
                nth_element=0,  # Select the leading b-jet
                bins=70, xlim=(0, 400), #ylim=(0, 0.54*1e3), #logy=True,
                xlabel='Leading b-jet $p_T$ [GeV]', title=f'leading_bjet_pt_{cut_id}_{cut_names[cut_id]}'
                )

    #Number of jets
    plot_variable(branch="n_jets", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
                labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], cut_id=cut_id,
                bins=np.arange(-0.5, 10.5, 1), #ylim=(0, 3.5*1e3), #logy=True,
                xlabel='Number of jets', title=f'n_jets_{cut_id}_{cut_names[cut_id]}'
                )
    

    #leading jet pt
    plot_variable(branch="Jet_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
                labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], cut_id=cut_id,
                nth_element=0,  # Select the leading jet
                bins=95, xlim=(0, 600), #ylim=(0, 0.54*1e3), #logy=True,
                xlabel='Leading jet $p_T$ [GeV]', title=f'leading_jet_pt_{cut_id}_{cut_names[cut_id]}'
                )

    #leading bjet eta
    plot_variable(branch="bjet_eta", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
                labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], cut_id=cut_id,
                nth_element=0,  # Select the leading b-jet
                bins=15, xlim=(-3, 3), #ylim=(0, 0.54*1e3), #logy=True,
                xlabel='Leading b-jet $\eta$', title=f'leading_bjet_eta_{cut_id}_{cut_names[cut_id]}'
                )
    
    #leading jet eta
    plot_variable(branch="Jet_eta", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
                labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], cut_id=cut_id,
                nth_element=0,  # Select the leading jet
                bins=15, xlim=(-3, 3), #ylim=(0, 0.54*1e3), #logy=True,
                xlabel='Leading jet $\eta$', title=f'leading_jet_eta_{cut_id}_{cut_names[cut_id]}'
                )

    #Number of good subjets in topjet
    plot_variable(branch="num_subjets_in_topjet_pt30_eta2p5", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
                labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], cut_id=cut_id,
                bins=np.arange(-0.5, 10.5, 1), #ylim=(0, 3.5*1e3), #logy=True,
                xlabel='Number of good subjets in topjet', title=f'n_goodTopSubjets_{cut_id}_{cut_names[cut_id]}'
                )

    # Topjet pt
    plot_variable(branch="topjets.pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
                labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], cut_id=cut_id,
                bins=70, xlim=(0, 800), #ylim=(0, 0.54*1e3), #logy=True,
                xlabel='Topjet $p_T$ [GeV]', title=f'topjet_pt_{cut_id}_{cut_names[cut_id]}'
                )

    # Topjet eta
    plot_variable(branch="topjets.eta", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
                labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], cut_id=cut_id,
                bins=15, xlim=(-3, 3), #ylim=(0, 0.54*1e3), #logy=True,
                xlabel='Topjet $\eta$', title=f'topjet_eta_{cut_id}_{cut_names[cut_id]}'
                )

    # Topjet mass
    plot_variable(branch="topjets.mass", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
                labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], cut_id=cut_id,
                bins=70, xlim=(0, 300), #ylim=(0, 0.54*1e3), #logy=True,
                xlabel='Topjet mass [GeV]', title=f'topjet_mass_{cut_id}_{cut_names[cut_id]}'
                )

    # Leading subjet pt
    plot_variable(branch="subjets.pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
                labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], cut_id=cut_id,
                nth_element=0,  # Select the leading subjet
                bins=95, xlim=(0, 600), #ylim=(0, 0.54*1e3), #logy=True,
                xlabel='Leading subjet $p_T$ [GeV]', title=f'leading_subjet_pt_{cut_id}_{cut_names[cut_id]}'
                )

    # Leading subjet eta
    plot_variable(branch="subjets.eta", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
                labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], cut_id=cut_id,
                nth_element=0,  # Select the leading subjet
                bins=15, xlim=(-3, 3), #ylim=(0, 0.54*1e3), #logy=True,
                xlabel='Leading subjet $\eta$', title=f'leading_subjet_eta_{cut_id}_{cut_names[cut_id]}'
                )

    #MET pt
    plot_variable(branch="PuppiMET_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
                labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], cut_id=cut_id,
                bins=95, xlim=(0, 600), #ylim=(0, 0.54*1e3), #logy=True,
                xlabel='Puppi MET $p_T$ [GeV]', title=f'PuppiMET_pt_{cut_id}_{cut_names[cut_id]}'
                )



# cut_id = 11  # Example cut_id for final selection after all cuts
# df_data_muon = filter_by_cut(df_data, cut_id)
# df_ttbar_muon = filter_by_cut(df_ttbar, cut_id)
# df_ttbar_bck_muon = filter_by_cut(df_ttbar_bck, cut_id)
# df_singletop_muon = filter_by_cut(df_singletop, cut_id)
# df_wjets_muon = filter_by_cut(df_wjets, cut_id)
# df_qcdmu_muon = filter_by_cut(df_qcdmu, cut_id)
# df_dy_muon = filter_by_cut(df_dy, cut_id)
# df_vv_muon = filter_by_cut(df_vv, cut_id)

# columns_mc = ["eventWeight", "lepton_trg_pt", "lepton_tuneP_pt", "lepton_trg_eta", "lepton_trg_jet_pt", "pass_trigger", "PUReweighting", "pass_MET_filters_2023", "Num_lep_above_15", "num_muon_selected", "is_muon_trg", "PuppiMET_pt", "n_bjets", "bjet_pt", "bjet_eta", "n_jets", "TopPtWeight_NNLOpNLOEW", "muoWeight", "btagWeight", "topjets.pt", "topjets.eta", "topjets.mass", "num_subjets_in_topjet_pt30_eta2p5", "topjets.n_subjets", "subjets.pt", "subjets.eta", "Jet_pt", "Jet_eta", "Jet_hadronFlavour", "Jet_btagDeepFlavB", "Pileup_nTrueInt"]

# #LEPTONS
# #Numer of leptons
# plot_variable( branch="lepton_trg_n_lep", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
#               bins=np.arange(-0.5, 3.5, 1), print_percentages=True, #ylim=(0, 3.5*1e3), #logy=True,)
#               xlabel='Number of muons', title='lepton_trg_n_lep'
#              )

# #Lepton pdgId
# plot_variable(branch="lepton_trg_pdgId", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
#               bins=np.arange(-14.5, 14.5, 1),logy=False, #ylim=(0, 1.85*1e3), #logy=True,)
#               xlabel='Lepton pdgID', title='lepton_trg_pdgId'
#              )

# #Lepton pt
# plot_variable( branch="lepton_trg_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
#               bins=70, logy=False, xlim=(0, 400), #ylim=(0.75, 1e3),
#               xlabel=r'$p_T$ muon [GeV]', title='lepton_trg_pt'
#              )

# #Lepton TuneP pt              
# plot_variable( branch="lepton_tuneP_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
#               bins=70, logy=False, xlim=(0, 400), #ylim=(0.75, 1e3),
#               xlabel=r'$p_T$ muon TuneP [GeV]', title='lepton_tuneP_pt'
#              )

# #Lepton eta
# plot_variable(branch="lepton_trg_eta", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
#               bins=15, xlim=(-3, 3), #ylim=(0, 0.54*1e3), #logy=True,
#               xlabel=r'$\eta$ muon', title='lepton_trg_eta'
#              )

# #Lepton dR to jet
# plot_variable(branch="lepton_trg_dR", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
#               bins=20, xlabel=r"$\Delta R$ (muon, closest AK4 jet)", xlim=(0, 4), logy=False, #ylim=(0, 0.8*1e3),
#               title='lepton_trg_dR'
#              )

# #Lepton ptrel to jet
# plot_variable(branch="lepton_trg_ptrel", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
#               bins=20, xlabel=r"$p_T^{rel}$ (muon, closest AK4 jet) [GeV]", xlim=(0, 200), logy = False, #ylim=(0, 0.35*1e3),
#               title='lepton_trg_ptrel'
#              )

# #XCONE JETS
# #Mass top hadronic decay
# plot_variable(branch="topjets.mass", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
#              bins=20, xlim=(80, 320), logy=False, #ylim=(0, 0.64*1e3),
#              title='topjets_mass', xlabel=r'$m_{top-jet}$ [GeV]'
#              )

# #Mass W hadronic decay
# plot_variable(branch="topjets.mass_W", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
#              bins=20, xlim=(40, 120), logy=False, #ylim=(0, 0.6*1e3),
#              title='topjets_Wmass', xlabel=r'$m_{W-jet}$ [GeV]'
#              )

# #MET
# #Missing transverse energy pt
# plot_variable(branch="PuppiMET_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
#               bins=25, xlim=(0, 600), logy=False, #ylim=(0.75, 1e3),
#               xlabel=r'$p_T^{MET}$ [GeV]', title='PuppiMET_pt'
#               )

# #BJETS
# #Numer of bjets
# plot_variable( branch="n_bjets", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
#               bins=np.arange(-0.5, 6.5, 1), logy=False, #ylim=(0.7, 2.5*1e3), #logy=True,)
#               xlabel='Number of b-jets', title='n_bjets'
#              )

# # #bjets btag discriminator
# plot_variable(branch="bjet_btag", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=-1,
#               bins=20, xlim=(0., 1.),logy=False, #ylim=(0.95, 4*1e3), #logy=True,)
#               xlabel='DeepBTag b-jets', title='bjet_btag'
#              )

# #bjets pt
# plot_variable( branch="bjet_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=-1,
#               bins=20, logy=False, xlim=(0, 600), #ylim=(0.75, 4*1e3),
#               xlabel=r'$p_T$ b-jets [GeV]', title='bjet_pt'
#              )

# #bjets eta
# plot_variable(branch="bjet_eta", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=-1,
#               bins=15, xlim=(-3, 3), #ylim=(0, 0.9*1e3), #logy=True,
#               xlabel=r'$\eta$ b-jets', title='bjet_eta'
#              )

# #LEADDING BJET
# #Leading btag discriminator
# plot_variable(branch="bjet_btag", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=0,
#               bins=20, xlim=(0., 1.),logy=False, #ylim=(0.95, 4.85*1e3), #logy=True,)
#               xlabel='DeepBTag leading b-jet', title='bjet_btag_leading'
#              )

# #Leading bjet pt
# plot_variable( branch="bjet_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=0,
#               bins=20, logy=False, xlim=(0, 600), #ylim=(0.75, 4*1e3),
#               xlabel=r'$p_T$ leading b-jet [GeV]', title='bjet_pt_leading'
#              )

# #Leading bjet eta
# plot_variable(branch="bjet_eta", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=0,
#               bins=15, xlim=(-3, 3), #ylim=(0, 6.6*1e2), #logy=True,
#               xlabel=r'$\eta$ leading b-jet', title='bjet_eta_leading'
#              )

# #SUBLEADING BJET
# #Subleading btag discriminator
# plot_variable(branch="bjet_btag", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=1,
#               bins=20, xlim=(0., 1.),logy=False, #ylim=(0.95, 1*1e3), #logy=True,)
#               xlabel='DeepBTag subleading b-jet', title='bjet_btag_subleading'
#              )

# #Subleading bjet pt
# plot_variable( branch="bjet_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=1,
#               bins=20, logy=False, xlim=(0, 300), #ylim=(0.75, 1*1e3),
#               xlabel=r'$p_T$ subleading b-jet [GeV]', title='bjet_pt_subleading'
#              )

# #Subleading bjet eta
# plot_variable(branch="bjet_eta", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=1,
#               bins=15, xlim=(-3, 3), #ylim=(0, 2.05*1e2), #logy=True,
#               xlabel=r'$\eta$ subleading b-jet', title='bjet_eta_subleading'
#              )

# #AK4 JETS
# #AK4 jets pt
# plot_variable( branch="Jet_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=-1,
#               bins=20, xlim=(0, 800), logy=False, #ylim=(0.7, 1.1*1e4),
#               xlabel=r'$p_T$ AK4 jets [GeV]', title='AK4Jet_pt'
#              )

# plot_variable( branch="fatjets.pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=-1,
#               bins=20, xlim=(0, 800), logy=False, #ylim=(0.7, 1.1*1e4),
#               xlabel=r'$p_T$ Fat-jets [GeV]', title='fatjets_pt'
#              )  