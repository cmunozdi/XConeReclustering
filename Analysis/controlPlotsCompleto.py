# . /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc11-opt/setup.sh
import ROOT
import numpy as np
import pandas as pd
import os
from glob import glob

#Bool variables to know over which sample should we run the code
is_ttbar_semi = True
is_ttbar_bck = True
is_singletop = True
is_wjets = True
is_qcdmu = True
is_dy = True
is_vv = True
is_data0 = True
is_data1 = True



# Charging header (adjust path if necessary)
ROOT.gInterpreter.Declare('#include "./btagSFs/header_btag_scaleFactors.h"')
# Initializing CorrectionSet ONLY ONCE
ROOT.gInterpreter.ProcessLine("initializeBTagCorrectionSet();")
ROOT.gInterpreter.ProcessLine("initializeMUOCorrectionSet();")
ROOT.gInterpreter.ProcessLine("initializePUReweightingCorrectionSet();")

ROOT.ROOT.EnableImplicitMT()

def get_pattern(input_dir):
    if any(keyword in input_dir for keyword in ["wjets", "ttbar"]):
        return "**/250818*/**/*.root"  # Specific pattern for wjets and ttbar
    else:
        return "**/250*/**/*.root"  # General pattern for other directories

def create_rdf_with_weights(input_dir, cross_section, luminosity, max_files=None, isTTbar=False):
    # Choose the pattern depending on whether it is QCDMu or not
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


luminosity = 18.062658998*1000 #This new value is using /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json instead of /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json 18.084440726*1000 #7.229453396*1000 #18.083517794*1000 pb^-1

cross_sections_singletop = [87.200000, 4.534000, 4.663095, 19.970880, 19.302840, 145.000000, 7.244000, 4.663095, 19.970880, 19.302840] #[87.200000, 1.477177, 19.302840, 145.000000, 2.360095, 19.302840] # #
cross_sections_wjets = [368.200000, 421.900000, 25.600000, 54.770000, 0.878500, 3.119000, 4427.000000, 1598.000000, 0.105300, 0.526100]
cross_sections_ttbar = [923.6*0.4392] #July 2022 value: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO#Updated_reference_cross_sections
cross_sections_ttbar_bck = [923.6*0.1061, 923.6*0.4544]
cross_sections_qcdmu = [7763, 699.1, 68.24, 21.37, 3.913, 1.323] #[1.323000, 23240.000000, 7763.000000, 699.100000, 68.240000, 404400.000000, 21.370000, 3.913000, 95900.000000]   #
cross_sections_dy = [6688] #[141500, 20950, 6688] #[20950, 141500, 6688] #
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

        #Full production samples
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/SingleTop/TbarBQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/SingleTop/TbarBtoLminusNuB-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/SingleTop/TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/SingleTop/TbarWplusto4Q_TuneCP5_13p6TeV_powheg-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/SingleTop/TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/SingleTop/TBbarQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/SingleTop/TBbartoLplusNuBbar-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/SingleTop/TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/SingleTop/TWminusto4Q_TuneCP5_13p6TeV_powheg-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/SingleTop/TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",


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

        #Full production samples
        # "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-200to400_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-400to600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-400to600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-40to100_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-40to100_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",



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

        #Full production samples
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/TTbar/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",

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

        #Full production samples
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/TTbar/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/TTbar/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8",


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

        #Full production samples
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/QCDMu/QCD_PT-170to300_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/QCDMu/QCD_PT-300to470_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/QCDMu/QCD_PT-470to600_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/QCDMu/QCD_PT-600to800_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/QCDMu/QCD_PT-800to1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/QCDMu/QCD_PT-1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/QCDMu/QCD_PT-120to170_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/QCDMu/QCD_PT-80to120_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",



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

        #Full production samples
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/DY/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",


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

        #Full production samples
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/VV/WW_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/VV/WZ_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2023preBPix/VV/ZZ_TuneCP5_13p6TeV_pythia8"
    ]

print("Iniciando la conversión de RDataFrames a DataFrames.")



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


def apply_event_selection(rdfs, isMC=True, TWPbtag2023=0.6553):
    filtered_rdfs = []

    for rdf in rdfs:
        branches = set(rdf.GetColumnNames())
        if "Num_lep_above_15" in branches:
            rdf = rdf.Redefine("Num_lep_above_15", "Sum(Muon_pt > 15 && abs(Muon_eta) < 2.4) + Sum(Electron_pt > 15 && abs(Electron_eta) < 2.4)")
        else:
            rdf = rdf.Define("Num_lep_above_15", "Sum(Muon_pt > 15 && abs(Muon_eta) < 2.4) + Sum(Electron_pt > 15 && abs(Electron_eta) < 2.4)")


        rdf_filtered = (
            rdf
            # 1. Detector selection
            .Filter("pass_detector_selection", "Detector selection") #n_leptons == 1 && n_jets > 0 && pt_miss > 37.5 && n_fatjets > 0
#             .Filter("lepton.n_lep==1")

            # 2. Lepton selection
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
            .Filter('is_muon_trg') \
            # .Filter('(lepton.dR_to_jet[selected_lepton_idx] > 3.)&&(lepton.dR_to_jet[selected_lepton_idx] < 3.2)') \
            .Define('lepton_trg_pt', 'lepton.pt[selected_lepton_idx]') \
            .Define('lepton_trg_eta', 'lepton.eta[selected_lepton_idx]') \
            .Define('lepton_trg_phi', 'lepton.phi[selected_lepton_idx]') \
            .Define('lepton_trg_pdgId', 'lepton.pdgId[selected_lepton_idx]') \
            .Define('lepton_trg_dR', 'lepton.dR_to_jet[selected_lepton_idx]') \
            .Define('lepton_trg_ptrel', 'lepton.pt_rel_to_jet[selected_lepton_idx]') \
            .Define('num_muon_selected', 'Sum((Muon_pt == lepton_trg_pt) && (Muon_eta == lepton_trg_eta) && (Muon_phi == lepton_trg_phi) && (Muon_highPtId == 2) && (Muon_tkIsoId >= 1))') \
            .Filter('num_muon_selected == 1', "Just 1 selected muon") \
            .Define(
               "lepton_tuneP_pt",
               "Muon_tunepRelPt[ROOT::VecOps::ArgMax((Muon_pt == lepton_trg_pt) && (Muon_eta == lepton_trg_eta) && (Muon_phi == lepton_trg_phi) && (Muon_highPtId == 2) && (Muon_tkIsoId >= 1))] * lepton_trg_pt"
            ) \
            # .Redefine('lepton_trg_pt', 'lepton_tuneP_pt') \
            .Define('lepton_trg_n_lep', 'Sum((lepton.pt > 55) && (abs(lepton.eta) < 2.4))') \
            #Si ya existe la rama Num_lep_above_15, redefinela, si no, créala
            # .Redefine('Num_lep_above_15', "Sum(Muon_pt > 15 && abs(Muon_eta) < 2.4) + Sum(Electron_pt > 15 && abs(Electron_eta) < 2.4)")
            .Filter('Num_lep_above_15 == 1', "Only one lepton above 15 GeV") \
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
                f"Sum((Jet_pt > 30) &&( abs(Jet_eta) < 2.4) && (Jet_btagDeepFlavB > {TWPbtag2023}) && Jet_passJetIdTightLepVeto)"
            )
            .Define(
                "bjet_pt",
                f"Jet_pt[(Jet_pt > 30) &&( abs(Jet_eta) < 2.4) && (Jet_btagDeepFlavB > {TWPbtag2023}) && Jet_passJetIdTightLepVeto]"
            )
            .Define(
                "bjet_eta",
                f"Jet_eta[(Jet_pt > 30) &&( abs(Jet_eta) < 2.4) && (Jet_btagDeepFlavB > {TWPbtag2023}) && Jet_passJetIdTightLepVeto]"
            )
            .Define(
                "bjet_btag",
                f"Jet_btagDeepFlavB[(Jet_pt > 30) &&( abs(Jet_eta) < 2.4) && (Jet_btagDeepFlavB > {TWPbtag2023}) && Jet_passJetIdTightLepVeto]"
            )
            .Filter("n_bjets >= 1", "At least one b-tagged jet with pt > 30 GeV and |eta| < 2.4 passing tight ID with lepton veto")
            # .Filter(
            #     f"Sum(Jet_pt > 15)>=4",
            #     "4 jets at least"
            # )


            # 3. pt_miss > 50
#             .Filter("pt_miss > 50", "Missing pt selection")
            .Filter("PuppiMET_pt > 50", "Missing pt selection")

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

####################################APARENTLY THIS CAUSES THE ~0.73 OFFSET BETWEEN MC AND DATA###################################################
#             # 5. Subjets dentro del topjet con pt>30 y eta<2.5 
            # .Filter(
            #     f"Sum(subjets.pt[topjets.subjets_in_topjet] > 30 && abs(subjets.eta[topjets.subjets_in_topjet]) < 2.5) >= 3",
            #     "Subjet pt and eta inside topjet"
            # )
#################################################################################################################################################            

            # 6. There must be 3 subjets inside the topjet
            # .Filter("topjets.n_subjets == 3", "Three subjets inside the topjet")
            # 7. The XCone jet from hadronic side must have pt > 350 GeV
            .Filter("topjets.pt > 350")
            # .Filter("fatjets.n_jets > 0", "At least one fatjet from XCone")

#             .Filter("HLT_Mu50", "triggerMuons")
            .Define("pass_trigger",
                    f"(abs(lepton.pdgId[selected_lepton_idx]) == 13 && HLT_Mu50)") # || \
#                      (abs(lepton.pdgId[0]) == 11 && lepton.pt[0] > 120 && (HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Photon175)) || \
#                      (abs(lepton.pdgId[0]) == 11 && lepton.pt[0] < 120 && HLT_Ele30_WPTight_Gsf)")
            .Filter("pass_trigger", "Trigger selection")

            # MET filters
            .Define("pass_MET_filters_2023", "Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_BadPFMuonDzFilter && Flag_hfNoisyHitsFilter && Flag_eeBadScFilter")
            .Filter("pass_MET_filters_2023", "MET filters selection")

        )
        if isMC:
            # Aplicar corrección de b-tagging si es MC
            branches = set(rdf_filtered.GetColumnNames())
            if "btagWeight" in branches:
                rdf_filtered = rdf_filtered.Redefine("btagWeight", "compute_btagWeight(Jet_pt, Jet_eta, Jet_hadronFlavour, Jet_btagDeepFlavB)") \
                                           .Filter("std::isfinite(btagWeight)", "btagWeight is finite") \
                                           .Redefine("eventWeight", "eventWeight * btagWeight")
            else:
                rdf_filtered = rdf_filtered.Define("btagWeight", "compute_btagWeight(Jet_pt, Jet_eta, Jet_hadronFlavour, Jet_btagDeepFlavB)") \
                                           .Filter("std::isfinite(btagWeight)", "btagWeight is finite") \
                                           .Redefine("eventWeight", "eventWeight * btagWeight")

            if "TopPtWeight_NNLOpNLOEW" in branches:
                rdf_filtered = rdf_filtered.Filter("std::isfinite(TopPtWeight_NNLOpNLOEW)", "TopPtWeight_NNLOpNLOEW is finite") \
                                           .Redefine("eventWeight", "eventWeight * TopPtWeight_NNLOpNLOEW")
            else:
                rdf_filtered = rdf_filtered.Define("TopPtWeight_NNLOpNLOEW", "1.")
                

            rdf_filtered = rdf_filtered.Define("muoWeight", "compute_MUOWeight(lepton_trg_pt, lepton_trg_eta)") \
                                       .Filter("std::isfinite(muoWeight)", "muoWeight is finite") \
                                       .Redefine("eventWeight", "eventWeight * muoWeight")

            rdf_filtered = rdf_filtered.Define("PUReweighting", "compute_PUReweight(Pileup_nTrueInt)") \
                                       .Filter("std::isfinite(PUReweighting)", "PUReweighting is finite") \
                                       .Redefine("eventWeight", "eventWeight * PUReweighting")

        filtered_rdfs.append(rdf_filtered)

    return filtered_rdfs


columns_mc = ["eventWeight", "btagWeight", "muoWeight", "PUReweighting", #"event", "run", "luminosityBlock", "lepton.pt", "lepton.eta", "lepton.phi", "lepton.pdgId", "lepton.n_lep", "lepton.dR_to_jet", "lepton.pt_rel_to_jet", "lepton.closest_jet_idx", #Uncomment this when lepton already has dR and ptrel
#               "muon.pt", "muon.eta", "muon.phi", "muon.pdgId", "electron.pt", "electron.eta", "electron.phi", "electron.pdgId",
            #   "Muon_pt", "Muon_eta", "Muon_phi", "Muon_pdgId", "Electron_pt", "Electron_eta", "Electron_phi", "Electron_pdgId", "Num_lep_above_15", "Muon_jetIdx",
              "lepton_trg_pt", "lepton_tuneP_pt", "lepton_trg_eta", "lepton_trg_phi", "lepton_trg_pdgId", "lepton_trg_n_lep", "lepton_trg_dR", "lepton_trg_ptrel",
            #   "PuppiMET_pt", "Jet_pt", "HLT_Mu50", "HLT_Ele30_WPTight_Gsf", "HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Photon175",
#               "electron_trg.pt", "muon.pt",
              "Jet_pt", "Jet_eta", "Jet_btagDeepFlavB", "Jet_hadronFlavour", "Jet_passJetIdTightLepVeto", "Jet_passJetIdTight",
              "Jet_phi", "Jet_rawFactor", "Jet_area","Rho_fixedGridRhoFastjetAll", #"Jet_passJetIdTightLepVeto",#
              "bjet_pt", "bjet_eta", "bjet_btag", "n_bjets", #"pt_miss",
              "PuppiMET_pt",
              "fatjets.n_jets", "fatjets.pt", "fatjets.eta",
              "fatjets.phi", "fatjets.mass",
              "subjets.n_jets", "subjets.pt", "subjets.eta", "subjets.phi", "subjets.mass", "subjets.findLepton",
              "topjets.n_subjets", "topjets.pt", "topjets.eta", "topjets.phi", "topjets.mass",
              "topjets.subjets_in_topjet", "topjets.pt_W", "topjets.eta_W", "topjets.phi_W",
              "topjets.mass_W", "topjets.subjets_in_W", "fastjet_CandsList", "subjet_CandsList", "fatjets.findLepton",
              "pass_detector_selection", ]#"PFCands_pt", "PFCands_eta", "PFCands_phi", "PFCands_pdgId"]
columns_data = [col for col in columns_mc if col not in ["eventWeight", "btagWeight", "muoWeight", "PUReweighting", "Jet_hadronFlavour"]]
# columns_data = columns_data + ["Jet_passJetIdTightLepVeto"]

# DATA: MUON 0
if is_data0:
    print("Creating RDataFrames for Data 0 samples...")
    input_dirs_data0 = "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/Data/2023preBPix/Muon0/Run*"
    root_files_data0 = glob(os.path.join(input_dirs_data0, "**/**/*.root"))#"250818*/**/*.root")) #"**/**/*.root"))
    rdf_data0 = ROOT.RDataFrame("Events", root_files_data0)

    print("Applying event selection to Data 0...")
    rdf_data0 = apply_event_selection([rdf_data0], isMC=False)[0]

    print("Iniciando la conversión de RDataFrames a DataFrames para Data 0...")
    df_data0 = pd.DataFrame(rdf_data0.AsNumpy(columns=columns_data))
    print("Conversión completada para Data 0.")

    print("Guardando DataFrame de Data 0 en formato pickle...")
    df_data0.to_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/data0.pkl")
    print("Archivo guardado: data0.pkl")
    print("Number of events: ", len(df_data0))
    del rdf_data0  # Liberar memoria
else:
    print("Cargando DataFrame de Data 0 desde archivo pickle existente...")
    df_data0 = pd.read_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/data0.pkl")
    print("Carga completada para Data 0.")

# DATA: MUON 1
if is_data1:
    print("Creating RDataFrames for Data 1 samples...")
    input_dirs_data1 = "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/Data/2023preBPix/Muon1/Run*"
    root_files_data1 = glob(os.path.join(input_dirs_data1, "**/**/*.root"))#"250818*/**/*.root")) #"**/**/*.root"))
    rdf_data1 = ROOT.RDataFrame("Events", root_files_data1)

    print("Applying event selection to Data 1...")
    rdf_data1 = apply_event_selection([rdf_data1], isMC=False)[0]

    print("Iniciando la conversión de RDataFrames a DataFrames para Data 1...")
    df_data1 = pd.DataFrame(rdf_data1.AsNumpy(columns=columns_data))
    print("Conversión completada para Data.")

    print("Guardando DataFrame de Data 1 en formato pickle...")
    df_data1.to_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/data1.pkl")
    print("Archivo guardado: data1.pkl")
    print("Number of events: ", len(df_data1))
    del rdf_data1  # Liberar memoria
else:
    print("Cargando DataFrame de Data 1 desde archivo pickle existente...")
    df_data1 = pd.read_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/data1.pkl")
    print("Carga completada para Data 1.")

df_data = pd.concat([df_data0, df_data1], ignore_index=True)
del df_data0
del df_data1

# MC: SINGLE TOP
if is_singletop:
    print("Creating RDataFrames for Single Top samples...")
    rdfs_singletop = [
        create_rdf_with_weights(path, xs, luminosity)#, max_files=5)
        for path, xs in zip(input_dirs_singletop, cross_sections_singletop)
    ]

    print("Applying event selection to Single Top samples...")
    rdfs_singletop = apply_event_selection(rdfs_singletop)

    print("Iniciando la conversión de RDataFrames a DataFrames para SingleTop...")
    df_singletop = pd.concat([
        pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        for rdf in rdfs_singletop
    ], ignore_index=True)
    print("Conversión completada para SingleTop.")
    print("Guardando DataFrame de SingleTop en formato pickle...")
    df_singletop.to_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/singletop.pkl")
    print("Archivo guardado: singletop.pkl")
    print("Number of events: ", len(df_singletop))
    del rdfs_singletop  # Liberar memoria
else:
    print("Cargando DataFrame de Single Top desde archivo pickle existente...")
    df_singletop = pd.read_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/singletop.pkl")
    print("Carga completada para Single Top.")

# MC: W+JETS
if is_wjets:
    print("Creating RDataFrames for W+jets samples...")
    rdfs_wjets = [
        create_rdf_with_weights(path, xs, luminosity)#, max_files=5)
        for path, xs in zip(input_dirs_wjets, cross_sections_wjets)
    ]

    print("Applying event selection to W+jets samples...")
    rdfs_wjets = apply_event_selection(rdfs_wjets)

    print("Iniciando la conversión de RDataFrames a DataFrames para WJets...")
    df_wjets = pd.concat([
        pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        for rdf in rdfs_wjets
    ], ignore_index=True)
    print("Conversión completada para WJets.")
    df_wjets = df_wjets[
        (abs(df_wjets['lepton_trg_pdgId']) == 13) &
    #     (
    #         df_wjets['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
    #         df_wjets['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
    #         == 1
    #     )
        (
            df_wjets.apply( lambda row: sum(pt > 30 and abs(eta) < 2.4 and btag > 0.6553 and passjetid
                                                for pt, eta, btag, passjetid in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'], row['Jet_passJetIdTightLepVeto'])) >= 1, axis=1)
        )
    ]
    print("Guardando DataFrame de WJets en formato pickle...")
    df_wjets.to_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/wjets.pkl")
    print("Archivo guardado: wjets.pkl")
    print("Number of events: ", len(df_wjets))
    del rdfs_wjets  # Liberar memoria
else:
    print("Cargando DataFrame de W+jets desde archivo pickle existente...")
    df_wjets = pd.read_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/wjets.pkl")
    print("Carga completada para W+jets.")

# MC: TTBAR SEMILEPTONIC
if is_ttbar_semi:
    print("Creating RDataFrames for ttbar semileptonic samples...")
    rdfs_ttbar = [
        create_rdf_with_weights(path, xs, luminosity, isTTbar=True)#, max_files=5)
        for path, xs in zip(input_dirs_ttbar, cross_sections_ttbar)
    ]

    print("Applying event selection to ttbar semileptonic samples...")
    rdfs_ttbar= apply_event_selection(rdfs_ttbar)

    print("Iniciando la conversión de RDataFrames a DataFrames para TTBar...")
    df_ttbar = pd.concat([
        pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        for rdf in rdfs_ttbar
    ], ignore_index=True)
    print("Conversión completada para TTBar.")
    print("Guardando DataFrame de TTBar en formato pickle...")
    df_ttbar.to_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/ttbar_semi.pkl")
    print("Archivo guardado: ttbar_semi.pkl")
    print("Number of events: ", len(df_ttbar))
    del rdfs_ttbar  # Liberar memoria
else:
    print("Cargando DataFrame de TTBar semileptonic desde archivo pickle existente...")
    df_ttbar = pd.read_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/ttbar_semi.pkl")
    print("Carga completada para TTBar semileptonic.")

    # # Suponiendo que el peso total se llama "weight" (ajusta si es distinto)
    # df_ttbar["eventWeight"] = df_ttbar["eventWeight"] / 0.73
    # print("Guardando ttbar signal sin extra weight")
    # # Guarda un nuevo pickle sin el factor
    # df_ttbar.to_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/ttbar_semi_beforeBTagMoreStats.pkl")

# MC: TTBAR BACKGROUND
if is_ttbar_bck:
    print("Creating RDataFrames for ttbar background samples...")
    rdfs_ttbar_bck = [
        create_rdf_with_weights(path, xs, luminosity, isTTbar=True)#, max_files=5)
        for path, xs in zip(input_dirs_ttbar_bck, cross_sections_ttbar_bck)
    ]

    print("Applying event selection to ttbar background samples...")
    rdfs_ttbar_bck= apply_event_selection(rdfs_ttbar_bck)

    print("Iniciando la conversión de RDataFrames a DataFrames para TTBar Background...")
    df_ttbar_bck = pd.concat([
        pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        for rdf in rdfs_ttbar_bck
    ], ignore_index=True)
    print("Conversión completada para TTBar Background.")
    print("Guardando DataFrame de TTBar Background en formato pickle...")
    df_ttbar_bck.to_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/ttbar_bck.pkl")
    print("Archivo guardado: ttbar_bck.pkl")
    print("Number of events: ", len(df_ttbar_bck))
    del rdfs_ttbar_bck  # Liberar memoria
else:
    print("Cargando DataFrame de TTBar Background desde archivo pickle existente...")
    df_ttbar_bck = pd.read_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/ttbar_bck.pkl")
    print("Carga completada para TTBar Background.")

    # df_ttbar_bck["eventWeight"] = df_ttbar_bck["eventWeight"] / 0.73
    # print("Guardando ttbar bck sin extra weight")
    # #Guarda un nuevo pickle sin el factor
    # df_ttbar_bck.to_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/ttbar_bck_beforeBTag.pkl")

# MC: QCD MUON
if is_qcdmu:
    print("Creating RDataFrames for QCD muon samples...")
    rdfs_qcdmu = [
        create_rdf_with_weights(path, xs, luminosity)#, max_files=5)
        for path, xs in zip(input_dirs_qcdmu, cross_sections_qcdmu)
    ]

    print("Applying event selection to QCD muon samples...")
    rdfs_qcdmu = apply_event_selection(rdfs_qcdmu)

    print("Iniciando la conversión de RDataFrames a DataFrames para QCDMu...")
    df_qcdmu = pd.concat([
        pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        for rdf in rdfs_qcdmu
    ], ignore_index=True)
    print("Conversión completada para QCDMu.")
    print("Guardando DataFrame de QCDMu en formato pickle...")
    df_qcdmu.to_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/qcdmu.pkl")
    print("Archivo guardado: qcdmu.pkl")
    print("Number of events: ", len(df_qcdmu))
    del rdfs_qcdmu  # Liberar memoria
else:
    print("Cargando DataFrame de QCDMu desde archivo pickle existente...")
    df_qcdmu = pd.read_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/qcdmu.pkl")
    print("Carga completada para QCDMu.")

# MC: DY
if is_dy:
    print("Creating RDataFrames for DY samples...")
    rdfs_dy = [
        create_rdf_with_weights(path, xs, luminosity)#, max_files=5)
        for path, xs in zip(input_dirs_dy, cross_sections_dy)
    ]

    print("Applying event selection to DY samples...")
    rdfs_dy = apply_event_selection(rdfs_dy)

    print("Iniciando la conversión de RDataFrames a DataFrames para DY...")
    df_dy = pd.concat([
        pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        for rdf in rdfs_dy
    ], ignore_index=True)
    print("Conversión completada para DY.")
    print("Guardando DataFrame de DY en formato pickle...")
    df_dy.to_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/dy.pkl")
    print("Archivo guardado: dy.pkl")
    print("Number of events dy: ", len(df_dy))
    del rdfs_dy  # Liberar memoria
else:
    print("Cargando DataFrame de DY desde archivo pickle existente...")
    df_dy = pd.read_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/dy.pkl")
    print("Carga completada para DY.")

# MC: VV
if is_vv:
    print("Creating RDataFrames for VV samples...")
    rdfs_vv = [
        create_rdf_with_weights(path, xs, luminosity)#, max_files=5)
        for path, xs in zip(input_dirs_vv, cross_sections_vv)
    ]

    print("Applying event selection to VV samples...")
    rdfs_vv = apply_event_selection(rdfs_vv)

    print("Iniciando la conversión de RDataFrames a DataFrames para VV...")
    df_vv = pd.concat([
        pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        for rdf in rdfs_vv
    ], ignore_index=True)
    print("Conversión completada para VV.")
    print("Guardando DataFrame de VV en formato pickle...")
    df_vv.to_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/vv.pkl")
    print("Archivo guardado: vv.pkl")
    print("Number of events vv: ", len(df_vv))
    del rdfs_vv  # Liberar memoria
else:
    print("Cargando DataFrame de VV desde archivo pickle existente...")
    df_vv = pd.read_pickle("/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/vv.pkl")
    print("Carga completada para VV.")


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
print("DataFrame de Data convertido a listas.")
df_ttbar = convert_cpp_vectors_to_list(df_ttbar)
print("DataFrame de TTBar convertido a listas.")
df_ttbar_bck = convert_cpp_vectors_to_list(df_ttbar_bck)
print("DataFrame de TTBar Background convertido a listas.")
df_singletop = convert_cpp_vectors_to_list(df_singletop)
print("DataFrame de SingleTop convertido a listas.")
df_wjets = convert_cpp_vectors_to_list(df_wjets)
print("DataFrame de WJets convertido a listas.")
df_qcdmu = convert_cpp_vectors_to_list(df_qcdmu)
print("DataFrame de QCDMu convertido a listas.")
df_dy = convert_cpp_vectors_to_list(df_dy)
print("DataFrame de DY convertido a listas.")
df_vv = convert_cpp_vectors_to_list(df_vv)
print("DataFrame de VV convertido a listas.")

import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.ticker import ScalarFormatter
from matplotlib.gridspec import GridSpec

hep.style.use("CMS")

def save_histograms_to_root(hist_data, hist_mc, bin_edges, labels_mc, output_file):
    """
    Guarda los histogramas en un archivo ROOT.

    Args:
        hist_data (array): Histograma de datos.
        hist_mc (list of arrays): Histogramas de Monte Carlo.
        bin_edges (array): Bordes de los bins.
        labels_mc (list): Etiquetas de los canales MC.
        output_file (str): Nombre del archivo ROOT.
    """
    # Crear el archivo ROOT
    root_file = ROOT.TFile(output_file, "RECREATE")

    # Guardar el histograma de datos
    h_data = ROOT.TH1F("Data", "Data", len(bin_edges) - 1, bin_edges)
    for i, value in enumerate(hist_data):
        h_data.SetBinContent(i + 1, value)
    h_data.Write()

    # Guardar los histogramas de MC
    for i, (hist, label) in enumerate(zip(hist_mc, labels_mc)):
        h_mc = ROOT.TH1F(f"MC_{label}", label, len(bin_edges) - 1, bin_edges)
        for j, value in enumerate(hist):
            h_mc.SetBinContent(j + 1, value)
        h_mc.Write()

    # Cerrar el archivo ROOT
    root_file.Close()
    print(f"Histogramas guardados en: {output_file}")


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
    output_dir="/eos/project/r/rtu-topanalysis/cmunozdi/DataFramesPKL/2023preBPix/AllSelectionCuts_noTightKinHadSubJets/plots",
    nth_element=None  # Nuevo argumento
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
                # Si es iterable, asumimos que puede contener jets
                for val in item:
                    if isinstance(val, (int, float, np.integer, np.floating)):
                        flat_data.append(val)
            except TypeError:
                # No iterable, pero puede ser número suelto
                if isinstance(item, (int, float, np.integer, np.floating)):
                    flat_data.append(item)
                # De lo contrario, se ignora

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
                # jets no iterable → ignoramos
                continue

        return np.array(flat_data), np.array(flat_weights)

    # Verificar que todos los bins sean >= 0
    def verify_bins_non_negative(hist, label):
        if np.any(hist < 0):
            print(f"Advertencia: El histograma '{label}' contiene bins con valores negativos.")
            print(f"Bins negativos: {hist[hist < 0]}")
            # Opcional: Ajustar los valores negativos a 0
            hist[hist < 0] = 0
        return hist
    # Crear el directorio de salida si no existe
    os.makedirs(output_dir, exist_ok=True)
    output_dir_hist = output_dir.replace("plots", "Histograms")
    os.makedirs(output_dir_hist, exist_ok=True)

    # Crear figura y subplots
    fig = plt.figure(figsize=(8, 12))
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.1)
    gs = GridSpec(3, 1, height_ratios=[3, 1, 1], hspace=0.1)
    ax_main = fig.add_subplot(gs[0])
    ax_fraction = fig.add_subplot(gs[1], sharex=ax_main)
    ax_ratio = fig.add_subplot(gs[2], sharex=ax_main)

    # Labels and luminosity
    ax_main.text(
        0, 1.05,
        "Private work (CMS Data)",
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

    # Preparar los datos
    bin_edges = np.linspace(xlim[0], xlim[1], bins + 1) if xlim else bins
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Filtrar el i-ésimo elemento de la rama en los datos
    if (nth_element is not None) and (nth_element > -1):
        df_data_filtered = df_data[branch].apply(lambda x: x[nth_element] if len(x) > nth_element else np.nan).dropna()
    elif (nth_element == -1):
        print(f"Using all elements in branch '{branch}'")
        df_data_filtered = df_to_numpy(df_data, branch)
    else:
        # print("Estoy aquiiiiiiii")
        df_data_filtered = df_data[branch]

    hist_data, _ = np.histogram(df_data_filtered, bins=bin_edges)
    errors_data = np.sqrt(hist_data)

    # Histogramas de MC
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
        # Verificar que los valores del histograma sean >= 0
        # hist = verify_bins_non_negative(hist, label)
        hist_mc.append(hist)
        mc_distributions.append(mc_filtered)
        mc_weights.append(weights_filtered)

    root_output_file = os.path.join(output_dir_hist, f"{title}.root")
    save_histograms_to_root(hist_data, hist_mc, bin_edges, labels_mc, root_output_file)
    # print(f"number of elements in ttbar: {len(dfs_mc[0])}")
    # print(f"number of elements in weights: {len(mc_weights[0])}")
    # print(f"shape of ttbar DataFrame: {dfs_mc[0].shape}")
    # print(f"shape of weights for ttbar: {mc_weights[0].shape}")
    # print("\n")
    # Dibujar los histogramas de MC apilados
    ax_main.hist(
        [bin_edges[:-1]] * len(hist_mc),  # Usar los bordes de los bins para cada histograma
        bins=bin_edges,
        weights=hist_mc,  # Usar los histogramas corregidos
        stacked=True,
        label=labels_mc,
        color=colors_mc,
        alpha=0.7,
    )

    # Configurar el gráfico principal
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
    handles, labels = ax_main.get_legend_handles_labels()
    ax_main.legend(handles[::-1], labels[::-1], frameon=False, fontsize=16) #loc="upper right", bbox_to_anchor=(1,1), borderaxespad=0,
    ax_main.grid(True, linestyle="--", alpha=0.5)
    ax_main.tick_params(axis='x', labelbottom=False)
    # if title:
    #     ax_main.set_title(title)

    # Calcular fracciones para el subplot de fracciones
    total_mc = np.sum(hist_mc, axis=0)
    total_mc[total_mc == 0] = 1  # Evitar divisiones por cero
    fractions = [hist / total_mc for hist in hist_mc]

    # Dibujar las fracciones en el subplot de fracciones (stacked)
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

    # Calcular el ratio Data/MC para el subplot de ratios
    total_mc_sum = np.sum(hist_mc, axis=0)
    # total_mc_sum[total_mc_sum == 0] = -999  # Evitar divisiones por cero
    ratio = np.where(total_mc_sum > 0, (hist_data) / total_mc_sum, np.nan)
    ratio_err = np.where(total_mc_sum > 0, abs(errors_data / total_mc_sum), np.nan)

    # Dibujar el ratio en el subplot de ratios
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

    # Guardar el gráfico como archivo
    output_file = os.path.join(output_dir, f"{title}_JetTIDAndXConeLepIso.png")
    # hep.cms.label("Private", loc=0, rlabel=f"{luminosity/1000:.2f} fb$^{{-1}}$", ax=ax_main, fontsize=22)
    # plt.tight_layout()
    plt.savefig(output_file)
    print(f"Gráfico guardado en: {output_file}")

    # Cerrar el gráfico para liberar memoria
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


print("\n\nApplying tight b-tagging ...")
# Filtrar los DataFrames para muones
df_data_muon = df_data[
    (abs(df_data['lepton_trg_pdgId']) == 13) &
    (df_data['lepton_trg_pt'] > 55) &
    (abs(df_data['lepton_trg_eta']) < 2.4) &
#     (
#         df_data['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
#         df_data['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
#         == 1
#     )
    (
        df_data.apply( lambda row: sum(pt > 30 and abs(eta) < 2.4 and btag > 0.6553 and passjetid
                                              for pt, eta, btag, passjetid in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'], row['Jet_passJetIdTightLepVeto'])) >= 1, axis=1)
    )
 ]
print("Number of events in data muon: ", len(df_data_muon))

df_ttbar_muon = df_ttbar[
    (abs(df_ttbar['lepton_trg_pdgId']) == 13) &
    (df_ttbar['lepton_trg_pt'] > 55) &
    (abs(df_ttbar['lepton_trg_eta']) < 2.4) &
#     (
#         df_ttbar['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
#         df_ttbar['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
#         == 1
#     )
    (
        df_ttbar.apply( lambda row: sum(pt > 30 and abs(eta) < 2.4 and btag > 0.6553 and passjetid
                                              for pt, eta, btag, passjetid in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'], row['Jet_passJetIdTightLepVeto'])) >= 1, axis=1)
    )
 ]
print("Number of events in ttbar muon: ", len(df_ttbar_muon))

df_ttbar_bck_muon = df_ttbar_bck[
    (abs(df_ttbar_bck['lepton_trg_pdgId']) == 13) &
    (df_ttbar_bck['lepton_trg_pt'] > 55) &
    (abs(df_ttbar_bck['lepton_trg_eta']) < 2.4) &
#     (
#         df_ttbar_bck['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
#         df_ttbar_bck['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
#         == 1
#     )
    (
        df_ttbar_bck.apply( lambda row: sum(pt > 30 and abs(eta) < 2.4 and btag > 0.6553 and passjetid
                                              for pt, eta, btag, passjetid in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'], row['Jet_passJetIdTightLepVeto'])) >= 1, axis=1)
    )
 ]
print("Number of events in ttbar background muon: ", len(df_ttbar_bck_muon))

df_singletop_muon = df_singletop[
    (abs(df_singletop['lepton_trg_pdgId']) == 13) &
    (df_singletop['lepton_trg_pt'] > 55) &
    (abs(df_singletop['lepton_trg_eta']) < 2.4) &
#     (
#         df_singletop['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
#         df_singletop['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
#         == 1
#     )
    (
        df_singletop.apply( lambda row: sum(pt > 30 and abs(eta) < 2.4 and btag > 0.6553 and passjetid
                                              for pt, eta, btag, passjetid in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'], row['Jet_passJetIdTightLepVeto'])) >= 1, axis=1)
    )
 ]
print("Number of events in single top muon: ", len(df_singletop_muon))

df_wjets_muon = df_wjets[
    (abs(df_wjets['lepton_trg_pdgId']) == 13) &
    (df_wjets['lepton_trg_pt'] > 55) &
    (abs(df_wjets['lepton_trg_eta']) < 2.4) &
#     (
#         df_wjets['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
#         df_wjets['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
#         == 1
#     )
    (
        df_wjets.apply( lambda row: sum(pt > 30 and abs(eta) < 2.4 and btag > 0.6553 and passjetid
                                              for pt, eta, btag, passjetid in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'], row['Jet_passJetIdTightLepVeto'])) >= 1, axis=1)
    )
 ]
print("Number of events in W+jets muon: ", len(df_wjets_muon))

df_qcdmu_muon = df_qcdmu[
    (abs(df_qcdmu['lepton_trg_pdgId']) == 13) &
    (df_qcdmu['lepton_trg_pt'] > 55) &
    (abs(df_qcdmu['lepton_trg_eta']) < 2.4) &
#     (
#         df_qcdmu['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
#         df_qcdmu['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
#         == 1
#     )
    (
        df_qcdmu.apply( lambda row: sum(pt > 30 and abs(eta) < 2.4 and btag > 0.6553 and passjetid
                                              for pt, eta, btag, passjetid in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'], row['Jet_passJetIdTightLepVeto'])) >= 1, axis=1)
    )
 ]
print("Number of events in QCD muon: ", len(df_qcdmu_muon))

df_dy_muon = df_dy[
    (abs(df_dy['lepton_trg_pdgId']) == 13) &
    (df_dy['lepton_trg_pt'] > 55) &
    (abs(df_dy['lepton_trg_eta']) < 2.4) &
#     (
#         df_others['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
#         df_others['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
#         == 1
#     )
    (
        df_dy.apply( lambda row: sum(pt > 30 and abs(eta) < 2.4 and btag > 0.6553 and passjetid
                                              for pt, eta, btag, passjetid in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'], row['Jet_passJetIdTightLepVeto'])) >= 1, axis=1)
    )
 ]
print("Number of events in DY processes muon: ", len(df_dy_muon))

df_vv_muon = df_vv[
    (abs(df_vv['lepton_trg_pdgId']) == 13) &
    (df_vv['lepton_trg_pt'] > 55) &
    (abs(df_vv['lepton_trg_eta']) < 2.4) &
#     (
#         df_others['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
#         df_others['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
#         == 1
#     )
    (
        df_vv.apply( lambda row: sum(pt > 30 and abs(eta) < 2.4 and btag > 0.6553 and passjetid
                                              for pt, eta, btag, passjetid in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'], row['Jet_passJetIdTightLepVeto'])) >= 1, axis=1)
    )
 ]
print("Number of events in VV processes muon: ", len(df_vv_muon))

print("Plots after tight b-tagging ...")

#LEPTONS
#Numer of leptons
plot_variable( branch="lepton_trg_n_lep", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
              bins=np.arange(-0.5, 3.5, 1), print_percentages=True, #ylim=(0, 3.5*1e3), #logy=True,)
              xlabel='Number of muons', title='tight_lepton_trg_n_lep'
             )

#Lepton pdgId
plot_variable(branch="lepton_trg_pdgId", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
              bins=np.arange(-14.5, 14.5, 1),logy=False, #ylim=(0, 1.85*1e3), #logy=True,)
              xlabel='Lepton pdgID', title='tight_lepton_trg_pdgId'
             )

#Lepton pt
plot_variable( branch="lepton_trg_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
              bins=70, logy=False, xlim=(0, 400), #ylim=(0.75, 1e3),
              xlabel=r'$p_T$ muon [GeV]', title='tight_lepton_trg_pt'
             )

#Lepton TuneP pt              
plot_variable( branch="lepton_tuneP_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
              bins=70, logy=False, xlim=(0, 400), #ylim=(0.75, 1e3),
              xlabel=r'$p_T$ muon TuneP [GeV]', title='tight_lepton_trg_pt'
             )

#Lepton eta
plot_variable(branch="lepton_trg_eta", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
              bins=15, xlim=(-3, 3), #ylim=(0, 0.54*1e3), #logy=True,
              xlabel=r'$\eta$ muon', title='tight_lepton_trg_eta'
             )

#Lepton dR to jet
plot_variable(branch="lepton_trg_dR", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
              bins=20, xlabel=r"$\Delta R$ (muon, closest AK4 jet)", xlim=(0, 4), logy=False, #ylim=(0, 0.8*1e3),
              title='tight_lepton_trg_dR'
             )

#Lepton ptrel to jet
plot_variable(branch="lepton_trg_ptrel", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
              bins=20, xlabel=r"$p_T^{rel}$ (muon, closest AK4 jet) [GeV]", xlim=(0, 200), logy = False, #ylim=(0, 0.35*1e3),
              title='tight_lepton_trg_ptrel'
             )

#XCONE JETS
#Mass top hadronic decay
plot_variable(branch="topjets.mass", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
             labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
             bins=20, xlim=(80, 320), logy=False, #ylim=(0, 0.64*1e3),
             title='tight_topjets_mass', xlabel=r'$m_{top-jet}$ [GeV]'
             )

#Mass W hadronic decay
plot_variable(branch="topjets.mass_W", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
             labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
             bins=20, xlim=(40, 120), logy=False, #ylim=(0, 0.6*1e3),
             title='tight_topjets_Wmass', xlabel=r'$m_{W-jet}$ [GeV]'
             )

#MET
#Missing transverse energy pt
plot_variable(branch="PuppiMET_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
              bins=25, xlim=(0, 600), logy=False, #ylim=(0.75, 1e3),
              xlabel=r'$p_T^{MET}$ [GeV]', title='tight_PuppiMET_pt'
              )

#BJETS
#Numer of bjets
plot_variable( branch="n_bjets", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
              bins=np.arange(-0.5, 6.5, 1), logy=False, #ylim=(0.7, 2.5*1e3), #logy=True,)
              xlabel='Number of b-jets', title='tight_n_bjets'
             )

# #bjets btag discriminator
plot_variable(branch="bjet_btag", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=-1,
              bins=20, xlim=(0., 1.),logy=False, #ylim=(0.95, 4*1e3), #logy=True,)
              xlabel='DeepBTag b-jets', title='tight_bjet_btag'
             )

#bjets pt
plot_variable( branch="bjet_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=-1,
              bins=20, logy=False, xlim=(0, 600), #ylim=(0.75, 4*1e3),
              xlabel=r'$p_T$ b-jets [GeV]', title='tight_bjet_pt'
             )

#bjets eta
plot_variable(branch="bjet_eta", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=-1,
              bins=15, xlim=(-3, 3), #ylim=(0, 0.9*1e3), #logy=True,
              xlabel=r'$\eta$ b-jets', title='tight_bjet_eta'
             )

#LEADDING BJET
#Leading btag discriminator
plot_variable(branch="bjet_btag", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=0,
              bins=20, xlim=(0., 1.),logy=False, #ylim=(0.95, 4.85*1e3), #logy=True,)
              xlabel='DeepBTag leading b-jet', title='tight_bjet_btag_leading'
             )

#Leading bjet pt
plot_variable( branch="bjet_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=0,
              bins=20, logy=False, xlim=(0, 600), #ylim=(0.75, 4*1e3),
              xlabel=r'$p_T$ leading b-jet [GeV]', title='tight_bjet_pt_leading'
             )

#Leading bjet eta
plot_variable(branch="bjet_eta", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=0,
              bins=15, xlim=(-3, 3), #ylim=(0, 6.6*1e2), #logy=True,
              xlabel=r'$\eta$ leading b-jet', title='tight_bjet_eta_leading'
             )

#SUBLEADING BJET
#Subleading btag discriminator
plot_variable(branch="bjet_btag", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=1,
              bins=20, xlim=(0., 1.),logy=False, #ylim=(0.95, 1*1e3), #logy=True,)
              xlabel='DeepBTag subleading b-jet', title='tight_bjet_btag_subleading'
             )

#Subleading bjet pt
plot_variable( branch="bjet_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=1,
              bins=20, logy=False, xlim=(0, 300), #ylim=(0.75, 1*1e3),
              xlabel=r'$p_T$ subleading b-jet [GeV]', title='tight_bjet_pt_subleading'
             )

#Subleading bjet eta
plot_variable(branch="bjet_eta", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=1,
              bins=15, xlim=(-3, 3), #ylim=(0, 2.05*1e2), #logy=True,
              xlabel=r'$\eta$ subleading b-jet', title='tight_bjet_eta_subleading'
             )

#AK4 JETS
#AK4 jets pt
plot_variable( branch="Jet_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=-1,
              bins=20, xlim=(0, 800), logy=False, #ylim=(0.7, 1.1*1e4),
              xlabel=r'$p_T$ AK4 jets [GeV]', title='tight_AK4Jet_pt'
             )

plot_variable( branch="fatjets.pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon,
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1], nth_element=-1,
              bins=20, xlim=(0, 800), logy=False, #ylim=(0.7, 1.1*1e4),
              xlabel=r'$p_T$ Fat-jets [GeV]', title='tight_fatjets_pt'
             )             
