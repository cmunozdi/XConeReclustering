# . /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc11-opt/setup.sh
import ROOT
import numpy as np
import pandas as pd
import os
from glob import glob
#open as a pandas dataframe
# cmssw_base = os.environ['CMSSW_BASE']
# ROOT.gInterpreter.Declare(f'#include "{cmssw_base}/src/XConeReclustering/selection_helpers_BoostedTopQuark.h"')
# ROOT.gSystem.Load(f"{cmssw_base}/src/XConeReclustering/selection_helpers_BoostedTopQuark.so")

#Bool variables to know over which sample should we run the code
is_ttbar_semi = True
is_ttbar_bck = True
is_singletop = True
is_wjets = True
is_qcdmu = True
is_other = True
is_data0 = False
is_data1 = False

ROOT.ROOT.EnableImplicitMT()

def create_rdf_with_weights(input_dir, cross_section, luminosity, max_files=None):
    # Elegir el patrón según si es QCDMu o no
#     pattern = "**/250622*/**/*.root"  if "QCDMu" in input_dir else  "**/250625*/**/*.root" if "TTbar" in input_dir else "**/250625*/**/*.root"
    pattern = "**/250*/**/*.root"
    all_root_files = glob(os.path.join(input_dir, pattern), recursive=True)

    # Filtrar archivos que no terminan en coor.root, lumi.root, o runs.root
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
        
        # Verificar que tenga TTree "Events"
        if f.Get("Events"):
            valid_files.append(file)
        
        # Contar eventos generados si hay TTree "Runs"
        runs_tree = f.Get("Runs")
        if runs_tree:
            for entry in runs_tree:
                total_gen_events_sumW += entry.genEventSumw

        f.Close()
        
        # Si max_files está definido y ya tenemos suficientes archivos válidos, salimos
        if max_files is not None and len(valid_files) >= max_files:
            break

    if not valid_files:
        print(f"Advertencia: No se encontraron archivos válidos con TTree 'Events' en el directorio {input_dir}.")
        return None  # O devolver un RDataFrame vacío si es necesario

    if total_gen_events_sumW == 0:
        print(f"Advertencia: No se encontraron eventos generados en {input_dir}.")
        return None  # O devolver un RDataFrame vacío si es necesario

    event_weight = (luminosity * cross_section) / total_gen_events_sumW

    # Crear el RDataFrame solo con los archivos válidos
    rdf = ROOT.RDataFrame("Events", valid_files)
    rdf = rdf.Define("eventWeight", f"{event_weight}*genWeight")
    return rdf


luminosity = 18.083517794*1000 #7.229453396*1000 # pb^-1  18.084440726

cross_sections_singletop = [87.200000, 1.477177, 4.663095, 19.970880, 19.302840, 145.000000, 2.360095, 4.663095, 19.970880, 19.302840] #[145.0, 87.2, 2.3600952, 1.4771772, 19.30280, 19.30280] 
cross_sections_wjets = [368.200000, 421.900000, 25.600000, 54.770000, 0.878500, 3.119000, 4427.000000, 1598.000000, 0.105300, 0.526100] #[368.200000, 421.900000, 25.600000, 54.770000, 0.878500, 3.119000, 4427.000000, 1598.000000, 0.105300, 0.526100] #[55760, 9529, 3532]#[2994.4278, 959.481, 284.7492, 134.52282, 17948.322, 55740*3*0.1086, 9558*3*0.1086, 3602*3*0.1086, 68360*3*0.1086]
cross_sections_ttbar = [923.6*0.4392] #, 923.6/3*0.4392, 923.6/3*0.4392] #762.1                                            68360
cross_sections_ttbar_bck = [923.6*0.1061, 923.6*0.4544]
cross_sections_qcdmu = [1.323000, 23240.000000, 7763.000000, 699.100000, 1458000.000000, 68.240000, 404400.000000, 21.370000, 3.913000, 95900.000000] #[3.913, 95900]#[7768, 695.2, 67.32, 21.4, 3.929, 1.301] #2964000.000000, 2688000.000000, 
cross_sections_dy = [20950, 141500, 6688] #[6744.000000000000, 21.650000000000, 0.001111000000, 3.058000000000, 0.000059490000, 0.000001558000, 0.269100000000, 2219.000000000000, 0.000000035190, 0.019150000000] #[141500, 20950, 6688, 0.01915, 6744.000000000000, 21.650000000000, 0.001111000000, 3.058000000000, 0.000059490000, 0.000001558000, 0.269100000000, 2219.000000000000, 0.000000035190, 0.019150000000, 80.23, 29.1, 12.75] 
cross_sections_vv = [80.23, 29.1, 12.75]

if is_singletop:
    input_dirs_singletop = [
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/SingleTop/TBbarQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/SingleTop/TbarBQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/SingleTop/TBbartoLplusNuBbar-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/SingleTop/TbarBtoLminusNuB-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/SingleTop/TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/SingleTop/TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",

        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop/TbarBQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop/TbarBtoLminusNuB-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop/TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop/TbarWplusto4Q_TuneCP5_13p6TeV_powheg-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop/TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop/TBbarQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop/TBbartoLplusNuBbar-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop/TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop/TWminusto4Q_TuneCP5_13p6TeV_powheg-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop/TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",


    ]

if is_wjets:
    input_dirs_wjets = [
        # "/eos/user/c/cmunozdi/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-4Jets_1J_TuneCP5_13p6TeV_madgraphMLM-pythia8",
        # "/eos/user/c/cmunozdi/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-4Jets_2J_TuneCP5_13p6TeV_madgraphMLM-pythia8",
        # "/eos/user/c/cmunozdi/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-4Jets_3J_TuneCP5_13p6TeV_madgraphMLM-pythia8",
        # "/eos/user/c/cmunozdi/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-4Jets_4J_TuneCP5_13p6TeV_madgraphMLM-pythia8",
        # "/eos/user/c/cmunozdi/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-4Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8",
        # "/eos/cms/store/group/phys_top/RTU_topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_0J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/cms/store/group/phys_top/RTU_topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/cms/store/group/phys_top/RTU_topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/cms/store/group/phys_top/RTU_topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-200to400_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-400to600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-400to600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-40to100_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-40to100_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8", 

        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-200to400_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-400to600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-400to600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-40to100_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-40to100_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/WJets/WtoLNu-2Jets_PTLNu-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",


   
    ]

if is_ttbar_semi:
    input_dirs_ttbar = [
        # "/eos/user/c/cmunozdi/AnalysisSamples/MC/2023preBPix/TTbar/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/TTbar/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8_moreStats",

        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/TTbar/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",

    ]

if is_ttbar_bck:
    input_dirs_ttbar_bck = [
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/TTbar/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/TTbar/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8"

        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/TTbar/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/TTbar/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8",

    ]

if is_qcdmu:
    input_dirs_qcdmu = [
    #     "/eos/cms/store/group/phys_top/RTU_topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/QCDMu/QCD_PT-80to120_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
    #     "/eos/cms/store/group/phys_top/RTU_topanalysis/cmunozdi/AnalysisSamples/MC/2023preBPix/QCDMu/QCD_PT-120to170_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/QCDMu/QCD_PT-170to300_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/QCDMu/QCD_PT-300to470_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/QCDMu/QCD_PT-470to600_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/QCDMu/QCD_PT-600to800_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/QCDMu/QCD_PT-800to1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/QCDMu/QCD_PT-1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",

        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-120to170_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-15to20_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-170to300_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-20to30_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-300to470_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-30to50_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-470to600_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-50to80_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-600to800_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-800to1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu/QCD_PT-80to120_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",


    ]

if is_other:
    input_dirs_dy = [
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/DY/DYto2L-2Jets_MLL-4to10_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/DY/DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/DY/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu-2Jets_MLL-105To160_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu_MLL-10to50_TuneCP5_13p6TeV_powheg-pythia8_moreStats",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu_MLL-120to200_TuneCP5_13p6TeV_powheg-pythia8_moreStats",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu_MLL-1500to2500_TuneCP5_13p6TeV_powheg-pythia8_moreStats",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu_MLL-200to400_TuneCP5_13p6TeV_powheg-pythia8_moreStats",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu_MLL-2500to4000_TuneCP5_13p6TeV_powheg-pythia8_moreStats",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu_MLL-4000to6000_TuneCP5_13p6TeV_powheg-pythia8_moreStats",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu_MLL-400to800_TuneCP5_13p6TeV_powheg-pythia8_moreStats",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu_MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8_moreStats",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu_MLL-6000_TuneCP5_13p6TeV_powheg-pythia8_moreStats",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/DY/DYto2Mu_MLL-800to1500_TuneCP5_13p6TeV_powheg-pythia8_moreStats",

        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/DY/DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/DY/DYto2L-2Jets_MLL-4to10_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/DY/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",


    ]
    input_dirs_vv = [
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/VV/WW_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/VV/WZ_TuneCP5_13p6TeV_pythia8",
        # "/eos/project/r/rtu-topanalysis/AnalysisSamples/MC/2023preBPix/VV/ZZ_TuneCP5_13p6TeV_pythia8",

        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/VV/WW_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/VV/WZ_TuneCP5_13p6TeV_pythia8",
        "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/VV/ZZ_TuneCP5_13p6TeV_pythia8",
    ]

print("Iniciando la conversión de RDataFrames a DataFrames.")

if is_singletop:
    rdfs_singletop = [
        create_rdf_with_weights(path, xs, luminosity)#, max_files=1)
        for path, xs in zip(input_dirs_singletop, cross_sections_singletop)
    ]

if is_wjets:
    rdfs_wjets = [
        create_rdf_with_weights(path, xs, luminosity)#, max_files=1)
        for path, xs in zip(input_dirs_wjets, cross_sections_wjets)
    ]

if is_ttbar_semi:
    rdfs_ttbar = [
        create_rdf_with_weights(path, xs, luminosity)#, max_files=1)
        for path, xs in zip(input_dirs_ttbar, cross_sections_ttbar)
    ]

if is_ttbar_bck:
    rdfs_ttbar_bck = [
        create_rdf_with_weights(path, xs, luminosity)#, max_files=1)
        for path, xs in zip(input_dirs_ttbar_bck, cross_sections_ttbar_bck)
    ]

if is_qcdmu:
    rdfs_qcdmu = [
        create_rdf_with_weights(path, xs, luminosity)#, max_files=1)
        for path, xs in zip(input_dirs_qcdmu, cross_sections_qcdmu)
    ]

if is_other:
    rdfs_dy = [
        create_rdf_with_weights(path, xs, luminosity)#, max_files=1)
        for path, xs in zip(input_dirs_dy, cross_sections_dy)
    ]
    rdfs_vv = [
        create_rdf_with_weights(path, xs, luminosity)#, max_files=1)
        for path, xs in zip(input_dirs_vv, cross_sections_vv)
    ]

if is_data0:
    input_dirs_data0 = "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/Data/2023preBPix/Muon0/**"
    root_files_data0 = glob(os.path.join(input_dirs_data0, "250818*/**/*.root"))
    rdf_data0 = ROOT.RDataFrame("Events", root_files_data0)

if is_data1:
    input_dirs_data1 = "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/Data/2023preBPix/Muon1/**"
    root_files_data1 = glob(os.path.join(input_dirs_data1, "250818*/**/*.root"))
    rdf_data1 = ROOT.RDataFrame("Events", root_files_data1)

# Imprime el numero de eventos de rdf_data para comprobar:
# print(f"Number of events in rdf_data: {rdf_data0.Count().GetValue()} + {rdf_data1.Count().GetValue()}")

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

def apply_event_selection(rdfs, TWPbtag2023=0.6553):
    filtered_rdfs = []

    for rdf in rdfs:
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
            .Filter('Sum(lepton.pt>55 && abs(lepton.eta)<2.4) == 1')
# #             .Filter("lepton.n_lep>=1")
#             .Define('lepton_trg_pt', 'lepton.pt[0]')
#             .Define('lepton_trg_eta', 'lepton.eta[0]')
#             .Define('lepton_trg_phi', 'lepton.phi[0]')
#             .Define('lepton_trg_pdgId', 'lepton.pdgId[0]')
#             .Define('lepton_trg_dR', 'lepton.dR_to_jet[0]')
#             .Define('lepton_trg_ptrel', 'lepton.pt_rel_to_jet[0]')
#             .Define('lepton_trg_n_lep', 'lepton.n_lep')
            .Define('selected_lepton_idx', 'ROOT::VecOps::ArgMax(lepton.pt > 55 && abs(lepton.eta) < 2.4)') \
            .Filter('abs(lepton.pdgId[selected_lepton_idx]) == 13') \
            .Define('lepton_trg_pt', 'lepton.pt[selected_lepton_idx]') \
            .Define('lepton_trg_eta', 'lepton.eta[selected_lepton_idx]') \
            .Define('lepton_trg_phi', 'lepton.phi[selected_lepton_idx]') \
            .Define('lepton_trg_pdgId', 'lepton.pdgId[selected_lepton_idx]') \
            .Define('lepton_trg_dR', 'lepton.dR_to_jet[selected_lepton_idx]') \
            .Define('lepton_trg_ptrel', 'lepton.pt_rel_to_jet[selected_lepton_idx]') \
            .Define('lepton_trg_n_lep', 'Sum((lepton.pt > 55) && (abs(lepton.eta) < 2.4))') \
            .Define('Num_lep_above_15', "Sum(Muon_pt > 15 && abs(Muon_eta) < 2.4) + Sum(Electron_pt > 15 && abs(Electron_eta) < 2.4)")
            .Filter('Num_lep_above_15 == 1', "Only one lepton above 15 GeV") \
            # .Define('lepton_muon_idx',
            #    'match_lepton_to_muon(lepton_trg_pt, lepton_trg_eta, lepton_trg_phi, Muon_pt, Muon_eta, Muon_phi)') \
            # .Filter('lepton_muon_idx > -1 && Muon_tkIsoId[lepton_muon_idx] >= 1',
            #    'Muon_tkIsoId == 1 if lepton matches a muon')
            # .Filter('lepton_trg_n_lep == 1')

            # 2. Al menos un jet con las condiciones
            # .Filter(
            #     f"Sum(Jet_pt > 30 && abs(Jet_eta) < 2.4 && Jet_btagDeepFlavB > {TWPbtag2023}) >= 1",
            #     "Jet selection"
            # )
            .Define(
                "n_bjets",
                f"Sum(Jet_pt > 30 && abs(Jet_eta) < 2.4 && Jet_btagDeepFlavB > {TWPbtag2023} && Jet_passJetIdTightLepVeto)"
            )
            .Define(
                "bjet_pt",
                f"Jet_pt[Jet_pt > 30 && abs(Jet_eta) < 2.4 && Jet_btagDeepFlavB > {TWPbtag2023} && Jet_passJetIdTightLepVeto]"
            )
            .Define(
                "bjet_eta",
                f"Jet_eta[Jet_pt > 30 && abs(Jet_eta) < 2.4 && Jet_btagDeepFlavB > {TWPbtag2023} && Jet_passJetIdTightLepVeto]"
            )
            .Define(
                "bjet_btag",
                f"Jet_btagDeepFlavB[Jet_pt > 30 && abs(Jet_eta) < 2.4 && Jet_btagDeepFlavB > {TWPbtag2023} && Jet_passJetIdTightLepVeto]"
            )
            # .Filter("n_bjets >= 1", "At least one b-jet with pt > 30, |eta| < 2.4 and btagDeepFlavB > TWPbtag2023") \
            # .Filter(
            #     f"Sum(Jet_pt > 50)>=2",
            #     "2 jets with at least 50 GeV"
            # )
            

            # 3. pt_miss > 50
#             .Filter("pt_miss > 50", "Missing pt selection")
            .Filter("PuppiMET_pt > 50", "Missing pt selction")

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
#             .Filter(
#                 f"Sum(subjets.pt[topjets.subjets_in_topjet] > 30 && abs(subjets.eta[topjets.subjets_in_topjet]) < 2.5) >= 1",
#                 "Subjet pt and eta inside topjet"
#             )
            
#             .Filter("HLT_Mu50", "triggerMuons")
            .Define("pass_trigger", 
                    f"(abs(lepton.pdgId[selected_lepton_idx]) == 13 && HLT_Mu50)") # || \
#                      (abs(lepton.pdgId[0]) == 11 && lepton.pt[0] > 120 && (HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Photon175)) || \
#                      (abs(lepton.pdgId[0]) == 11 && lepton.pt[0] < 120 && HLT_Ele30_WPTight_Gsf)")
            .Filter("pass_trigger", "Trigger selection")

        )

        filtered_rdfs.append(rdf_filtered)

    return filtered_rdfs

if is_data0:
    rdf_data0 = apply_event_selection([rdf_data0])[0]
    # print("Number of events in rdf_data0 after selection: ", rdf_data0.Count().GetValue())
    # rdf_data0 = rdf_data0.Define('lepton_muon_idx',
    #                              'match_lepton_to_muon(lepton_trg_pt, lepton_trg_eta, lepton_trg_phi, Muon_pt, Muon_eta, Muon_phi)') \
    #                      .Filter('lepton_muon_idx > -1 && Muon_tkIsoId[lepton_muon_idx] >= 1',
    #                              'Muon_tkIsoId == 1 if lepton matches a muon')
    # print("Number of events in rdf_data0 after lepton selection: ", rdf_data0.Count().GetValue())

if is_data1:
    rdf_data1 = apply_event_selection([rdf_data1])[0]
    # print("Number of events in rdf_data1 after selection: ", rdf_data1.Count().GetValue())
    # rdf_data1 = rdf_data1.Define('lepton_muon_idx',
    #                              'match_lepton_to_muon(lepton_trg_pt, lepton_trg_eta, lepton_trg_phi, Muon_pt, Muon_eta, Muon_phi)') \
    #                      .Filter('lepton_muon_idx > -1 && Muon_tkIsoId[lepton_muon_idx] >= 1',
    #                              'Muon_tkIsoId == 1 if lepton matches a muon')
    # print("Number of events in rdf_data1 after lepton selection: ", rdf_data1.Count().GetValue())

if is_ttbar_semi:
    rdfs_ttbar= apply_event_selection(rdfs_ttbar)
    # print("Number of events in rdfs_ttbar after selection: ", sum(rdf.Count().GetValue() for rdf in rdfs_ttbar))
    # rdfs_ttbar = (rdf.Define('lepton_muon_idx',
    #                          'match_lepton_to_muon(lepton_trg_pt, lepton_trg_eta, lepton_trg_phi, Muon_pt, Muon_eta, Muon_phi)') \
    #                  .Filter('lepton_muon_idx > -1 && Muon_tkIsoId[lepton_muon_idx] >= 1',
    #                          'Muon_tkIsoId == 1 if lepton matches a muon') for rdf in rdfs_ttbar)
    # print("Number of events in rdf_ttbar after lepton selection: ", sum(rdf.Count().GetValue() for rdf in rdfs_ttbar))

if is_ttbar_bck:
    rdfs_ttbar_bck= apply_event_selection(rdfs_ttbar_bck)
    # print("Number of events in rdfs_ttbar_bck after selection: ", sum(rdf.Count().GetValue() for rdf in rdfs_ttbar_bck))
    # rdfs_ttbar_bck = (rdf.Define('lepton_muon_idx',
    #                          'match_lepton_to_muon(lepton_trg_pt, lepton_trg_eta, lepton_trg_phi, Muon_pt, Muon_eta, Muon_phi)') \
    #                  .Filter('lepton_muon_idx > -1 && Muon_tkIsoId[lepton_muon_idx] >= 1',
    #                          'Muon_tkIsoId == 1 if lepton matches a muon') for rdf in rdfs_ttbar_bck)
    # print("Number of events in rdf_ttbar_bck after lepton selection: ", sum(rdf.Count().GetValue() for rdf in rdfs_ttbar_bck))

if is_singletop:
    rdfs_singletop = apply_event_selection(rdfs_singletop)
    # print("Number of events in rdfs_singletop after selection: ", sum(rdf.Count().GetValue() for rdf in rdfs_singletop))
    # rdfs_singletop = (rdf.Define('lepton_muon_idx',
    #                          'match_lepton_to_muon(lepton_trg_pt, lepton_trg_eta, lepton_trg_phi, Muon_pt, Muon_eta, Muon_phi)') \
    #                  .Filter('lepton_muon_idx > -1 && Muon_tkIsoId[lepton_muon_idx] >= 1',
    #                          'Muon_tkIsoId == 1 if lepton matches a muon') for rdf in rdfs_singletop)
    # print("Number of events in rdf_singletop after lepton selection: ", sum(rdf.Count().GetValue() for rdf in rdfs_singletop))

if is_wjets:
    rdfs_wjets = apply_event_selection(rdfs_wjets)
    # print("Number of events in rdfs_wjets after selection: ", sum(rdf.Count().GetValue() for rdf in rdfs_wjets))
    # rdfs_wjets = (rdf.Define('lepton_muon_idx',
    #                          'match_lepton_to_muon(lepton_trg_pt, lepton_trg_eta, lepton_trg_phi, Muon_pt, Muon_eta, Muon_phi)') \
    #                  .Filter('lepton_muon_idx > -1 && Muon_tkIsoId[lepton_muon_idx] >= 1',
    #                          'Muon_tkIsoId == 1 if lepton matches a muon') for rdf in rdfs_wjets)
    # print("Number of events in rdf_wjets after lepton selection: ", sum(rdf.Count().GetValue() for rdf in rdfs_wjets))

if is_qcdmu:
    rdfs_qcdmu = apply_event_selection(rdfs_qcdmu)
    # print("Number of events in rdfs_qcdmu after selection: ", sum(rdf.Count().GetValue() for rdf in rdfs_qcdmu))
    # rdfs_qcdmu = (rdf.Define('lepton_muon_idx',
    #                          'match_lepton_to_muon(lepton_trg_pt, lepton_trg_eta, lepton_trg_phi, Muon_pt, Muon_eta, Muon_phi)') \
    #                  .Filter('lepton_muon_idx > -1 && Muon_tkIsoId[lepton_muon_idx] >= 1',
    #                          'Muon_tkIsoId == 1 if lepton matches a muon') for rdf in rdfs_qcdmu)
    # print("Number of events in rdf_qcdmu after lepton selection: ", sum(rdf.Count().GetValue() for rdf in rdfs_qcdmu))

if is_other:
    rdfs_dy = apply_event_selection(rdfs_dy)
    # print("Number of events in rdfs_others after selection: ", sum(rdf.Count().GetValue() for rdf in rdfs_others))
    # rdfs_others = (rdf.Define('lepton_muon_idx',
    #                          'match_lepton_to_muon(lepton_trg_pt, lepton_trg_eta, lepton_trg_phi, Muon_pt, Muon_eta, Muon_phi)') \
    #                  .Filter('lepton_muon_idx > -1 && Muon_tkIsoId[lepton_muon_idx] >= 1',
    #                          'Muon_tkIsoId == 1 if lepton matches a muon') for rdf in rdfs_others)
    # print("Number of events in rdf_others after lepton selection: ", sum(rdf.Count().GetValue() for rdf in rdfs_others))
    rdfs_vv = apply_event_selection(rdfs_vv)



columns_mc = ["eventWeight", "lepton.pt", "lepton.eta", "lepton.phi", "lepton.pdgId", "lepton.n_lep", "lepton.dR_to_jet", "lepton.pt_rel_to_jet", #Uncomment this when lepton already has dR and ptrel
#               "muon.pt", "muon.eta", "muon.phi", "muon.pdgId", "electron.pt", "electron.eta", "electron.phi", "electron.pdgId",
              "Muon_pt", "Muon_eta", "Muon_phi", "Muon_pdgId", "Electron_pt", "Electron_eta", "Electron_phi", "Electron_pdgId", "Num_lep_above_15",
              "lepton_trg_pt", "lepton_trg_eta", "lepton_trg_phi", "lepton_trg_pdgId", "lepton_trg_n_lep", "lepton_trg_dR", "lepton_trg_ptrel",
              "PuppiMET_pt", "Jet_pt", "HLT_Mu50", "HLT_Ele30_WPTight_Gsf", "HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Photon175",
#               "electron_trg.pt", "muon.pt",
              "Jet_pt", "Jet_eta", "Jet_btagDeepFlavB", "Jet_passJetIdTightLepVeto", "Jet_hadronFlavour",
              "bjet_pt", "bjet_eta", "bjet_btag", "n_bjets", #"pt_miss", 
            #   "fatjets.n_jets", "fatjets.pt", "fatjets.eta",
            #   "fatjets.phi", "fatjets.mass", "subjets.n_jets", "subjets.pt", "subjets.eta", "subjets.phi", "subjets.mass", 
            #   "topjets.n_subjets", "topjets.pt", "topjets.eta", "topjets.phi", "topjets.mass",
            #   "topjets.subjets_in_topjet", "topjets.pt_W", "topjets.eta_W", "topjets.phi_W", 
            #   "topjets.mass_W", "topjets.subjets_in_W", "fastjet_CandsList", "subjet_CandsList", "fatjets.findLepton", 
              "pass_detector_selection", ]#"PFCands_pt", "PFCands_eta", "PFCands_phi", "PFCands_pdgId"]
columns_data = [col for col in columns_mc if col not in ["eventWeight", "Jet_hadronFlavour"]]

if is_singletop:
    print("Iniciando la conversión de RDataFrames a DataFrames para SingleTop...")
    df_singletop = pd.concat([
        pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        for rdf in rdfs_singletop
    ], ignore_index=True)
    print("Conversión completada para SingleTop.")
    print("Guardando DataFrame de SingleTop en formato pickle...")
    df_singletop.to_pickle("/eos/project/r/rtu-topanalysis/DataFramesPKL/2023preBPix/baselineTightBTag/singletop_tightBTag.pkl")
    print("Archivo guardado: singletop_tightBTag.pkl")
    print("Number of events: ", len(df_singletop))
    # del df_singletop  # Liberar memoria
    del rdfs_singletop  # Liberar memoria

if is_wjets:
    print("Iniciando la conversión de RDataFrames a DataFrames para WJets...")
    df_wjets = pd.concat([
        pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        for rdf in rdfs_wjets
    ], ignore_index=True)
    print("Conversión completada para WJets.")
    print("Guardando DataFrame de WJets en formato pickle...")
    df_wjets.to_pickle("/eos/project/r/rtu-topanalysis/DataFramesPKL/2023preBPix/baselineTightBTag/wjets_tightBTag.pkl")
    print("Archivo guardado: wjets_tightBTag.pkl")
    print("Number of events: ", len(df_wjets))
    # del df_wjets  # Liberar memoria
    del rdfs_wjets  # Liberar memoria

if is_ttbar_semi:
    print("Iniciando la conversión de RDataFrames a DataFrames para TTBar...")
    df_ttbar = pd.concat([
        pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        for rdf in rdfs_ttbar
    ], ignore_index=True)
    print("Conversión completada para TTBar.")
    print("Guardando DataFrame de TTBar en formato pickle...")
    df_ttbar.to_pickle("/eos/project/r/rtu-topanalysis/DataFramesPKL/2023preBPix/baselineTightBTag/ttbar_semi_tightBTag.pkl")
    print("Archivo guardado: ttbar_semi_tightBTag.pkl")
    print("Number of events: ", len(df_ttbar))
    # del df_ttbar  # Liberar memoria
    del rdfs_ttbar  # Liberar memoria

if is_ttbar_bck:
    print("Iniciando la conversión de RDataFrames a DataFrames para TTBar Background...")
    df_ttbar_bck = pd.concat([
        pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        for rdf in rdfs_ttbar_bck
    ], ignore_index=True)
    print("Conversión completada para TTBar Background.")
    print("Guardando DataFrame de TTBar Background en formato pickle...")
    df_ttbar_bck.to_pickle("/eos/project/r/rtu-topanalysis/DataFramesPKL/2023preBPix/baselineTightBTag/ttbar_bck_tightBTag.pkl")
    print("Archivo guardado: ttbar_bck_tightBTag.pkl")
    print("Number of events: ", len(df_ttbar_bck))
    # del df_ttbar_bck  # Liberar memoria
    del rdfs_ttbar_bck  # Liberar memoria

if is_qcdmu:
    print("Iniciando la conversión de RDataFrames a DataFrames para QCDMu...")
    df_qcdmu = pd.concat([
        pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        for rdf in rdfs_qcdmu
    ], ignore_index=True)
    print("Conversión completada para QCDMu.")
    print("Guardando DataFrame de QCDMu en formato pickle...")
    df_qcdmu.to_pickle("/eos/project/r/rtu-topanalysis/DataFramesPKL/2023preBPix/baselineTightBTag/qcdmu_tightBTag.pkl")
    print("Archivo guardado: qcdmu_tightBTag.pkl")
    print("Number of events: ", len(df_qcdmu))
    # del df_qcdmu  # Liberar memoria
    del rdfs_qcdmu  # Liberar memoria

if is_other:
    print("Iniciando la conversión de RDataFrames a DataFrames para DY...")
    df_dy = pd.concat([
        pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        for rdf in rdfs_dy
    ], ignore_index=True)
    print("Conversión completada para DY.")
    print("Guardando DataFrame de DY en formato pickle...")
    df_dy.to_pickle("/eos/project/r/rtu-topanalysis/DataFramesPKL/2023preBPix/baselineTightBTag/dy_tightBTag.pkl")
    print("Archivo guardado: dy_tightBTag.pkl")
    print("Number of events: ", len(df_dy))
    # del df_others  # Liberar memoria
    del rdfs_dy  # Liberar memoria

    print("Iniciando la conversión de RDataFrames a DataFrames para VV...")
    df_vv = pd.concat([
        pd.DataFrame(rdf.AsNumpy(columns=columns_mc))
        for rdf in rdfs_vv
    ], ignore_index=True)
    print("Conversión completada para VV.")
    print("Guardando DataFrame de VV en formato pickle...")
    df_vv.to_pickle("/eos/project/r/rtu-topanalysis/DataFramesPKL/2023preBPix/baselineTightBTag/vv_tightBTag.pkl")
    print("Archivo guardado: vv_tightBTag.pkl")
    print("Number of events: ", len(df_vv))
    # del df_others  # Liberar memoria
    del rdfs_vv  # Liberar memoria    

if is_data0:
    print("Iniciando la conversión de RDataFrames a DataFrames para Data 0...")
    df_data0 = pd.DataFrame(rdf_data0.AsNumpy(columns=columns_data))
    print("Conversión completada para Data 0.")
    
    print("Guardando DataFrame de Data 0 en formato pickle...")
    df_data0.to_pickle("/eos/project/r/rtu-topanalysis/DataFramesPKL/2023preBPix/baselineTightBTag/data0_tightBTag.pkl")
    print("Archivo guardado: data0_tightBTag.pkl")
    print("Number of events: ", len(df_data0))
    # del df_data0  # Liberar memoria
    del rdf_data0  # Liberar memoria

if is_data1:
    print("Iniciando la conversión de RDataFrames a DataFrames para Data 1...")
    df_data1 = pd.DataFrame(rdf_data1.AsNumpy(columns=columns_data))
    print("Conversión completada para Data.")

    print("Guardando DataFrame de Data en formato pickle...")
    df_data1.to_pickle("/eos/project/r/rtu-topanalysis/DataFramesPKL/2023preBPix/baselineTightBTag/data1_tightBTag.pkl")
    print("Archivo guardado: data_tightBTag.pkl")
    print("Number of events: ", len(df_data1))
    # del df_data1  # Liberar memoria
    del rdf_data1  # Liberar memoria