import argparse
import os
import ROOT
import subprocess


# Verify if FASTJET_CONTRIB_BASE is set
fjc_base = os.environ.get("FASTJET_CONTRIB_BASE")
if not fjc_base:
    raise RuntimeError("FASTJET_CONTRIB_BASE is not set. Ensure the environment is configured correctly.")

# Agregar la ruta de los headers al intérprete de ROOT
ROOT.gSystem.AddIncludePath(f'-I{fjc_base}/include')




# Input argument parsing
parser = argparse.ArgumentParser(description="Process PFNANO with XCone algo for MC or Data.")
parser.add_argument("--input", required=True, help="Input file or directory")
parser.add_argument("--output", required=True, help="Output file")
parser.add_argument("--isMC", action="store_true", help="Are we processing MC? (buy default is Data)")

args = parser.parse_args()

# Verify if the input is a file or a directory
if os.path.isdir(args.input):
    input_files = [os.path.join(args.input, f) for f in os.listdir(args.input) if f.endswith(".root")]
elif os.path.isfile(args.input):
    input_files = [args.input]
else:
    raise ValueError(f"Input is not a valid file or direcroty: {args.input}")

output_file_events = args.output.replace(".root", "_events.root")
isMC = args.isMC

# Crear un nombre de archivo para rdf_runs basado en args.output
output_file_runs = args.output.replace(".root", "_runs.root")
output_file_lumi = args.output.replace(".root", "_lumi.root")

# ROOT.gSystem.Load("$LCIO/lib/libRivet.so")
# ROOT.gSystem.Load("$LCIO/lib/libHepMC.so")
# ROOT.gSystem.Load("$LCIO/lib/libfastjet.so")

ROOT.gInterpreter.Declare('#include "selection_helpers_BoostedTopQuark.h"') ##RUN CRAB
# Llamar a la función para inicializar `cset`
ROOT.gInterpreter.ProcessLine(f"initializeCorrectionSet({str(isMC).lower()});")
# verbosity = ROOT.Experimental.RLogScopedVerbosity(ROOT.Detail.RDF.RDFLogChannel(), ROOT.Experimental.ELogLevel.kInfo)

ROOT.ROOT.EnableImplicitMT()
# Creates the RDataFrame
rdf_runs = ROOT.RDataFrame("Runs", input_files)
rdf_runs.Snapshot('Runs', output_file_runs)
rdf_lumi = ROOT.RDataFrame("LuminosityBlocks", input_files)
rdf_lumi.Snapshot('LuminosityBlocks', output_file_lumi)
rdf = ROOT.RDataFrame("Events", input_files)

# Define Jet_passTightID
rdf = rdf.Define("Jet_passJetIdTight", "pass_JetIdTightLepVeto(Jet_eta, Jet_neHEF, Jet_neEmEF, Jet_chMultiplicity, Jet_neMultiplicity, \
                 Jet_chHEF, Jet_muEF, Jet_chEmEF, false)") \
         .Define("Jet_passJetIdTightLepVeto", "pass_JetIdTightLepVeto(Jet_eta, Jet_neHEF, Jet_neEmEF, Jet_chMultiplicity, Jet_neMultiplicity, \
                 Jet_chHEF, Jet_muEF, Jet_chEmEF, true)") \
         .Define("Jet_jetId", "(Jet_passJetIdTight * 2) + (Jet_passJetIdTightLepVeto * 4)")

# # Define PFCands_pt_puppiWeighted and PFCands_mass_puppiWeighted
rdf = rdf.Redefine("PFCands_pt", "PFCands_pt * PFCands_puppiWeight") \
         .Redefine("PFCands_mass", "PFCands_mass * PFCands_puppiWeight")

# Lepton selection: select in the kinematic region of interest. Only events with exactly one lepton are kept
rdf = rdf.Define('lepton', 'pickLepton(Muon_pt, Muon_eta, Muon_phi, Muon_pdgId, Electron_pt, Electron_eta, Electron_phi, Electron_pdgId)') \
         .Define('one_lepton', 'lepton.n_lep == 1')

# Apply trigger selection and lepton ID for those events with one lepton
# According to https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIIISummary#2023 the trigger for Photon should be HLT_Photon200 instead of HLT_Photon175
rdf = rdf.Define('pass_trigger', 'one_lepton && ((abs(lepton.pdgId[0]) == 13 && HLT_Mu50) || (abs(lepton.pdgId[0]) == 11 && ((lepton.pt[0] > 120 && (HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Photon200)) || (lepton.pt[0] > 55 && lepton.pt[0] <= 120 && HLT_Ele30_WPTight_Gsf))))')

# Apply cut to select events with good PVs
rdf = rdf.Define('good_PV', 'PV_npvsGood > 0')

# Apply cut MET filters
rdf = rdf.Define("pass_MET_filters", "Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_BadPFMuonDzFilter && Flag_hfNoisyHitsFilter && Flag_eeBadScFilter")

# Apply cut over PuppiMET_pt
rdf = rdf.Define('pass_MET_cut', 'PuppiMET_pt > 35.0')

# Add to the lepton struct the information related to the closest jet, and define one good lepton variable for further selections
rdf = rdf.Redefine('lepton', 'addLeptonJetInfo(lepton, Muon_jetIdx, Electron_jetIdx, PFCands_pt, PFCands_eta, PFCands_phi, PFCands_pdgId, Jet_pt, Jet_eta, Jet_phi, Jet_mass, Jet_jetId, Jet_rawFactor, Jet_area, Rho_fixedGridRhoFastjetAll, run)') \
         .Define('one_good_lepton', 'one_lepton && lepton.n_lep_after2Dcuts == 1')

# Define potential b-jet from AK4 jets
rdf = rdf.Define('potential_bjet_num', f"Sum((Jet_pt>22.5)&&(abs(Jet_eta)<3)&&(Jet_jetId>5))") \
         .Define('at_least_one_potential_bjet', 'potential_bjet_num > 0')
# CHANGE BTAGGING WORKING POINT AS NEEDED DEPENDING ON THE YEAR
rdf = rdf.Define('bjet_num', f"Sum((Jet_pt>22.5)&&(abs(Jet_eta)<3)&&(Jet_jetId>5)&&(Jet_btagPNetB>0.1917))") \
         .Define('at_least_one_bjet', 'bjet_num > 0')

# XCone jet clustering for PFCands
rdf = rdf.Define('jet_XConefromPFCands','buildXConeJets(lepton, PFCands_pt, PFCands_eta, PFCands_phi, PFCands_mass, PFCands_pdgId)') \
         .Define('n_fatjets', 'jet_XConefromPFCands.fatjets.n_jets')
# To reduce a little bit more the size of the output files, we can request two extra conditions: having exactly 3 subjets on the hadronic side, and having a hadtopjet with a pt avobe 325 GeV:
rdf = rdf.Define('ThreeHadSubjets', 'jet_XConefromPFCands.topjets.n_subjets == 3') \
         .Define('HadTopJet_pt_above_325', 'jet_XConefromPFCands.topjets.pt > 325')



# # XCone for PFCands: pass detector selection flag: #Not final cuts, because some observables need to be corrected before cutting definitelly. At least one b-tagged jet is possible for datasets that are not going to be used for b-tagging efficiencies calculation (e.g. Data + MC samples for background estimation: W+jets, Single top, QCD, VV, DY)
rdf = rdf.Define('pass_detector_selection', 'one_good_lepton && pass_trigger && good_PV && pass_MET_filters && pass_MET_cut && at_least_one_potential_bjet && n_fatjets > 0 && ThreeHadSubjets && HadTopJet_pt_above_325 && at_least_one_bjet')
#&& at_least_one_bjet --> Not to use on ttbar MC samples for b-tagging efficiency calculation

if isMC:
    # Only execute this line if we are processing MC
    rdf = rdf.Define('GenCands_jetIdx', 'calculateGenCandsJetIdx(GenJetCands_genCandsIdx, GenJetCands_jetIdx, nGenCands, nGenJet)') \
             .Define('gen_lepton', 'GenLepton(GenCands_pt, GenCands_eta, GenCands_phi, GenCands_pdgId, GenCands_jetIdx, GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass)') \
             .Define('one_good_gen_lepton', 'gen_lepton.n_lep == 1 && gen_lepton.n_lep_after2Dcuts == 1') \
             .Define('potential_gen_bjet_num', f"Sum((GenJet_pt>22.5)&&(abs(GenJet_eta)<3))") \
             .Define('at_least_one_potential_gen_bjet', 'potential_gen_bjet_num > 0') \
             .Define('jet_XConefromGenCands', 'buildXConeJets(gen_lepton, GenCands_pt, GenCands_eta, GenCands_phi, GenCands_mass, GenCands_pdgId)') \
             .Define('n_gen_fatjets', 'jet_XConefromGenCands.fatjets.n_jets') \
             .Define('ThreeHadSubjets_gen', 'jet_XConefromGenCands.topjets.n_subjets == 3') \
             .Define('HadTopJet_pt_above_325_gen', 'jet_XConefromGenCands.topjets.pt > 325') \
             .Define('pass_particle_selection', 'one_good_gen_lepton && at_least_one_potential_gen_bjet && n_gen_fatjets > 0 && ThreeHadSubjets_gen && HadTopJet_pt_above_325_gen') \
             .Redefine("GenCands_jetIdx", "std::vector<Short_t>(GenCands_jetIdx.begin(), GenCands_jetIdx.end())")
else:
    # If not MC, we don't need to define gen_* variables
    rdf = rdf.Define('pass_particle_selection', 'false')

# Filter the events with at least one flag set to true
rdf = rdf.Filter('pass_detector_selection || pass_particle_selection') 

# Save the columns needed (if running this on crab, maybe better to save all the columns, since selecting the interesting events will notably reduce the size of the file)
# columns = ['event', 'lepton', 'n_jets', 'Jet_btagDeepFlavB', 'pt_miss']
# columns += ['jet_XConefromPFCands']
# columns += ['PFCands_pt', 'PFCands_phi', 'PFCands_eta', 'PFCands_pdgId']

# if isMC:
#     # Only add these columns if we are processing MC
#     columns += ['jet_XConefromGenCands']
#     columns += ['GenCands_pt', 'GenCands_phi', 'GenCands_eta', 'GenCands_pdgId']


# # Verify that columns is a list of strings
# if not isinstance(columns, list):
#     raise TypeError("columns must be a list of strings")
# for col in columns:
#     if not isinstance(col, str):
#         raise TypeError(f"Each element in columns must be a string. Found: {type(col)}")


# opts = ROOT.RDF.RSnapshotOptions()
# opts.fMode = "RECREATE"
# opts.fOverwriteIfExists = True

# Create a snapshot with the selected columns
# Check if rdf has any events before snapshotting
num_events = rdf.Count().GetValue()
print(f"Number of events after filtering: {num_events}")

if num_events > 0:
    rdf.Snapshot('Events', output_file_events)
else:
    # If no events, create a minimal empty TTree with a single dummy branch
    print("Warning: No events passed the selection. Creating minimal empty Events tree.")
    out_file = ROOT.TFile(output_file_events, "RECREATE")
    tree = ROOT.TTree("Events", "Empty Events tree")
    
    # Create a single dummy branch (e.g., a counter variable)
    empty_var = ROOT.std.vector('float')()
    tree.Branch("empty", empty_var)
    
    # Write empty tree (no entries, but structure exists)
    out_file.Write()
    out_file.Close()
    print(f"Created minimal empty Events tree in {output_file_events}")
# output_file.WriteObject(rdf, "Events")
# output_file.WriteObject(rdf_runs, "Runs")
# Copiar el árbol "Runs" desde el primer input file
# f_in = ROOT.TFile.Open(input_files[0])
# t_runs = f_in.Get("Runs")
# f_out = ROOT.TFile.Open(output_file, "UPDATE")
# if t_runs:
#     f_out.cd()
#     t_runs.CloneTree(-1, "fast").Write("Runs")
# f_out.Close()
# f_in.Close()

# rdf_runs.Snapshot('Runs', output_file, "", opts)

ROOT.ROOT.DisableImplicitMT()

print(f"Processing completed. Output file: {args.output}")

# Merge the output files if there are multiple input files
output_files = [output_file_events, output_file_runs, output_file_lumi]

combined_output_file = args.output

hadd_command = ["hadd", "-f", combined_output_file] + output_files

try: 
    subprocess.run(hadd_command, check=True)
    print(f"Successfully merged output files into: {combined_output_file}")
    subprocess.run(["rm"] + output_files, check=True)
except subprocess.CalledProcessError as e:
    print(f"Error merging output files: {e}")
    

