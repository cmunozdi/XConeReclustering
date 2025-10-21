import ROOT
import os
import numpy as np
import json

import time

# Adding time execution measurement
start_time = time.time()

#Input files TTbar semi directory
input_files = "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/MC/2024/TTbar/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2_BTV_Run3_2024_Comm_MINIAODv6/251007_070711/000*/*root"

#Output directory for the generated files (json and root)
output_dir = "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/output_efficiencies_rdf_topSemi_2024"
os.makedirs(output_dir, exist_ok=True)

print(f"🔍 Processing ROOT files from: {input_files}")
print(f"📂 Output files will be saved in: {output_dir}")

#Defining pt and eta edges
pt_edges = np.array([30, 50, 70, 100, 140, 200, 300, 600, 1000, float('inf')], dtype="float64")
eta_edges = np.array([-2.5, -1.6, -0.8, 0. , 0.8, 1.6, 2.5], dtype="float64")

# Btag WP para 2023 pre-BPix (Tight)
# Btag Tight WP for 2023 pre-BPix: 0.6553 (DeepJet: https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer23/#ak4-b-tagging)
# Btag Tight WP for 2024: 0.4648 (UParTAK4: https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer24/#ak4-b-tagging)
btag_WP = 0.4648

# Create the RDataFrame
rdf = ROOT.RDataFrame("Events", input_files)

# Apply kinematic and WP masks directly in the filters

        #  .Define('Num_lep_above_15', "Sum(Muon_pt > 15 && abs(Muon_eta) < 2.4) + Sum(Electron_pt > 15 && abs(Electron_eta) < 2.4)") \
print("🎭 Applying kinematic and WP masks...")
#rdf = rdf.Define("Jet_btag_pass", f"Jet_btagDeepFlavB > {btag_WP}")  # Evaluate the btag WP directly. For 2024 onwards, use Jet_btagUParTAK4B branch
rdf = rdf.Define("Jet_btag_pass", f"Jet_btagUParTAK4B > {btag_WP}")  # Evaluate the btag WP directly. For 2024 onwards, use Jet_btagUParTAK4B branch
rdf = rdf.Filter("pass_detector_selection", "Detector selection") \
         .Filter('Sum(lepton.pt>55 && abs(lepton.eta)<2.4) == 1') \
         .Define('selected_lepton_idx', 'ROOT::VecOps::ArgMax(lepton.pt > 55 && abs(lepton.eta) < 2.4)') \
         .Filter('abs(lepton.pdgId[selected_lepton_idx]) == 13') \
         .Define('lepton_trg_pt', 'lepton.pt[selected_lepton_idx]') \
         .Define('lepton_trg_eta', 'lepton.eta[selected_lepton_idx]') \
         .Define('lepton_trg_phi', 'lepton.phi[selected_lepton_idx]') \
         .Filter("Sum((Muon_pt == lepton_trg_pt) && (Muon_eta == lepton_trg_eta) && (Muon_phi == lepton_trg_phi) && (Muon_highPtId == 2) && (Muon_tkIsoId >= 1)) == 1", "Muon seleccionado con highPtId==2") \
         .Filter('Num_lep_above_15 == 1', "Only one lepton above 15 GeV") \
         .Filter("PuppiMET_pt > 50", "Missing pt selction") \
         .Define("pass_trigger", 
                    f"(abs(lepton.pdgId[selected_lepton_idx]) == 13 && HLT_Mu50)") \
         .Filter("pass_trigger", "Trigger selection") \
         .Filter(
                f"Sum(subjets.pt[topjets.subjets_in_topjet] > 30 && abs(subjets.eta[topjets.subjets_in_topjet]) < 2.5) >= 3",
                "Subjet pt and eta inside topjet"
                ) \
         .Filter("topjets.n_subjets == 3", "Three subjets inside the topjet") \
         .Filter("topjets.pt > 350") \


# Calculate efficiencies in pt and eta bins for b, c, and light-flavour jets
def calculate_efficiency(rdf, flavour_filter, eta_edges, pt_edges):
    print(f"📊 Calculating efficiency for {flavour_filter}...")
    # denom = rdf.Filter(f"ROOT::VecOps::Any((Jet_pt > 30) && (abs(Jet_eta) < 2.4) && Jet_passJetIdTightLepVeto && {flavour_filter})") \
    #            .Histo2D(
    #                ("denom", "Denominator", len(eta_edges)-1, eta_edges, len(pt_edges)-1, pt_edges),
    #                "Jet_eta", "Jet_pt"
    #            )
    # num = rdf.Filter(f"ROOT::VecOps::Any((Jet_pt > 30) && (abs(Jet_eta) < 2.4) && Jet_passJetIdTightLepVeto && Jet_btag_pass && {flavour_filter})") \
    #           .Histo2D(
    #               ("num", "Numerator", len(eta_edges)-1, eta_edges, len(pt_edges)-1, pt_edges),
    #               "Jet_eta", "Jet_pt"
    #           )
    # eff = num.Clone("eff")
    # eff.Divide(denom.GetValue())
    # print(f"✅ Efficiency calculated for {flavour_filter}.")
    # Filter jets that meet the conditions
    rdf = rdf.Define("selected_jets", f"(Jet_pt > 30) && (abs(Jet_eta) < 2.4) && Jet_passJetIdTightLepVeto && {flavour_filter}")
    rdf = rdf.Define("selected_btag_jets", "selected_jets && Jet_btag_pass")

    # Create histograms directly from the selected jets
    denom = rdf.Filter("Sum(selected_jets) > 0") \
               .Histo2D(
                   ("denom", "Denominator", len(eta_edges)-1, eta_edges, len(pt_edges)-1, pt_edges),
                   "Jet_eta", "Jet_pt",
                   "selected_jets"
               )
    num = rdf.Filter("Sum(selected_btag_jets) > 0") \
              .Histo2D(
                  ("num", "Numerator", len(eta_edges)-1, eta_edges, len(pt_edges)-1, pt_edges),
                  "Jet_eta", "Jet_pt",
                  "selected_btag_jets"
              )

    # Calculate efficiency
    eff = num.Clone("eff")
    eff.Divide(denom.GetValue())
    print(f"✅ Efficiency calculated for {flavour_filter}.")
    return eff

# Calculate efficiencies for each flavour
eff_b = calculate_efficiency(rdf, "Jet_hadronFlavour == 5", eta_edges, pt_edges)
eff_c = calculate_efficiency(rdf, "Jet_hadronFlavour == 4", eta_edges, pt_edges)
eff_lf = calculate_efficiency(rdf, "Jet_hadronFlavour == 0", eta_edges, pt_edges)

# Save to ROOT for validation
print("💾 Saving efficiencies to ROOT...")
outfile = ROOT.TFile(os.path.join(output_dir, "btag_efficiencies_rdf.root"), "RECREATE")
eff_b.Write("eff_b")
eff_c.Write("eff_c")
eff_lf.Write("eff_lf")
outfile.Close()
print("✅ ROOT file saved.")


def save_json_correctionlib_combined(filename, efficiencies, eta_edges, pt_edges):
    print(f"   📄 Saving combined JSON with category in flavour: {filename}")

    # Node category content
    category_content = []
    for eff, flavour, name, description in efficiencies:
        # Extract the contents of the TH2D histogram in 2D format
        content = []
        for i in range(1, len(eta_edges)):  # bins in eta
            # row = []
            for j in range(1, len(pt_edges)):  # bins in pt
                # row.append(eff.GetBinContent(i, j))
                content.append(eff.GetBinContent(i, j))
            # content.append(row)

        # Add the corresponding category for this flavour
        category_content.append({
            "key": int(flavour),  # e.g. 5=b, 4=c, 0=light
            "value": {
                "nodetype": "multibinning",
                "inputs": ["eta", "pt"],
                "edges": [
                    eta_edges.tolist(),
                    pt_edges.tolist()
                ],
                "content": content,
                "flow": "clamp"
            }
        })

    # Main structure
    corrections = [{
        "name": "btag_efficiency",
        "description": "B-tagging efficiency for different jet flavours",
        "version": 1,
        "inputs": [
            {"name": "eta", "type": "real"},
            {"name": "pt", "type": "real"},
            {"name": "flavour", "type": "int"}
        ],
        "output": {"name": "efficiency", "type": "real"},
        "data": {
            "nodetype": "category",
            "input": "flavour",
            "content": category_content
        }
    }]

    # Save the JSON to a file
    combined_json = {
        "schema_version": 2,
        "corrections": corrections
    }
    with open(filename, "w") as f:
        json.dump(combined_json, f, indent=2)

# Save combined JSON
print("💾 Saving efficiencies to a single JSON file...")
efficiencies = [
    (eff_b, 5, "btag_eff_b", "b"),
    (eff_c, 4, "btag_eff_c", "c"),
    (eff_lf, 0, "btag_eff_lf", "light-flavour")
]
save_json_correctionlib_combined(os.path.join(output_dir, "btag_efficiencies_combined.json"), efficiencies, eta_edges, pt_edges)

# Printing time execution
end_time = time.time()
execution_time = end_time - start_time
print(f"⏱️ Execution time: {execution_time:.2f} seconds")# --- IGNORE ---
# In hours, minutes, seconds
hours, rem = divmod(execution_time, 3600)
minutes, seconds = divmod(rem, 60)
print(f"⏱️ Execution time: {int(hours)}h {int(minutes)}m {int(seconds)}s")
# --- IGNORE ---