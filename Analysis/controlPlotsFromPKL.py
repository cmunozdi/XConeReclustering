# . /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc11-opt/setup.sh
import ROOT
import cppyy  # Importar cppyy para acceder a ROOT.VecOps.RVec
import numpy as np
import pandas as pd
import os
from glob import glob
import subprocess
import pickle

luminosity = 18.083517794*1000 #7.229453396*1000 # pb^-1   18.084440726*1000

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

# columns_to_convert = ['Jet_pt', 'Jet_eta', 'Jet_btagDeepFlavB']

print("Loading data from pickle files...")
df_ttbar = pd.read_pickle("/eos/project/r/rtu-topanalysis/DataFramesPKL/2023preBPix/baselineTightBTag/ttbar_semi_tightBTag.pkl")
print("TTbar data loaded successfully.")
print("Number of events in ttbar: ", len(df_ttbar))
df_qcdmu = pd.read_pickle("/eos/project/r/rtu-topanalysis/DataFramesPKL/2023preBPix/baselineTightBTag/qcdmu_tightBTag.pkl")
print("QCD muon data loaded successfully.")
print("Number of events in QCD muon: ", len(df_qcdmu))
df_data0 = pd.read_pickle("/eos/project/r/rtu-topanalysis/DataFramesPKL/2023preBPix/baselineTightBTag/data0_TightBTag.pkl")
print("Data0 loaded successfully.")
df_data1 = pd.read_pickle("/eos/project/r/rtu-topanalysis/DataFramesPKL/2023preBPix/baselineTightBTag/data1_TightBTag.pkl")
print("Data1 loaded successfully.")
df_data = pd.concat([df_data0, df_data1], ignore_index=True)
del df_data0
del df_data1
print("Number of events in data: ", len(df_data))
df_ttbar_bck = pd.read_pickle("/eos/project/r/rtu-topanalysis/DataFramesPKL/2023preBPix/baselineTightBTag/ttbar_bck_tightBTag.pkl")
print("TTbar background data loaded successfully.")
print("Number of events in ttbar background: ", len(df_ttbar_bck))
df_singletop = pd.read_pickle("/eos/project/r/rtu-topanalysis/DataFramesPKL/2023preBPix/baselineTightBTag/singletop_tightBTag.pkl")
print("Single top data loaded successfully.")
print("Number of events in single top: ", len(df_singletop))
df_wjets = pd.read_pickle("/eos/project/r/rtu-topanalysis/DataFramesPKL/2023preBPix/baselineTightBTag/wjets_tightBTag.pkl")
print("W+jets data loaded successfully.")
print("Number of events in W+jets: ", len(df_wjets))
df_dy = pd.read_pickle("/eos/project/r/rtu-topanalysis/DataFramesPKL/2023preBPix/baselineTightBTag/dy_tightBTag.pkl")
print("DY processes data loaded successfully.")
print("Number of events in DY processes: ", len(df_dy))
df_vv = pd.read_pickle("/eos/project/r/rtu-topanalysis/DataFramesPKL/2023preBPix/baselineTightBTag/vv_tightBTag.pkl")
print("VV processes data loaded successfully.")
print("Number of events in VV processes: ", len(df_vv))
print("Data loaded successfully.")

df_data = convert_cpp_vectors_to_list(df_data)
print("Data converted to list successfully.")
df_ttbar = convert_cpp_vectors_to_list(df_ttbar)
print("TTbar data converted to list successfully.")
df_ttbar_bck = convert_cpp_vectors_to_list(df_ttbar_bck)
print("TTbar background data converted to list successfully.")
df_singletop = convert_cpp_vectors_to_list(df_singletop)
print("Single top data converted to list successfully.")
df_wjets = convert_cpp_vectors_to_list(df_wjets)
print("W+jets data converted to list successfully.")
df_qcdmu = convert_cpp_vectors_to_list(df_qcdmu)
print("QCD muon data converted to list successfully.")
df_dy = convert_cpp_vectors_to_list(df_dy)
print("DY processes data converted to list successfully.")
df_vv = convert_cpp_vectors_to_list(df_vv)
print("VV processes data converted to list successfully.")

# print(type(df_data['Jet_pt'].iloc[0]))
# print(type(df_data['Jet_eta'].iloc[0]))
# print(type(df_data['Jet_btagDeepFlavB'].iloc[0]))


import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import pandas as pd
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
    output_dir="/eos/project/r/rtu-topanalysis/DataFramesPKL/2023preBPix/baselineTightBTag/plots",
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

    # Crear figura y subplots
    fig = plt.figure(figsize=(8, 12))
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.1)
    gs = GridSpec(3, 1, height_ratios=[3, 1, 1], hspace=0.1)
    ax_main = fig.add_subplot(gs[0])
    ax_fraction = fig.add_subplot(gs[1], sharex=ax_main)
    ax_ratio = fig.add_subplot(gs[2], sharex=ax_main)

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
        hist = verify_bins_non_negative(hist, label)
        hist_mc.append(hist)
        mc_distributions.append(mc_filtered)
        mc_weights.append(weights_filtered)

    root_output_file = os.path.join(output_dir, f"{title}_btagSFs.root")
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
        ax_main.get_yaxis().get_offset_text().set_y(0.5)
        ax_main.ticklabel_format(axis="y", style="sci", scilimits=(3, 3))
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
    ratio = np.where(total_mc_sum > 0, (hist_data - total_mc_sum) / total_mc_sum, np.nan)
    ratio_err = np.where(total_mc_sum > 0, abs(errors_data / total_mc_sum), np.nan)

    # Dibujar el ratio en el subplot de ratios
    ax_ratio.errorbar(
        bin_centers,
        ratio,
        yerr=ratio_err,
        fmt="o",
        color="black",
    )
    ax_ratio.axhline(0, color="red", linestyle="--", linewidth=1)
    ax_ratio.set_ylabel(r'$\frac{\mathrm{Data} - \mathrm{MC}}{\mathrm{MC}}$')
    ax_ratio.set_xlabel(xlabel or branch)
    ax_ratio.set_ylim(-1, 1)
    ax_ratio.grid(True, linestyle="--", alpha=0.5)

    # Guardar el gráfico como archivo
    output_file = os.path.join(output_dir, f"{title}_btagSFs.png")
    hep.cms.label("Preliminary", loc=0, rlabel=f"{luminosity/1000:.2f} fb$^{{-1}}$", ax=ax_main, fontsize=22)
    plt.tight_layout()
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


# Filtrar los DataFrames para muones
df_data_muon = df_data#[
#     (abs(df_data['lepton_trg_pdgId']) == 13) &
#     (
#         df_data['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
#         df_data['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
#         == 1
#     )
    # (
    #     df_data.apply( lambda row: any(pt > 30 and abs(eta) < 2.4 and btag > 0.6553
    #                                           for pt, eta, btag in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'])))   
    # )
# ]
print("Number of events in data muon: ", len(df_data_muon))

df_ttbar_muon = df_ttbar#[
#     (abs(df_ttbar['lepton_trg_pdgId']) == 13) &
#     (
#         df_ttbar['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
#         df_ttbar['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
#         == 1
#     )
    # (
    #     df_ttbar.apply( lambda row: any(pt > 30 and abs(eta) < 2.4 and btag > 0.6553
    #                                           for pt, eta, btag in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'])))   
    # )
# ]
print("Number of events in ttbar muon: ", len(df_ttbar_muon))

df_ttbar_bck_muon = df_ttbar_bck#[
#     (abs(df_ttbar_bck['lepton_trg_pdgId']) == 13) &
#     (
#         df_ttbar_bck['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
#         df_ttbar_bck['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
#         == 1
#     )
    # (
    #     df_ttbar_bck.apply( lambda row: any(pt > 30 and abs(eta) < 2.4 and btag > 0.6553
    #                                           for pt, eta, btag in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'])))   
    # )
# ]
print("Number of events in ttbar background muon: ", len(df_ttbar_bck_muon))

df_singletop_muon = df_singletop#[
#     (abs(df_singletop['lepton_trg_pdgId']) == 13) &
#     (
#         df_singletop['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
#         df_singletop['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
#         == 1
#     )
    # (
    #     df_singletop.apply( lambda row: any(pt > 30 and abs(eta) < 2.4 and btag > 0.6553
    #                                           for pt, eta, btag in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'])))   
    # )
# ]
print("Number of events in single top muon: ", len(df_singletop_muon))

df_wjets_muon = df_wjets#[
#     (abs(df_wjets['lepton_trg_pdgId']) == 13) &
#     (
#         df_wjets['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
#         df_wjets['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
#         == 1
#     )
    # (
    #     df_wjets.apply( lambda row: any(pt > 30 and abs(eta) < 2.4 and btag > 0.6553
    #                                           for pt, eta, btag in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'])))   
    # )
# ]
print("Number of events in W+jets muon: ", len(df_wjets_muon))

df_qcdmu_muon = df_qcdmu#[
#     (abs(df_qcdmu['lepton_trg_pdgId']) == 13) &
#     (
#         df_qcdmu['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
#         df_qcdmu['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
#         == 1
#     )
    # (
    #     df_qcdmu.apply( lambda row: any(pt > 30 and abs(eta) < 2.4 and btag > 0.6553
    #                                           for pt, eta, btag in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'])))   
    # )
# ]
print("Number of events in QCD muon: ", len(df_qcdmu_muon))

df_dy_muon = df_dy#[
#     (abs(df_others['lepton_trg_pdgId']) == 13) &
#     (
#         df_others['Muon_pt'].apply(lambda x: sum(pt > 15 for pt in x)) +
#         df_others['Electron_pt'].apply(lambda x: sum(pt > 15 for pt in x))
#         == 1
#     )
    # (
    #     df_others.apply( lambda row: any(pt > 30 and abs(eta) < 2.4 and btag > 0.6553
    #                                           for pt, eta, btag in zip(row['Jet_pt'], row['Jet_eta'], row['Jet_btagDeepFlavB'])))   
    # )
# ]
print("Number of events in DY processes muon: ", len(df_dy_muon))
df_vv_muon = df_vv
print("Number of events in VV processes muon: ", len(df_vv_muon))


# ##########################################LEPTONS######################################
# #Numer of leptons
# plot_variable( branch="lepton_trg_n_lep", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_others_muon][::-1], df_data=df_data_muon, 
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY&VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue"][::-1],
#               bins=np.arange(-0.5, 3.5, 1), print_percentages=True, ylim=(0, 280*1e3), #logy=True,)
#               xlabel='Number of muons'
#              )
              
# #Lepton pdgId
# plot_variable(branch="lepton_trg_pdgId", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_others_muon][::-1], df_data=df_data_muon, 
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY&VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue"][::-1],
#               bins=np.arange(-14.5, 14.5, 1),logy=False, ylim=(0, 140*1e3), #logy=True,)
#               xlabel='Lepton pdgID'
#              )
              
# #Lepton pt            
# plot_variable( branch="lepton_trg_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_others_muon][::-1], df_data=df_data_muon, 
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY&VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue"][::-1],
#               bins=45, logy=True, xlim=(0, 600), ylim=(0.7, 5*1e5), 
#               xlabel=r'$p_T$ muon [GeV]',
#              )

# #Lepton eta
# plot_variable(branch="lepton_trg_eta", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_others_muon][::-1], df_data=df_data_muon, 
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY&VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue"][::-1],
#               bins=15, xlim=(-3, 3), ylim=(0, 38*1e3), #logy=True,
#               xlabel=r'$\eta$ muon',
#              )

# #Lepton dR to jet
# plot_variable(branch="lepton_trg_dR", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_others_muon][::-1], df_data=df_data_muon, 
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY&VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue"][::-1],
#               bins=40, xlabel=r"$\Delta R$ (muon, closest AK4 jet)", xlim=(0, 4), logy=False, ylim=(0, 16*1e3)#, logy=True,
#              )

# #Lepton ptrel to jet
# plot_variable(branch="lepton_trg_ptrel", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_others_muon][::-1], df_data=df_data_muon, 
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY&VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue"][::-1],
#               bins=20, xlabel=r"$p_T^{rel}$ (muon, closest AK4 jet) [GeV]", xlim=(0, 200), logy = False, ylim=(0, 35*1e3)#, logy=True,
#              )

# #########################################AK4JETS########################################
# # All ak4-jets pt            
# plot_variable( branch="Jet_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_others_muon][::-1], df_data=df_data_muon, 
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY&VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue"][::-1],
#               bins=50, xlim=(0, 1000), logy=True, ylim=(0.7, 3*1e6),
#               xlabel=r'$p_T$ jet [GeV]', print_percentages=True, nth_element=-1
#              )

#####################################B-JETS######################################

#Number of b-jets
plot_variable(branch="n_bjets", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon, 
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
              bins=np.arange(-0.5, 7.5, 1), logy=True, ylim=(0.7, 5*1e5),# xlim=(-3, 3), ylim=(0.7, 2*1e2), logy=True,
             )

#ALL bjet
#pt
plot_variable(branch="bjet_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon, 
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
              bins=50, xlim=(0,1000), xlabel="b-jets $p_T$", logy=True, nth_element=-1, ylim=(0.7, 5*1e5)
             )

#btag
plot_variable(branch="bjet_btag", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon, 
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
              bins=50, xlim=(0,1), xlabel="b-jets btag", logy=True, nth_element=-1, ylim=(0.7, 5*1e5)
             )

#eta
plot_variable(branch="bjet_eta", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_dy_muon, df_vv_muon][::-1], df_data=df_data_muon, 
              labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY", "VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue", "indigo"][::-1],
              bins=20, xlim=(-3,3), xlabel="b-jets $\eta$", nth_element=-1, ylim=(0, 40*1e3)#logy=True, ylim=(0.7, 5*1e2)
             )

# #####################################LEADING B-JET######################################
# #Leading bjet
# #pt
# plot_variable(branch="bjet_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_others_muon][::-1], df_data=df_data_muon, 
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY&VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue"][::-1],
#               title="leadingBjet_pt",
#               bins=50, xlim=(0,1000), nth_element=0, xlabel="Leading b-jet $p_T$", logy=True, ylim=(0.7, 5*1e5)
#              )

# #btag
# plot_variable(branch="bjet_btag", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_others_muon][::-1], df_data=df_data_muon, 
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY&VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue"][::-1],
#               title="leadingBjet_btag",
#               bins=50, xlim=(0,1), nth_element=0, xlabel="Leading b-jet btag", logy=True, ylim=(0.7, 5*1e5)
#              )

# #eta
# plot_variable(branch="bjet_eta", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_others_muon][::-1], df_data=df_data_muon, 
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY&VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue"][::-1],
#                 title="leadingBjet_eta",
#               bins=20, nth_element=0, xlim=(-3,3), xlabel="Leading b-jet $\eta$", ylim=(0, 27*1e3)#logy=True, ylim=(0.7, 5*1e2)
#              )

# #####################################SUBLEADING B-JET######################################
# #Subleading bjet
# #pt
# plot_variable(branch="bjet_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_others_muon][::-1], df_data=df_data_muon, 
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY&VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue"][::-1],
#                 title="subleadingBjet_pt",
#               bins=50, xlim=(0,1000), nth_element=1, xlabel="Subleading b-jet $p_T$", logy=True, ylim=(0.7, 5*1e5)
#              )

# #btag
# plot_variable(branch="bjet_btag", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_others_muon][::-1], df_data=df_data_muon, 
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY&VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue"][::-1],
#                 title="subleadingBjet_btag",
#               bins=50, xlim=(0,1), nth_element=1, xlabel="Subleading b-jet btag", logy=True, ylim=(0.7, 5*1e5)
#              )

# #eta
# plot_variable(branch="bjet_eta", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_others_muon][::-1], df_data=df_data_muon, 
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY&VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue"][::-1],
#                 title="subleadingBjet_eta",
#               bins=20, nth_element=1, xlim=(-3,3), xlabel="Subleading b-jet $\eta$", ylim=(0, 27*1e3)#logy=True, ylim=(0.7, 5*1e2)
#              )

# ###########################################MET pt ########################################
# #MET pt
# plot_variable(branch="PuppiMET_pt", dfs_mc=[df_ttbar_muon, df_ttbar_bck_muon, df_singletop_muon, df_wjets_muon, df_qcdmu_muon, df_others_muon][::-1], df_data=df_data_muon, 
#               labels_mc=[r'$t\overline{t}\rightarrow l\nu2q$', r"$t\overline{t}\rightarrow others$",  "SingleTop", "W+jets", "QCD", "DY&VV"][::-1], colors_mc=["red", "tomato", "gold", "lime", "deepskyblue", "blue"][::-1],
#               bins=50, xlim=(0,1000), logy=True, ylim=(0.7, 5*1e5)
#              )


del df_data
del df_ttbar
del df_ttbar_bck
del df_singletop
del df_wjets
del df_qcdmu
del df_dy
del df_vv