//Ejecutar con: root -l -b -q plot2D_by_threshold.c

void AddFilesWithWildcards(TChain& chain, const TString& path, const std::vector<TString>& exclude_suffixes) {
    // Ejecutar el comando `find` para buscar archivos con comodines
    TString command = "find " + path + " -name \"*.root\"";
    TString file_list = gSystem->GetFromPipe(command);

    std::istringstream stream(file_list.Data());
    std::string line;
    while (std::getline(stream, line)) {
        TString filename(line);

        // Verificar si el archivo debe ser excluido
        bool exclude = false;
        for (const auto& suffix : exclude_suffixes) {
            if (filename.EndsWith(suffix)) {
                exclude = true;
                break;
            }
        }
        if (!exclude) {
            chain.Add(filename);
        }
    }
}

// Función para agregar archivos excluyendo los que terminan en _events.root, _runs.root, _lumi.root
void AddFilesExcludingRecursive(TChain& chain, const TString& path, const std::vector<TString>& exclude_suffixes) {
    void* dir = gSystem->OpenDirectory(path);
    if (!dir) {
        std::cerr << "Error: No se pudo abrir el directorio " << path << std::endl;
        return;
    }

    const char* file;
    while ((file = gSystem->GetDirEntry(dir))) {
        TString filename(file);

        // Ignorar entradas especiales "." y ".."
        if (filename == "." || filename == "..") {
            continue;
        }

        TString full_path = path + "/" + filename;

        // Verificar si es un subdirectorio
        if (gSystem->AccessPathName(full_path, kFileExists) == 0 && gSystem->OpenDirectory(full_path)) {
            // Llamada recursiva para explorar el subdirectorio
            AddFilesExcludingRecursive(chain, full_path, exclude_suffixes);
        } else if (filename.EndsWith(".root")) {
            // Verificar si el archivo debe ser excluido
            bool exclude = false;
            for (const auto& suffix : exclude_suffixes) {
                if (filename.EndsWith(suffix)) {
                    exclude = true;
                    break;
                }
            }
            if (!exclude) {
                chain.Add(full_path);
            }
        }
    }
    gSystem->FreeDirectory(dir);
}

void plot2D_by_threshold_25() {
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  gROOT->ForceStyle();

  const double ZMIN = 1e-9;
  const double ZMAX = 1e-1;

  const int NX = 320;
  const double XMIN = -8, XMAX = 8;
  const int NY = 320;
  const double YMIN = -200, YMAX = 1400;

  // Rango visual para el PNG
  double VISUAL_XMIN = 0.0, VISUAL_XMAX = 4.0;   // Rango visual para dR
  double VISUAL_YMIN = 0.0, VISUAL_YMAX = 400.0; // Rango visual para ptrel

  int nMaxEvents = -1;

  std::vector<int> thresholds = {25};

  TChain data("Events");

  TChain ttbar("Events");
  TChain ttbar_bck("Events");
  TChain singletop("Events");
  TChain wjets("Events");
  TChain qcd("Events");
  TChain dy("Events");
  TChain vv("Events");



  // Paths para ttbar y qcd
  TString data_path = "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/Data/2023preBPix/Muon*/**/250818*";
  TString ttbar_path = "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/TTbar/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8";
  TString ttbar_bck_path = "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/TTbar/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8";
  TString singletop_path = "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/SingleTop";
  TString wjets_path = "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/WJets";
  TString qcd_path = "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/QCDMu";
  TString dy_path = "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/DY";
  TString vv_path = "/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightID/MC/2023preBPix/VV";

  // Sufijos a excluir
  std::vector<TString> exclude_suffixes = {"_events.root", "_runs.root", "_lumi.root"};

  // Agregar archivos excluyendo los sufijos
  AddFilesWithWildcards(data, data_path, exclude_suffixes);
  AddFilesWithWildcards(ttbar, ttbar_path, exclude_suffixes);
  AddFilesWithWildcards(ttbar_bck, ttbar_bck_path, exclude_suffixes);
  AddFilesWithWildcards(singletop, singletop_path, exclude_suffixes);
  AddFilesWithWildcards(wjets, wjets_path, exclude_suffixes);
  AddFilesWithWildcards(qcd, qcd_path, exclude_suffixes);
  AddFilesWithWildcards(dy, dy_path, exclude_suffixes);
  AddFilesWithWildcards(vv, vv_path, exclude_suffixes);


  std::cout << "[INFO] Número de archivos en data: " << data.GetEntries() << std::endl;
  std::cout << "[INFO] Número de archivos en ttbar: " << ttbar.GetEntries() << std::endl;
  std::cout << "[INFO] Número de archivos en ttbar_bck: " << ttbar_bck.GetEntries() << std::endl;
  std::cout << "[INFO] Número de archivos en singletop: " << singletop.GetEntries() << std::endl;
  std::cout << "[INFO] Número de archivos en wjets: " << wjets.GetEntries() << std::endl;
  std::cout << "[INFO] Número de archivos en qcd: " << qcd.GetEntries() << std::endl;
  std::cout << "[INFO] Número de archivos en dy: " << dy.GetEntries() << std::endl;
  std::cout << "[INFO] Número de archivos en vv: " << vv.GetEntries() << std::endl;


  TFile* fout = new TFile("/eos/home-c/cmunozdi/AnalysisSamples/2dQCDcuts/WithJetTightID_allContrib/histos_25.root", "RECREATE");
  for (int thr : thresholds) {
    TString thr_str = Form("%d", thr);
    std::cout << "\n[INFO] Processing threshold " << thr << std::endl;

    for (auto& [label, chain] : std::vector<std::pair<TString, TChain*>>{
           {"Data", &data},
           {"TTbar", &ttbar},
           {"TTbar_bck", &ttbar_bck},
           {"SingleTop", &singletop},
           {"WJets", &wjets},
           {"QCD", &qcd},
           {"DY", &dy},
           {"VV", &vv},}) {

      TString xvar = "lepton.dR_to_jet";
      TString yvar = "lepton.pt_rel_to_jet";
      TString cut  = Form("lepton.n_lep == 1 && HLT_Mu50");

      TString histname = Form("h2d_%s_thr%s", label.Data(), thr_str.Data());
      TH2F* h2d = new TH2F(histname, "", NX, XMIN, XMAX, NY, YMIN, YMAX);

      TString draw_expr = yvar + ":" + xvar + " >> " + histname;

      std::cout << "[INFO] Drawing for " << label << " at threshold " << thr << std::endl;
      std::cout << "       Branches: " << xvar << " vs " << yvar << std::endl;
      std::cout << "       Cut: " << cut << std::endl;

      Long64_t n_drawn = 0;
      if(nMaxEvents > 0){
        std::cout << "       Max events: " << nMaxEvents << std::endl;
        n_drawn = chain->Draw(draw_expr, cut, "goff", nMaxEvents);
      }else{
        std::cout << "       No limit on events." << std::endl;
        n_drawn = chain->Draw(draw_expr, cut, "goff");
      }
      std::cout << "       Events drawn: " << n_drawn << std::endl;

      // Sustituir ceros por mínimo
      for (int i = 1; i <= h2d->GetNbinsX(); ++i) {
        for (int j = 1; j <= h2d->GetNbinsY(); ++j) {
          if (h2d->GetBinContent(i, j) <= 0)
            h2d->SetBinContent(i, j, ZMIN);
        }
      }

      // Normalización
      double integral = h2d->Integral();
      std::cout << "       Histogram integral before norm: " << integral << std::endl;
      if (integral > 0)
        h2d->Scale(1.0 / integral);


      h2d->GetXaxis()->SetRangeUser(VISUAL_XMIN, VISUAL_XMAX);
      h2d->GetYaxis()->SetRangeUser(VISUAL_YMIN, VISUAL_YMAX);     

      // Dibujar
      TCanvas* c = new TCanvas(Form("c_%s", histname.Data()), histname, 800, 700);
      c->SetLogz();
      h2d->SetMinimum(ZMIN);
      // h2d->SetMaximum(ZMAX);
      h2d->GetXaxis()->SetTitle("#DeltaR (lepton, closest AK4 jet)");
      h2d->GetYaxis()->SetTitle("p_{T}^{rel} (lepton, closest AK4 jet) [GeV]");
      h2d->SetTitle(Form("%s, threshold %s", label.Data(), thr_str.Data()));

      c->SetRightMargin(0.15);  // Espacio para la barra de color COLZ
      c->SetLeftMargin(0.13);   // Para que el eje Y no se corte
      c->SetBottomMargin(0.12); // Más espacio para etiquetas

      h2d->Draw("colz");

      // Guardar
      TString fname = Form("/eos/home-c/cmunozdi/AnalysisSamples/2dQCDcuts/WithJetTightID_allContrib/%s_thr%s.png", label.Data(), thr_str.Data());
      c->SaveAs(fname);
      std::cout << "       Saved plot to " << fname << std::endl;

      fout->cd();
      h2d->Write();

      // =============================
      // Proyecciones 1D
      // =============================
      TH1D* hprojX = h2d->ProjectionX(Form("hprojX_%s_thr%s", label.Data(), thr_str.Data()));
      TH1D* hprojY = h2d->ProjectionY(Form("hprojY_%s_thr%s", label.Data(), thr_str.Data()));

      // Normalización opcional de las proyecciones
      if (hprojX->Integral() > 0) hprojX->Scale(1.0 / hprojX->Integral());
      if (hprojY->Integral() > 0) hprojY->Scale(1.0 / hprojY->Integral());

      // Plot ΔR
      TCanvas* cx = new TCanvas(Form("cx_%s_thr%s", label.Data(), thr_str.Data()), "", 800, 600);
      cx->SetLogy();   // <<< LOG SCALE
      hprojX->SetLineColor(kBlue+1);
      hprojX->SetLineWidth(2);
      hprojX->GetXaxis()->SetTitle("#DeltaR (lepton, closest AK4 jet)");
      hprojX->GetYaxis()->SetTitle("Normalized events");
      hprojX->SetTitle(Form("%s ΔR projection, thr %s", label.Data(), thr_str.Data()));
      hprojX->Draw("hist");
      TString fnameX = Form("/eos/home-c/cmunozdi/AnalysisSamples/2dQCDcuts/WithJetTightID_allContrib/%s_thr%s_projDeltaR.png", label.Data(), thr_str.Data());
      cx->SaveAs(fnameX);

      // Plot pTrel
      TCanvas* cy = new TCanvas(Form("cy_%s_thr%s", label.Data(), thr_str.Data()), "", 800, 600);
      cy->SetLogy();   // <<< LOG SCALE
      hprojY->SetLineColor(kRed+1);
      hprojY->SetLineWidth(2);
      hprojY->GetXaxis()->SetTitle("p_{T}^{rel} (lepton, closest AK4 jet) [GeV]");
      hprojY->GetYaxis()->SetTitle("Normalized events");
      hprojY->SetTitle(Form("%s p_{T}^{rel} projection, thr %s", label.Data(), thr_str.Data()));
      hprojY->Draw("hist");
      TString fnameY = Form("/eos/home-c/cmunozdi/AnalysisSamples/2dQCDcuts/WithJetTightID_allContrib/%s_thr%s_projPtrel.png", label.Data(), thr_str.Data());
      cy->SaveAs(fnameY);

      // Escribir en rootfile también
      hprojX->Write();
      hprojY->Write();      

    }
  }

  fout->Close();

  std::cout << "\n[INFO] Finished all plots." << std::endl;
}
