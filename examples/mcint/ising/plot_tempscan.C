
void plot_tempscan(TString rootfilename="tempscan.root",
		   TString savedir="NOSAVE")
{

  
  TFile* f = new TFile(rootfilename.Data());

  THnSparse* hist = (THnSparse*) f->Get("hist");

  
  // -------------------------------------------------------------------------
  cnv_magn_temp = new TCanvas("cnv_magn_temp","cnv_magn_temp",1000,800);

  TH2* h_magn_temp = hist->Projection(3,1);
  h_magn_temp->SetTitle(";k_{B}T / J;magnetisation");
  h_magn_temp->SetStats(0);
  h_magn_temp->Draw("colz");



  // -------------------------------------------------------------------------
  cnv_energy_temp = new TCanvas("cnv_energy_temp","cnv_energy_temp",1000,800);

  TH2* h_enrgy_temp = hist->Projection(2,1);
  h_enrgy_temp->SetTitle(";k_{B}T / J;E");
  h_enrgy_temp->SetStats(0);
  h_enrgy_temp->Draw("colz");
  

  // -------------------------------------------------------------------------
  cnv_magn_energy = new TCanvas("cnv_magn_energy","cnv_energy_temp",1000,800);

  TH2* h_magn_energy = hist->Projection(3,2);
  h_magn_energy->SetTitle(";E;magnetization");
  h_magn_energy->SetStats(0);
  h_magn_energy->Draw("colz");
  

  // -------------------------------------------------------------------------
  if (savedir != "NOSAVE") {
    cnv_magn_temp->SaveAs(savedir + "/tempscan_magn_temp.eps");
    cnv_energy_temp->SaveAs(savedir + "/tempscan_energy_temp.eps");
  }


}
