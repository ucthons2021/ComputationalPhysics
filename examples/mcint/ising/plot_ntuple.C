
plot_ntuple(TString infile="isingsteps.root",
	    TString savedir="NOSAVE")
{
  TFile f(infile.Data());

  TNtuple *nt = (TNtuple*) f.Get("nt");

  TCanvas *cnv1 = new TCanvas("ising1", "ising1", 
			      1000, 800);

  // determine axis scale
  double min = 0;
  double max = ( nt->GetMaximum("n") 
		 + 0.05*(nt->GetMaximum("n") - nt->GetMinimum("n")) );

  TH2F * frame = new TH2F("frame",";steps;magnetization (%)",
			  100, min, max,   100, -1.1, 1.1);
  
  frame->SetStats(0);
  frame->DrawCopy();
  
  int imax = nt->GetMaximum("r");
  imax = 1;
  for ( int i=0; i<imax; i++) {
    nt->Draw("m:n",
	     Form("n<20000000 && r==%d",i),"l,same");
  }

  if ( savedir != "NOSAVE" ) {
    cnv1->SaveAs(savedir + "/magn_vs_steps.eps");
  }

}
