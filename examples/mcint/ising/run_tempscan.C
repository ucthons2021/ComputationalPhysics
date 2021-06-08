
void run_tempscan(int nruns=5, double Tmax=5.0, int nsteps=500, 
		  TString rootfilename="tempscan.root")
{

  if ( gROOT->LoadMacro("randwalk_ising.C+") &&
       gROOT->LoadMacro("../examples/mcint/ising/randwalk_ising.C+") ) {

    cout << "loading randwalk_ising.C failed" << endl;
    return;
  }


  // ------------------------------------------------------------------------
  // run the metropolis algorithm
  //
  
  // create Ising/Metropolis random walk object
  randwalk_ising rwgen;
  
  // create n-dim histogram for data
  // Dimensions are:
  //   0: B field
  //   1: temperature
  //   2: energy
  //   3: magnetization
  int nbins[4] = {  1  , nsteps,   1000,   440};
  double min[4] = { -1. ,   0.0, -10000., -1.1};
  double max[4] = {  1. ,  Tmax,  10000.,  1.1};

    
  TFile of(rootfilename.Data(), "RECREATE");

  THnSparse* hist = new THnSparseI("hist","hist",4,nbins,min,max);

  TH2F* hmt = new TH2F("hmt",";k_{B}T / J;magnetisation", 
		       nsteps, 0., Tmax, 440, -1.1, 1.1);

  rwgen.next(10000000);

  cout << "initial conditions: m=" << rwgen.magnetisation()
       << "   E=" << rwgen.energy() << endl;

  //for (int r=0; r<10; r++) {

  double deltaT = Tmax/double(nsteps);

  randwalk_ising::stats_t stats;

  time_t start = time(NULL);
  
  for (int r=0; r<nruns; r++) {

    for (float T=0.5*deltaT; T<Tmax; T+= deltaT) {
      rwgen.setT(T);

      if ( T>2.2 && T< 2.5) {
	for (float TT = T; TT<T+deltaT; TT+=0.05*deltaT) {
	  rwgen.setT(T);
	  rwgen.next(200000);
	}

      } else if ( T>2.0 && T< 2.7) {
	rwgen.next(200000);
      } else {
	rwgen.next(10000);
      }

      stats = rwgen.run(0,10000,hist);
      hmt->Fill(T,stats.m_mean);
    }
    
    for (float T=Tmax - 0.5*deltaT; T>0; T-= deltaT) {
      rwgen.setT(T);
      if ( T>2.2 && T< 2.5) {
	for (float TT = T; TT<T+deltaT; TT+=0.05*deltaT) {
	  rwgen.setT(T);
	  rwgen.next(200000);
	}

      } else if ( T>2.0 && T< 2.7) {
	rwgen.next(200000);
      } else {
	rwgen.next(10000);
      }

      stats = rwgen.run(0,10000,hist);
      hmt->Fill(T,stats.m_mean);
    }

    time_t elapsed = time(NULL)-start;
    double eta = double(nruns - r - 1) / double(r+1) * elapsed;
      
    cout << "elapsed: " << elapsed
	 << " sec   ETA: " << eta
	 << " sec" << endl;

  }

  hist->Write();

  // hmt->Draw();

  // // -------------------------------------------------------------------
  // cnv_magn_temp = new TCanvas("cnv_magn_temp","cnv_magn_temp",1000,800);

  // h_magn_temp = hist->Projection(3,1);
  // h_magn_temp->SetTitle(";k_{B}T / J;magnetisation");
  // h_magn_temp->SetStats(0);
  // h_magn_temp->Draw("colz");



  // // -------------------------------------------------------------------
  // cnv_enrgy_temp = new TCanvas("cnv_enrgy_temp","cnv_enrgy_temp",1000,800);

  // h_enrgy_temp = hist->Projection(2,1);
  // h_enrgy_temp->SetTitle(";k_{B}T / J;E");
  // h_enrgy_temp->SetStats(0);
  // h_enrgy_temp->Draw("colz");
  

}
