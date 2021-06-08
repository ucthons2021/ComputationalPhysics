
void cooling(TString path)
{
  int nruns;
  int MCsteps;
  char savedir[1000];

  TString dir  = gSystem->DirName(path);
  TString file = gSystem->BaseName(path);

  sscanf(file.Data(), "cooling_%d_%d.eps", &nruns, &MCsteps);
  cooling(nruns, MCsteps, 5.0, 100, 0.0, dir);
}


void cooling(int runs=5, int MCsteps=10000, 
	     double Tmax=5.0, int Tbins=100, double B=0.0,
	     TString savedir="NOSAVE")
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
  rwgen.setB(B);
  
  // create n-dim histogram for data
  // Dimensions are:
  //   0: B field
  //   1: temperature
  //   2: energy
  //   3: magnetization
  int  nbins[4] = {  1  , Tbins,   1000,   440};
  double min[4] = { -1. ,   0.0, -10000., -1.1};
  double max[4] = {  1. ,  Tmax,  10000.,  1.1};

  THnSparse* hist = new THnSparseI("hist","hist",4,nbins,min,max);


  double deltaT = Tmax/double(Tbins);

  for (int r=0; r<runs; r++) {

    rwgen.reset();

    for (float T=Tmax - 0.5*deltaT; T>0; T-= deltaT) {
      rwgen.setT(T);
      rwgen.next(MCsteps);
      stats = rwgen.run(0,10000,hist);
    }
  }

  // -------------------------------------------------------------------
  cnv_magn_temp = new TCanvas("cnv_magn_temp","cnv_magn_temp",1000,800);

  h_magn_temp = hist->Projection(3,1);
  h_magn_temp->SetTitle(Form("%d iterations per temperature step;k_{B}T / J;magnetisation",
			     MCsteps));
  h_magn_temp->SetStats(0);
  h_magn_temp->Draw("colz");

  if (savedir != "NOSAVE") {
    cnv_magn_temp->SaveAs(Form("%s/cooling_%d_%d.eps",
			       savedir.Data(), runs, MCsteps));
  }


  // -------------------------------------------------------------------
  // cnv_enrgy_temp = new TCanvas("cnv_enrgy_temp","cnv_enrgy_temp",1000,800);

  // h_enrgy_temp = hist->Projection(2,1);
  // h_enrgy_temp->SetTitle(";k_{B}T / J;E");
  // h_enrgy_temp->SetStats(0);
  // h_enrgy_temp->Draw("colz");
  

}
