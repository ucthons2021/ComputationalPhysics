
void mixedstates(TString savedir="NOSAVE")
{

  gRandom->SetSeed(42);

  if ( gROOT->LoadMacro("randwalk_ising.C+") &&
       gROOT->LoadMacro("../examples/mcint/ising/randwalk_ising.C+") ) {

    cout << "loading randwalk_ising.C failed" << endl;
    return;
  }

  TCanvas* cnv = new TCanvas("ising_mixedstate", "ising_mixedstate",
			     1000, 1000);
  
  randwalk_ising rwgen;

  int nstates = 3;

  // this is in principle useless, but it takes me to a cool 
  // first picture
  rwgen.reset(0.5);
  rwgen.next(10000000);

  for (int i=0; i<nstates; i++) {

    do {
      rwgen.reset(0.5);
      rwgen.next(10000000);
    } while ( fabs(rwgen.magnetisation()) > 0.95);
    
    
    for (int j=0; j<5; j++) {

      rwgen.next(10000);
      rwgen.draw_state();

      cnv->Update();

      //char tmp[100];
      //cin.getline(tmp,100);

      if (savedir != "NOSAVE") {
	cnv->SaveAs(Form("%s/ising_mixedstate_%d_%d.eps",
			 savedir.Data(), i, j));
      }
    }

  }

}
