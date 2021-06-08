
void states(TString savedir="NOSAVE")
{

  if ( gROOT->LoadMacro("randwalk_ising.C+") &&
       gROOT->LoadMacro("../examples/mcint/ising/randwalk_ising.C+") ) {

    cout << "loading randwalk_ising.C failed" << endl;
    return;
  }


  gRandom->SetSeed(43);


  randwalk_ising rwgen;

  const int nstates = 11;
  TCanvas* cnv[nstates];

  for (int s=0; s<nstates; s++) {

    switch (s) {
    case  0:   rwgen.reset(0.0);     break;
    case  1:   rwgen.reset(1.0);     break;
    case  2:   rwgen.reset(0.5);     break;
    case  3:   rwgen.next(100);      break; // 100 sample
    case  4:   rwgen.next(1900);     break;
    case  5:   rwgen.next(28000);    break;
    case  6:   rwgen.next(70000);    break;
    case  6:   rwgen.next(400000);   break; // 500k samples
    case  7:
    case  8:
    case  9:
      rwgen.next(500000);  
      break;      

    case 10:
      gRandom->SetSeed(42);
      rwgen.next(50000000);
    }      
      
    cnv[s] = new TCanvas(Form("ising_state_%d", s),
			 Form("ising_state_%d", s),
			 1000, 1000);
    
    rwgen.draw_state();

    if (savedir != "NOSAVE") {
      cnv[s]->SaveAs(savedir + "/" + cnv[s]->GetTitle() + ".eps");
    }
  }

}
