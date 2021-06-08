
#include <TH1.h>
#include <TH2.h>
#include <TNtuple.h>
#include <TF2.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TRandom.h>
#include <TMath.h>
#include <TFile.h>
#include <TROOT.h>

#include <math.h>
#include <iostream>

using namespace std;





void gen_ntuple(TString filename="isingsteps.root", double T = 0.1, double B = 0.0, 
		int nruns = 50, int ntrials=50000000, int nentries=500)
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
  rwgen.setT(T);

  TFile of(filename.Data(), "RECREATE");
  TNtuple *nt = new TNtuple("nt","nt","r:n:b:t:e:m");
  

  int nfillstep = ntrials / nentries;

  for (int r=0; r<nruns; r++) {

    cout << "Run " << r << ":" << flush;

    //rwgen.reset(gRandom->Uniform(0.,1.));
    rwgen.reset();

    rwgen.run(r, nentries, NULL, nt, ntrials/nentries);
    
  //   // do the random walk and fill the histogram
  //   for (int i=0; i<ntrials; i+=nfillstep) {

  //     rwgen.next(nfillstep);
  //     nt->Fill(r,i,rwgen.energy(), rwgen.magnetisation());
  //     //      cout << "." << flush;

  //   }

    cout << "  E=" << rwgen.energy() 
   	 << "  m=" << rwgen.magnetisation() << endl;

  }
  nt->Write();

}

