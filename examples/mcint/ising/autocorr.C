
//#include <vector>

#include <TGraph.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TCanvas.h>


using namespace std;

void autocorr(TString savedir="NOSAVE")
{

  if ( gROOT->LoadMacro("randwalk_ising.C+") &&
       gROOT->LoadMacro("../examples/mcint/ising/randwalk_ising.C+") ) {

    cout << "loading randwalk_ising.C failed" << endl;
    return;
  }


  gRandom->SetSeed(43);


  randwalk_ising rwgen;
  rwgen.reset(0.5);
  rwgen.setB(0.01);
  rwgen.setT(2.3);
  //rwgen.next(1000000);

  int nauto = 100;
  int nruns = 100;
  TH1F* he  = new TH1F("he","autocorrelation", nauto, 0.5, nauto+0.5);
  //TH2F* he2 = new TH2F("he2","E vs n", 10*nauto, -0.5, nauto-0.5, 200, -8200, 8200);

  TNtuple *nt = new TNtuple("nt", "nt", "c:i:n");

  TGraph *ac = new TGraph(nauto);
  TGraph *gr = new TGraph(nruns*nauto);
  gr->SetMarkerStyle(20);

  float* de  = new float[nauto];
  float* snl = new float[nauto];
  //de.resize(nauto);

  for (int i=0; i<nauto; i++) {
    rwgen.next();
    de[i] = rwgen.energy();
  }
  

  for (int i=0; i<nruns*nauto; i++) {

    for (int m=0; m<nauto;m++) {
      snl[m] = 0.0;
    }

    float sn = 0.0;
    for (int m=0; m<nauto;m++) {
      sn  += de[m];

      for (int n=0; n<nauto;n++) {
	snl[n] += de[m%nauto]*de[(m+n)%nauto];
      }
    }
    
    for (int n=0; n<nauto;n++) {
      float c = ( ( snl[n]/nauto - sn*sn/nauto/nauto ) / 
		  (snl[0]/nauto / sn*sn/nauto/nauto) ); 
      nt->Fill(c,i,n);
      he->Fill(n,c);
    }


    
    rwgen.next(1);

    
    //he2->Fill(i,rwgen.energy());
    gr->SetPoint(i,i,rwgen.energy());
    de[i%nauto] = rwgen.energy();




    //cout << i << endl;
  }

  new TCanvas();
  he->Draw();
  

  new TCanvas();
  //he2->Draw("colz");
  gr->Draw("alp");

  if (savedir != "NOSAVE") {
    //cnv[s]->SaveAs(savedir + "/" + cnv[s]->GetTitle() + ".eps");
  }

}
