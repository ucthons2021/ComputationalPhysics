
// ===========================================================================
// coolheat - Mini-Howto
//
// Starting:
//  - start root
//  - compile and load randwalk_ising.C:   .L randwalk_ising.C+
//  - load coolheat.C:                     .L coolheat.C
//  - coolheat instance with 5 systems:    coolheat c(5)
//
// Initialize:
//  - set B-field and temperature:         c.setB(0.1); c.setT(0.0);
//  - reset state:
//              all spins up:              c.reset(1.0) 
//              all spins down:            c.reset(0.0) 
//              random spins:              c.reset(0.5) 
// 
// Play:
//  - ramp temperature in from current to 5.0, in 500 steps, with 100000
//    iterations of the Metropolis MC for each temperature step:
//      c.rampT(5.0, 500, 100000)
//
//  - same for B-field:
//      c.rampB(5.0, 500, 100000)
//
// ===========================================================================

// void coolheat()
// {
//   if ( gROOT->LoadMacro("randwalk_ising.C+") &&
//        gROOT->LoadMacro("../examples/mcint/ising/randwalk_ising.C+") ) {

//     cout << "loading randwalk_ising.C failed" << endl;
//     return;
//   }
// }

class coolheat
{

public:
  coolheat(int nsystems)
    : nsys(nsystems), cnvmag(new TCanvas), cnvstate(NULL), frame(NULL)
  {

    sys = new randwalk_ising[nsys];

    for (int i=0; i<nsys; i++) {
      sys[i].setB(0.0);
      sys[i].setT(0.3);

      if (i%2) {
	sys[i].reset(1.0);
      } else {
	sys[i].reset(0.0);
      }

      sys[i].next(1000000);
      sys[i].setT(0.0);
    }

    

  }

  double reset(double prob_up = 0.5)
  {
    for (int i=0; i<nsys; i++) {
      sys[i].reset(prob_up);
    }

  }
  void setB(double B)
  { for (int i=0; i<nsys; i++) sys[i].setB(B); }

  void setT(double T)
  { for (int i=0; i<nsys; i++) sys[i].setT(T); }

  rampT(double Tend, TString opt, int Tsteps=500, int MCsteps=100000);
  rampB(double Bend, TString opt, int Bsteps=500, int MCsteps=100000);

  rampT(double Tend, int Tsteps=500, int MCsteps=100000)
  { rampT(Tend, "", Tsteps, MCsteps); }

  rampB(double Bend, int Bsteps=500, int MCsteps=100000)
  { rampB(Bend, "", Bsteps, MCsteps); }
  

  void setStateCanvas(TCanvas* cnv) { cnvstate = cnv; }

protected:
  int nsys;
  randwalk_ising *sys;

  TCanvas* cnvmag;
  TH2F*    framemag;

  TCanvas* cnvstate;
};



void coolheat::rampT(double Tend, TString opt, int Tsteps, int MCsteps)
{



  double Tstart = sys[0].getT();

  if (Tstart == Tend) {
    cout << "error: Tstart = Tend" << endl;
    return;
  }

  TGraph *gr = new TGraph[nsys];
  for (int s=0; s<nsys; s++) {
    gr[s].SetMarkerStyle(7);
  }
  
  double Tmin, Tmax;
  if (Tstart < Tend) {
    Tmin = Tstart; Tmax = Tend;
  } else {
    Tmin = Tend; Tmax = Tstart;
  }


  if (framemag) {
    delete framemag;
  }
   
  cnvmag->cd();

  framemag = new TH2F("framemag", ";T;magnetisation", 
		      100, Tmin, Tmax, 110, -1.1, 1.1);

  framemag->SetStats(0);
  framemag->Draw();

  //if ( ! opt.Contains("same") ) {
  //}

  for (float i=0; i<Tsteps; i++) {
    
    double deltaT = (Tend-Tstart)/double(Tsteps-1);
    double T = Tstart + deltaT*double(i);
    
    for (int s=0; s<nsys; s++) {
      sys[s].setT(T);
      sys[s].next(MCsteps);
      //sys[i].run(0,10000,hist);
    }

    cnvmag->cd();
    framemag->Draw();
    for (int s=0; s<nsys; s++) {
      gr[s].Set(i+1);
      gr[s].SetPoint(i, sys[s].getT(), sys[s].magnetisation());
      gr[s].Draw("pl");
    }
    cnvmag->Update();

    if (cnvstate) {
      cnvstate->cd();
      sys[0].draw_state();
      cnvstate->Modified();
      cnvstate->Update();
    }
  }

  // // -------------------------------------------------------------------
  // cnv_magn_temp = new TCanvas("cnv_magn_temp","cnv_magn_temp",1000,800);

  // h_magn_temp = hist->Projection(3,1);
  // h_magn_temp->SetTitle(Form("%d iterations per temperature step;k_{B}T / J;magnetisation",
  // 			     MCsteps));
  // h_magn_temp->SetStats(0);
  // h_magn_temp->Draw("colz");

  // if (magn_plotname != "") {
  //   cnv_magn_temp->SaveAs(
  // }


  // -------------------------------------------------------------------
  // cnv_enrgy_temp = new TCanvas("cnv_enrgy_temp","cnv_enrgy_temp",1000,800);

  // h_enrgy_temp = hist->Projection(2,1);
  // h_enrgy_temp->SetTitle(";k_{B}T / J;E");
  // h_enrgy_temp->SetStats(0);
  // h_enrgy_temp->Draw("colz");
  

}


void coolheat::rampB(double Bend, TString opt, int Bsteps, int MCsteps)
{

  double Bstart = sys[0].getB();

  TGraph *gr = new TGraph[nsys];
  for (int s=0; s<nsys; s++) {
    gr[s].SetMarkerStyle(7);
  }
  
  double Bmin, Bmax;
  if (Bstart < Bend) {
    Bmin = Bstart; Bmax = Bend;
  } else {
    Bmin = Bend; Bmax = Bstart;
  }

  if (framemag) {
    delete framemag;
  }
   
  cnvmag->cd();

  framemag = new TH2F("framemag", ";B;magnetisation", 
		      100, Bmin, Bmax, 110, -1.1, 1.1);

  framemag->SetStats(0);
  framemag->Draw();

  // if ( ! opt.Contains("same") ) {
  //   TH2F* frame = new TH2F("frame", ";B;magnetisation", 100, Bmin, Bmax, 110, -1.1, 1.1);
  //   frame->SetStats(0);
  //   frame->Draw();
  // }

  for (float i=0; i<Bsteps; i++) {
    
    double deltaB = (Bend-Bstart)/double(Bsteps-1);
    double B = Bstart + deltaB*double(i);
    
    for (int s=0; s<nsys; s++) {
      sys[s].setB(B);
      sys[s].next(MCsteps);
      //sys[i].run(0,10000,hist);
    }

    cnvmag->cd();
    framemag->Draw();
    for (int s=0; s<nsys; s++) {
      gr[s].Set(i+1);
      gr[s].SetPoint(i, sys[s].getB(), sys[s].magnetisation());
      gr[s].Draw("pl");
    }
    gPad->Update();

    if (cnvstate) {
      cnvstate->cd();
      sys[0].draw_state();
      cnvstate->Modified();
      cnvstate->Update();
    }

  }

}
