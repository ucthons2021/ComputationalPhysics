
#include <TH1.h>
#include <TH2.h>
#include <TF2.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TRandom.h>
#include <TMath.h>


#include <math.h>
#include <iostream>
#include <unistd.h>

using namespace std;



TF2* create_probfct(TString probfct="ex1");
void draw_probfct(TF2* prob, TString savedir="NOSAVE");

void run(TF2* prob, 
	 int nsteps=100000, 
	 int nskip=1000, 
	 TString savedir="NOSAVE");


void trace(TF2* prob, 
	   int nsteps = 1000, 
	   bool singlestep = false, 
	   TString savedir="NOSAVE");




void metropolis(TString probfct, TString mode, TString savedir="NOSAVE")
{
  TF2* prob = create_probfct(probfct);

  if ( mode == "trace" ) {

    trace(prob, 50, true);

  } else if ( mode == "movie" ) {

    trace(prob, 500, false);

  } else if ( mode == "probfct" ) {

    draw_probfct(prob, savedir);

  } else if ( mode == "pos" ) {

    run(prob, 10000000, 100, savedir);

  }

}



// =========================================================================
// randwalk_metropolis: 
// class that implements a random walk using the Metropolis algorithm
class randwalk_metropolis
{

public:  

  // constructor: set initial conditions and PDF
  randwalk_metropolis(TF2* p, double x0=0.0, double y0=0.0, double s=0.5)
    : prob(p), stepsize(s)
  {
    pos[0] = x0;
    pos[1] = y0;

    trial[0] = 0.;
    trial[1] = 0.;
  }
  
  // proceed to the next step
  void next()
  {

    // generate trial position
    trial[0] = gRandom->Gaus(pos[0],stepsize);
    trial[1] = gRandom->Gaus(pos[1],stepsize);

    // calculate probability ratio between new and old position
    double w = prob->Eval(trial[0],trial[1]) / prob->Eval(pos[0],pos[1]);

    if (w >= 1) {
      // trial position is more likely -> ACCEPT
      pos[0] = trial[0];
      pos[1] = trial[1];
    } else {
      // trial position is less likely ->
      //   ACCEPT with probability w, else reuse current position
      double r = gRandom->Uniform(0.,1.);
      if ( r <= w ) {
	// accept
	pos[0] = trial[0];
	pos[1] = trial[1];
      } else {
	// reject - reuse current position
      }

    }

  }


  // get current x/y position
  double x() { return pos[0]; }
  double y() { return pos[1]; }

  // get last trial position
  double xt() { return trial[0]; }
  double yt() { return trial[1]; }

private:
  double pos[2];
  double trial[2];

  double stepsize;

  TF2*  prob;

};

TF2* create_probfct(TString probfct)
{

  // ------------------------------------------------------------------------
  // build the weight function - unnormalized PDF
  //
  TF2* prob = NULL;
  
  if ( probfct == "ex1" ) {
    // simple example
    prob = new TF2("ex1", "exp(-(x-1.0)*(x-1.0) - (y+2.0)*(y+2.0))", 
		   -5., 5., -5., 5.);
    prob->SetTitle("exp(-(x-1)^{2}-(y+2)^{2});x;y");
  } else if ( probfct == "ex2" ) {
    // not so simple example
    prob = new TF2("ex2", 
		   "exp(-(x-1.0)*(x-1.0)) * exp(-5.3*(y-0.2*x*x +2.3)*(y-0.2*x*x +2.3))", 
		   -5., 5., -5., 5.);
  } else {
    // user defined probability
    prob = new TF2("custom", probfct, -5., 5., -5., 5.);
  }

  return prob;
}




void draw_probfct(TF2* prob, TString savedir)
{

  TLatex txt;
  TMarker mrk;
  
  // set plot granularity
  prob->SetNpx(100);
  prob->SetNpy(100);

  // draw the function
  TCanvas *cnv_prob = new TCanvas("metropolis_prob", 
				  "randwalk prob distribution",
				  800,800);
  prob->Draw("colz");



  if ( savedir != "NOSAVE" ) {
    cnv_prob->SaveAs(savedir + "/metropolis_" + prob->GetName() + "_probfct.png");
  }
  
  if ( TString(prob->GetName()) == "ex1" ) {
    txt.SetTextAlign(11);
    mrk.SetMarkerStyle(2);
    mrk.SetMarkerSize(2.0);

    float x=1.0, y=-2.0;
    mrk.DrawMarker(x,y);
    txt.DrawText(x+0.1,y+0.1, Form("%4.2f",prob->Eval(x,y)));

    x = 0.7; y = -1.1; 
    mrk.DrawMarker(x,y);
    txt.DrawText(x+0.1,y+0.1, Form("%4.2f",prob->Eval(x,y)));

    x = 2.3; y = -1.2; 
    mrk.DrawMarker(x,y);
    txt.DrawText(x+0.1,y+0.1, Form("%4.2f",prob->Eval(x,y)));

    x = 1.0; y = 3.0; 
    mrk.DrawMarker(x,y);
    txt.DrawText(x+0.1,y+0.1, Form("%.2e",prob->Eval(x,y)));

    x = 2.0; y = 0.7; 
    mrk.DrawMarker(x,y);
    txt.DrawText(x+0.1,y+0.1, Form("%.2e",prob->Eval(x,y)));

    
    if ( savedir != "NOSAVE" ) {
      cnv_prob->SaveAs(savedir + "/metropolis_" + prob->GetName() + "_probfct_values.png");
    }
  }    

}


void run(TF2* prob, int nsteps, int nskip, TString savedir)
{
  // ------------------------------------------------------------------------
  // prepare canvas and histogram to be used for positions of the random walk
  TCanvas *cnv_pos = new TCanvas("metropolis_pos", "positions in random walk", 
				 800, 800);

  TH2F *hpos = new TH2F("frame", "positions in random walk",
			 200, -5., 5., 200, -5., 5.);
  hpos->SetStats(0);

  // ------------------------------------------------------------------------
  // run the metropolis algorithm
  //
  
  // create Metropolis random walk object
  randwalk_metropolis rwgen(prob, -4., 4.);

  // skip the first nskip steps
  for (int i=0; i<nskip; i++) {
    rwgen.next();
  }

  // now run for nsteps steps
  for (int i=0; i<nsteps; i++) {
    rwgen.next();
    hpos->Fill(rwgen.x(), rwgen.y());
  }

  // ------------------------------------------------------------------------

  hpos->Draw("colz");

  if ( savedir != "NOSAVE" ) {
    cnv_pos->SaveAs(Form("%s/metropolis_%s_pos.eps", 
			   savedir.Data(), prob->GetName()));
  }


  if ( TString(prob->GetName()) == "ex1" ) {
    TLatex txt;
    txt.SetTextAlign(11);

    TMarker mrk;
    mrk.SetMarkerStyle(2);
    mrk.SetMarkerSize(2.0);

    float x=1.0, y=-2.0;
    mrk.DrawMarker(x,y);
    txt.DrawText(x+0.1,y+0.1, Form("%4.2f",prob->Eval(x,y)));

    x = 0.7; y = -1.1; 
    mrk.DrawMarker(x,y);
    txt.DrawText(x+0.1,y+0.1, Form("%4.2f",prob->Eval(x,y)));

    x = 2.3; y = -1.2; 
    mrk.DrawMarker(x,y);
    txt.DrawText(x+0.1,y+0.1, Form("%4.2f",prob->Eval(x,y)));

    x = 1.0; y = 3.0; 
    mrk.DrawMarker(x,y);
    txt.DrawText(x+0.1,y+0.1, Form("%.2e",prob->Eval(x,y)));

    x = 2.0; y = 0.7; 
    mrk.DrawMarker(x,y);
    txt.DrawText(x+0.1,y+0.1, Form("%.2e",prob->Eval(x,y)));

    
    if ( savedir != "NOSAVE" ) {
      cnv_pos->SaveAs(Form("%s/metropolis_%s_pos_values.eps", 
			   savedir.Data(), prob->GetName()));
    }
  }    

}




void trace(TF2* prob, int nsteps, bool singlestep, TString savedir)
{
  //int argc = 0;
  //TApplication theApp("App",&argc, NULL);

  // ------------------------------------------------------------------------
  // prepare graphics output of Metropolis random walk
  //

  // create a canvas to be used for positions of the random walk
  TCanvas *cnv_trace = new TCanvas("metropolis_trace", "positions in random walk", 
				 800, 800);

  //cnv_trace->Draw();

  //cnv_trace->RaiseWindow();
  //gVirtualX->MapRaised(cnv_trace->GetId());

  // empty frame for random walk positions
  TH2F *frame = new TH2F("frame", "positions in random walk",
			 200, -5., 5., 200, -5., 5.);
  frame->SetStats(0);
  
  // TGraph with the first few positions
  TGraph* gr = new TGraph(1);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(kRed);
  
  TMarker mrk;
  mrk.SetMarkerStyle(21);
  mrk.SetMarkerColor(kBlue);

  // ------------------------------------------------------------------------
  // run the metropolis algorithm
  //

  // create Metropolis random walk object
  cout << "Generating random walk..." << flush;
  randwalk_metropolis rwgen(prob, -4., 4., 0.5);
  cout << " done." << endl;
  
  // save starting point in graph
  gr->SetPoint(0,rwgen.x(), rwgen.y());

  // do the random walk and trace the steps
  for (int i=0; i<nsteps; i++) {

    rwgen.next();

    gr->Set(i+1);
    gr->SetPoint(i,rwgen.x(), rwgen.y());

    frame->Draw();
    gr->Draw("pl");

    if ( rwgen.x() != rwgen.xt() && rwgen.y() != rwgen.yt() ) {
      mrk.DrawMarker(rwgen.xt(), rwgen.yt());
    }
    
    cnv_trace->Modified();
    cnv_trace->Update();
    //cnv_trace->Draw();
    //std::this_thread::sleep_for( std::chrono::seconds(3) );    

    usleep(100);
    
    if ( savedir != "NOSAVE" ) {
      //if ((i<20) || (i<100 && i%5==0) || (i%100 == 0) ) {
      cnv_trace->SaveAs(Form("%s/metropolis_trace_%s_%04d.eps", 
			     savedir.Data(), prob->GetName(), i));
      //}
    }

    if (singlestep) {
      cout << "wait for input" << endl;
      char tmp[100];
      cin.getline(tmp,100);
    }
  }

}






// void randwalk_r2avg()
// {

//   int ntrials = 100;
//   int nsteps = 10000;

//   double sum_mean_rsq = 0.;
  
//   TH1F* hrsq = new TH1F("hrsq", "< R^{2} >", 500, 0., 5000.);


//   for (int t=0; t<ntrials; t++) {

//     randwalk_generator rwgen;
//     double sumr2 = 0.;

//     for (int i=0; i<nsteps; i++) {
//       rwgen.next();
//       sumr2 += rwgen.x()*rwgen.x() + rwgen.y()*rwgen.y();
//     }
    
//     sum_mean_rsq += sumr2 / double(nsteps);
//     hrsq->Fill(sumr2 / double(nsteps));
//   }

//   cout << "R_avg^2 = " << ( sum_mean_rsq / double(ntrials) ) << endl;

//   TCanvas *cnv_r2avg = new TCanvas("randwalk_r2avg", "average R^2", 
// 				   1000, 800);
//   hrsq->Draw();
//   cnv_r2avg->SaveAs("randwalk_r2avg.pdf");
  
// }




// void randwalk_rrms()
// {


//   int ntrials = 100;
//   int nsteps = 10000;


//   TH1F* hrrms = new TH1F("hrrms", 
// 			 "R_{rms} vs. #sqrt{N};#sqrt{N};R_{rms}", 
// 			 int(sqrt(nsteps)), 0.5, sqrt(nsteps)+0.5);
  

//   for (int t=0; t<ntrials; t++) {

//     randwalk_generator rwgen;
//     double sumr2 = 0.;

//     for (int i=0; i<nsteps; i++) {
//       rwgen.next();
//       sumr2 += rwgen.x()*rwgen.x() + rwgen.y()*rwgen.y();
    
//       double s = floor(sqrt(i+1));
      
//       if ( fabs( s*s - double(i+1) ) < 0.0001 ) {
// 	//hrrms->Fill(s, sumr2 / double(i+1));
// 	hrrms->Fill(s, (rwgen.x()*rwgen.x() + rwgen.y()*rwgen.y()) );
//       }
//     }

//   }
  
//   hrrms->SetStats(0);
//   for (int b=1; b<=hrrms->GetXaxis()->GetNbins(); b++) {
//     hrrms->SetBinContent(b, sqrt(hrrms->GetBinContent(b) / ntrials));
//   }
  
//   TCanvas *cnv_rrms = new TCanvas("randwalk_rrms", "R_{rms}", 
//   				   1000, 800);
//   hrrms->Draw("l");
//   cnv_rrms->SaveAs("randwalk_rrms.pdf");
  
// }



// void randwalk()
// {
//   randwalk_pos();
//   randwalk_r2avg();
//   randwalk_rrms();
// }

