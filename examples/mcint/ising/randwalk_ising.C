
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TNtuple.h>
#include <TF2.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TRandom.h>
#include <TMath.h>
#include <TList.h>

#include <math.h>
#include <iostream>
#include <iomanip>

using namespace std;

// =========================================================================
// randwalk_ising: 
//
// Random walk for the Ising model using the Metropolis Monte-Carlo
// algorithm.
//
class randwalk_ising
{

public:  

    randwalk_ising(int LL=64, 
                   double JJ=1.0, double mumu=1.0, 
                   double BB=0.0, double TT=0.1);
  
    // set parameters for B-field and temperature
    void setB(double b) { B = b; calc_transprob(); } 
    void setT(double t) { T = t; calc_transprob(); }

    double getB() { return B; }
    double getT() { return T; }

    void set_spin(int i, double x) 
    { s[i] = x; calc_nup_nal(); }

    // initialize state with random spins
    void reset(double prob_up = 0.5);

    // proceed to next steps (update position nsteps times
    void next(int nsteps=1);

    // run simulation and fill histogram or ntuple
    struct stats_t {
        double E_mean;
        double E_variance;
        double m_mean;
        double m_variance;
    };

    stats_t run(int runid, int nsteps, 
                THnSparse* hist=NULL, TNtuple* nt=NULL,
                int stepsize=1);

  
    // calculate energy of current state
    double energy()
    { return J*double(2*n-2*nal) + mu*B*double(n-2*nup); }
    //{ return J*double(n-2*nup) + mu*B*double(n-2*nal); }


    // calculate magnetisation of current state
    double magnetisation()
    { return double(2*nup-n)/double(n); }

    // change of energy for a single spin-flip
    double delta_energy(bool spin, int nnnup);

    double delta_energy(int pos)
    { return delta_energy(s[pos], nnup(pos)); }

    //double delta_energy_old(int i);

    // draw the current state in 
    void draw_state(TString title="AUTO");


    // leftover from playing with changing more than one spin at a time
    // double delta_energy(int nf, int *idx);
    // void next(double mean);
    // void next(double mean, int nsteps);



private:

    // ---------------------------------------------------------------------
    // (macroscopic) parameters of simulation
    //
    int L;      // lattice size
    int n;      // number of nodes in lattice

    double mu;  // magnetic moment
    double J;   // magnetic coupling
    double B;   // magnetic field
    double T;   // temperature

    // cache of transition probabilities to speed up next
    double transprob_u[5]; // prob ratios for spins in up with 0..4 up neighbors
    double transprob_d[5]; // prob ratios for spins in up with 0..4 up neighbors
    void calc_transprob();

    // ---------------------------------------------------------------------
    // state information

    // vector of bools: true is spin up, false is spin down
    vector < bool > s;

    // calculate the number of nearest neighbors to point i with spin up
    int nnup(int i);

    // keep the numbers of up spins and aligned spin pairs
    int nup;   // number of spins in up state 
    int nal;   // number of aligned nearest neighbors <ij>: s[i]=s[j]

    void calc_nup_nal();  // (re-)calculate numbers 
    void check_nup_nal(); // check: recalculate and compare with old numbers


    
    // ---------------------------------------------------------------------
    // an empty frame
    TH2F* frame;

};


// -------------------------------------------------------------------------
// constructor: set initial conditions and PDF
randwalk_ising::randwalk_ising(int LL, 
			       double JJ, double mumu, 
			       double BB, double TT)
    : L(LL), n(LL*LL), mu(mumu), J(JJ), B(BB), T(TT)
{
    s.resize(n);
    reset();

    calc_transprob();
  
    // an empty histogram to be used as a frame for state plots
    frame = new TH2F("frame", ";x;y", 
                     L, -0.5, double(L)-0.5, 
                     L, -0.5, double(L)-0.5);
    frame->SetStats(0);
}


// -------------------------------------------------------------------------
void randwalk_ising::reset(double prob_up)
{
    for (int i=0; i<n; i++) {
        s[i] = (gRandom->Uniform(0.,1.) <= prob_up) ? true : false;
    }

    calc_nup_nal();
}



// -------------------------------------------------------------------------
void randwalk_ising::calc_transprob()
{
    // keep an array of transition probablities
 
    for (int i=0; i<5; i++) {
        transprob_u[i] = exp ( - delta_energy(true,  i) / T );
        transprob_d[i] = exp ( - delta_energy(false, i) / T );
    }
}

// -------------------------------------------------------------------------
void randwalk_ising::calc_nup_nal()
{
    // reset variables
    nup = 0; 
    nal = 0; 

    for(int i=0; i<n; i++) {
        if ( s[i] ) nup++;
        
        if ( s[i] == s[(i+1) % n] ) nal++; // right neighbor
        if ( s[i] == s[(i+L) % n] ) nal++; // neighbor below

        // pairs with right and top neighbors are calculated 
        // at the step of the respective neighbors
    }
}

// -------------------------------------------------------------------------
void randwalk_ising::check_nup_nal()
{
    int old_nup = nup; // number of up spins
    int old_nal = nal; // number of aligned spin pairs

    calc_nup_nal();

    if (old_nup != nup) {
        cerr << "Mismatch: number of up spins: " 
             << old_nup  << " != " << nup << endl;
    }
    
    if (old_nal != nal) {
        cerr << "Mismatch: number of aligned spin pairs: " 
             << old_nal  << " != " << nal << endl;
    }
}


// -------------------------------------------------------------------------
void randwalk_ising::draw_state(TString title)
{
    TMarker mrkup;
    mrkup.SetMarkerStyle(20);
    mrkup.SetMarkerColor(kRed);
    mrkup.SetMarkerSize(1.2);

    TMarker mrkdown;
    mrkdown.SetMarkerStyle(24);
    mrkdown.SetMarkerColor(kBlue);
    mrkdown.SetMarkerSize(1.0);

    if (title == "AUTO") {
        title = Form("N_{up} = %d    N_{down} = %d     E = %.0f     m = %.0f%%",
                     nup, n-nup, 
                     energy(), magnetisation()/double(n)*100.);
    }

    frame->SetTitle(title);
    frame->DrawCopy();

    for (int i=0; i<n; i++) {
        double x = i % L;
        double y = i / L;

        if (s[i]) {
            mrkup.DrawMarker(x,y);
        } else {
            mrkdown.DrawMarker(x,y);
        }
    }
      
}
  

// -------------------------------------------------------------------------
// determine the number of next neighbors with spin up
int randwalk_ising::nnup(int i)
{
    int nnu = 0; // number of neighbors in up

    if ( s[(i - L + n) % n] ) nnu++ ; // above
    if ( s[(i + L + n) % n] ) nnu++ ; // below
    if ( s[(i - 1 + n) % n] ) nnu++ ; // left 
    if ( s[(i + 1 + n) % n] ) nnu++ ; // right

    return nnu;
}


double randwalk_ising::delta_energy(bool spin, int nnnup)
{
    // spin flips up -> down
    // (all signs inverted for spin down -> up)

    // spin: original spin at position
    // nnnup: number of next neighbor spins in up

    // energy from -J*s[i]*s[j]
    // old state:   - J * nnup + J * nndown
    // new state:   + J * nnup - J * nndown
    // difference:  2J * nnup - 2J * nndown = 2J(2nnup-4) = 4J(nnup-2)
    
    // energy from -mu*B*s[i]
    // old state:  -mu*B
    // new state:  +mu*B
    // difference: +2*mu*B

    if ( spin ) {
        return + (J*double(4*(nnnup - 2)) + 2.*mu*B);
    } else {
        return - (J*double(4*(nnnup - 2)) + 2.*mu*B);
    } 
}


// -------------------------------------------------------------------------
// proceed to the next step - one flip at a time
void randwalk_ising::next(int nsteps)
{

    for (int i=0; i<nsteps; i++) {

        int f = gRandom->Integer(n);
      
        // calculate probability ratio between new and old position
        //double w = exp( - delta_energy(f) / T );


        // a lookup table for the transition probabilites can speed
        // things up
        double w;
        if (s[f] ) {
            w = transprob_u[nnup(f)];
        } else {
            w = transprob_d[nnup(f)];
        }
        
        // cout << "trial: flip s[" << f << "]: "
        //      << (s[f]?"UP":"DOWN") << " -> " << (!s[f]?"UP":"DOWN")
        //      << "   nnup=" << nnup(f)
        //      << "   dE=" << delta_energy(f)
        //      << "  w=" << setprecision(10) << w << endl;

        if (w < 1) {
            if ( gRandom->Uniform(0.,1.) > w ) {
                // w is not larger than one, and we did not accept
                // with probability w -> REJECT

                //cout << " -> REJECT" << endl;
                continue;
            }
        }

        // not rejected -> ACCEPT new state and update

        //cout << " -> ACCEPT" << endl;
        s[f] = !s[f];

        if ( s[f] ) {
            // flipped: down -> up
            nup++; // one more up now

            // aligned neighbors in old state:  nndown = 4-nnup
            // aligned neighbors in new state:  nnup
            // difference:                      2nnup - 4
            nal += 2*nnup(f) - 4;

        } else {
            // flipped: up -> down
            nup--;
            nal -= 2*nnup(f) - 4;
        }


    } // for i
}


// -------------------------------------------------------------------------
// run simulation and fill histogram or ntuple
randwalk_ising::stats_t randwalk_ising::run(int runid, int nsteps, 
                                            THnSparse* hist, TNtuple* nt,
                                            int stepsize)
{

    stats_t stats;
    

    double buf[4];
    buf[0] = B;
    buf[1] = T;

    double sum_e  = 0.0;
    double sum_e2 = 0.0;
    double sum_m  = 0.0;
    double sum_m2 = 0.0;
    
    for (int i=0; i<nsteps; i++) {

        next(stepsize-1);

        sum_e  += energy();
        sum_e2 += energy()*energy();
        sum_m  += magnetisation();
        sum_m2 += magnetisation()*magnetisation();

        if (hist) {
            buf[2] = energy();
            buf[3] = magnetisation();

            hist->Fill(buf);
        }

        if (nt) {
            nt->Fill(runid, i, B, T, energy(), magnetisation());
        }
    }

    stats.E_mean     = sum_e  / double(nsteps);
    stats.E_variance = sum_e2 / double(nsteps) - stats.E_mean*stats.E_mean;

    stats.m_mean     = sum_m  / double(nsteps);
    stats.m_variance = sum_m2 / double(nsteps) - stats.m_mean*stats.m_mean;

    return stats;
}



// Local Variables:
//   mode: c++
//     c-basic-offset: 4
//     indent-tabs-mode: nil
// End:

