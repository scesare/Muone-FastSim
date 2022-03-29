#include <iostream>
#include <sstream>
#include "MuEtree.h"
#include "Utils.h"

using namespace std;
using namespace MuE;

MCpara::MCpara(ifstream & input_file): MCpara() {
  // Read the header section
  string line, key, dump, str1, str2;
  istringstream stream;

  cout << "Start reading header section" <<endl;
  stream = input_line(input_file);
  stream >> key;
  if (key != "<header>") {
    cout << "*** ERROR: unexpected format for header section." << endl;
    exit(100);
  }
  
  stream = input_line(input_file);
  stream >> dump >> dump >> str1 >> str2;
  program_version = str1+" "+str2;
  
  stream = input_line(input_file);
  stream >> dump >> dump >> process_ID >> dump >> dump >> running_on;

  stream = input_line(input_file);
  start_time = stream.str();

  // Number of requested events
  stream = input_line(input_file);
  stream >> dump >> Nevreq;

  // Unweighted (wgt=1) events ?  
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> UNWGT;

  // generator Mode
  stream = input_line(input_file);
  stream >> dump >> dump >> Mode;

  // Initial seed for Random numbers 
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> rnd_ext >> rnd_int;

  stream = input_line(input_file);
  
  // Nominal Muon Beam energy
  stream = input_line(input_file);
  stream = input_line(input_file);  
  stream >> Ebeam; 

  // Beam energy Gaussian spread 
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> EbeamRMS;

  // muon charge
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> charge_mu;

  // muon mass
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> mass_mu;

  // electron mass
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> mass_e;

  // (fine structure constant)^-1  
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> invalfa0;  

  // photon cutoff (technical parameter)
  stream = input_line(input_file);
  stream >> dump >> dump >> dump >> dump >> k0cut;

  // Minimum Energy of outgoing electron
  stream = input_line(input_file);
  stream >> dump >> dump >> dump >> dump >> dump >> dump >> Emin_e;

  stream = input_line(input_file);
  
  // Cross section normalization factor 
  stream = input_line(input_file);
  stream = input_line(input_file);  
  stream >> Wnorm; 

  // Assumed maximum weight (for unweighted generation)
  stream = input_line(input_file);
  stream = input_line(input_file);  
  stream >> Wmax; 
  
  stream = input_line(input_file);
  stream >> key;
  if (key == "</header>") {
      cout << "End reading header section." <<endl;

      cout<<"\n"<<"========================================================================"<< endl;
      cout<<"MuE generator: "<< program_version << endl;
      string strwgt = UNWGT ? "Unweighted" : "Weighted";
      cout<< strwgt << " events generation" << endl;
      cout<<"Mode : "<< Mode << endl;
      cout<<"muon beam energy        = "<< Ebeam << " GeV" << endl;
      cout<<"RMS beam energy spread  = "<< EbeamRMS << " GeV" << endl;
      cout<<"muon charge             = "<< charge_mu << endl;
      cout<<"electron mass           = "<< mass_e << " GeV" << endl;
      cout<<"muon mass               = "<< mass_mu << " GeV" << endl;
      cout<<"1/alpha                 = "<< invalfa0 << endl;
      cout<<"Minimum electron energy = "<< Emin_e << " GeV" << endl;
      cout<<"initial random seeds    = "<< rnd_ext << " " << rnd_int << endl;
      cout<<"photon cutoff k_0       = "<< k0cut << " GeV" << endl;
      cout<<"initial Wmax            = "<< Wmax <<endl;
      cout<<"Initial normalization Wnorm = "<< Wnorm << " ub" << endl;

      cout<<"========================================================================"<< endl;

  } else {
    cout << "*** ERROR: unexpected format for the header section." << endl;
    exit(200);
  }
}

MCpara::MCpara(const Setup & s): 
  program_version(s.program_version),
  process_ID(s.process_ID),
  running_on(s.running_on),
  start_time(s.start_time),
  Nevreq(s.Nevreq),
  UNWGT(s.UNWGT),
  Mode(s.Mode),
  rnd_ext(s.rnd_ext), 
  rnd_int(s.rnd_int),
  Ebeam(s.Ebeam), 
  EbeamRMS(s.EbeamRMS),
  charge_mu(s.charge_mu),
  mass_mu(s.mass_mu), 
  mass_e(s.mass_e),
  invalfa0(s.invalfa0),
  k0cut(s.k0cut), 
  Emin_e(s.Emin_e), 
  Wnorm(s.Wnorm), 
  Wmax(s.Wmax)
{};


MCstat::MCstat(ifstream & input_file): MCstat() {
  
  // Read the footer section
  string line, key, dump, str1, str2;
  istringstream stream;
  cout <<endl<< "Start reading the footer section" <<endl;

  // estimated cross section for the generated process
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> Xsec >> dump >> XsecErr;

  // total number of weights
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> Nwgt;

  // true maximum weight at the end 
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> WmaxTrue;

  // number of weights greater than the assumed maximum (Wmax)
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> Nwgt_OverMax;

  // number of negative weights
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> Nwgt_Negative;

  // estimated bias on cross section due to weights greater than Wmax (in unweighted generation)
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> Xsec_OverMax >> dump >> Xsec_OverMax_Err;

  // estimated cross section contribution from negative weights
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> Xsec_Negative >> dump >> Xsec_Negative_Err;

  // final sum of weights and squared weights
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> Swgt >> SQwgt;

  // final sum of negative weights and squared weights
  stream = input_line(input_file);
  stream = input_line(input_file);
  stream >> Swgt_Negative >> SQwgt_Negative;

  stream = input_line(input_file);
  stream >> key;
  if (key == "</footer>") {
      cout << "End reading the footer section." <<endl;

      cout<<"\n"<<"========================================================================"<< endl;
      //      cout<<"N generated events        = "<< Nevgen << endl;
      cout<<"N weights                 = "<< Nwgt << endl;
      cout<<"N negative weights        = "<< Nwgt_Negative << endl;
      cout<<"N weights above Wmax      = "<< Nwgt_OverMax << endl;
      cout<<"True Max weight           = "<< WmaxTrue <<endl;
      cout<<"Cross section             = "<< Xsec << " +/- " << XsecErr << endl;
      cout<<"Cross section (negative)  = "<< Xsec_Negative << " +/- " << Xsec_Negative_Err << endl;
      cout<<"Cross section (above max) = "<< Xsec_OverMax << " +/- " << Xsec_OverMax_Err << endl;
      cout<<"========================================================================"<< endl;

  } else {
    cout << "*** ERROR: unexpected format for the footer section." << endl;
    exit(600);
  }

}


Setup::Setup(const MCpara & m, const MCstat & s):
  program_version(m.program_version),
  process_ID(m.process_ID),
  running_on(m.running_on),
  start_time(m.start_time),
  Nevreq(m.Nevreq),
  UNWGT(m.UNWGT),
  Mode(m.Mode),
  rnd_ext(m.rnd_ext), 
  rnd_int(m.rnd_int),
  Ebeam(m.Ebeam), 
  EbeamRMS(m.EbeamRMS),
  charge_mu(m.charge_mu),
  mass_mu(m.mass_mu), 
  mass_e(m.mass_e),
  invalfa0(m.invalfa0),
  k0cut(m.k0cut), 
  Emin_e(m.Emin_e), 
  Wnorm(m.Wnorm), 
  Wmax(m.Wmax),
  MCsums(s)
{}


bool Event::Read(ifstream & input_file) {

// Read a NLO MC event

  UInt_t nrows, nphot;

  string line, key, dump, str1, str2;
  istringstream stream;
  
  getline(input_file, line);
  stream = istringstream(line);
  stream >> key;
  if (key == "<footer>") return false;
  else if (key != "<event>") {
    cout << "*** ERROR: unexpected format for input event (begin)" << endl;
    exit(300);
  }
  
  // read run number
  stream = input_line(input_file);
  stream >> RunNr;
  
  // read event number
  stream = input_line(input_file);
  stream >> EventNr;
  
  // read number of final state particles (-> nr of photons)
  stream = input_line(input_file);
  stream >> nrows;
  nphot = nrows!=0 ? nrows-2 : 0;
  
  // read MC weights (with full/no/leptonic running)
  stream = input_line(input_file);
  stream  >> wgt_full >> wgt_norun >> wgt_lep;
  
  // read LO weight (with full running)
  stream = input_line(input_file);
  stream >> wgt_LO;
  
  // read incoming muon energy
  stream = input_line(input_file);
  stream >> E_mu_in;
  
  // read scattered muon
  stream = input_line(input_file);
  stream >> P_mu_out.E >> P_mu_out.px >> P_mu_out.py >> P_mu_out.pz ;
  
  // read scattered electron
  stream = input_line(input_file);
  stream >> P_e_out.E  >> P_e_out.px  >> P_e_out.py  >> P_e_out.pz ;
  
  // read photons (when nphot>0)
  photons.clear();
  for (UInt_t i=0; i<nphot; ++i) {
    MuE::P4 P_photon;
    stream = input_line(input_file);
    stream >> P_photon.E >> P_photon.px >> P_photon.py >> P_photon.pz;
    photons.push_back(P_photon);
  }
  
  stream = input_line(input_file);
  stream >> key;
  if (key != "</event>") {
    cout << "*** ERROR: unexpected format for input event (end)" << endl;
    exit(400);
  }

  return true;
}
