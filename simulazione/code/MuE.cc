///////////////////////////////////////////////////////////////
//
// read MuE MC events from root file format 
// and make some simple analysis
//
// G.Abbiendi  5/Jul/2018
///////////////////////////////////////////////////////////////

#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <utility>
#include "TROOT.h"
#include "TRint.h"
#include "TChain.h"
#include "TString.h"
#include "TSystem.h"
#include "MuEtree.h" 
#include "FastSim.h"
#include "Analysis.h"
#include "Utils.h"
#include "Inputs.h"
#include "TRandom.h"
#include "ECALProprieties.h"


using namespace std;
using namespace MuE;

int main(int argc, char* argv[]) {

  //  if (argc != 2) {
  //    cerr << "Usage : "<< argv[0] << " PATH_TO_INPUT_CFG_FILE \n";
  if (argc != 1) {
    //    cerr << "Usage : "<< argv[0] << endl;
    cerr << "argc = " << argc << endl;
    cerr << "argv[0] = "<< argv[0] << endl;
    cerr << "argv[1] = "<< argv[1] << endl;
    //    exit(100);
  }

  bool _debug_ = false;
  bool _debug_inp_ = true;

  MuE::One_Input cfg;
  //  cfg.configure(argv[1], _debug_inp_);
  cfg.configure("input.cfi", _debug_inp_);

  MuE::FS_Input fsi;
  fsi.configure(cfg.fastSim_ifname, _debug_inp_);

  MuE::AN_Input an_input;
  an_input.configure(cfg.analysis_ifname, _debug_inp_);

  // MC filenames to be opened are read in from list on an input file
  ifstream input_file_list(cfg.input_dirs_file);
  string root_ifname = "events-NLO.root";
  string line;
  vector<string> ifnames;
  vector<string> idirnames;
  vector<pair<string,Long64_t> > mc_inputs;

  cerr << "looking for the following input root files : "<<endl;
  unsigned ifile = 0;
  while(getline(input_file_list, line)) {
    idirnames.push_back(line);
    string input_rootfile = line+root_ifname;
    cerr << "\t" << "file "<< ifile+1 << ": " << input_rootfile << endl;
    ifnames.push_back(input_rootfile);
    ifile++;
  }
  Long64_t nfiles = ifnames.size();

  TString workdir(cfg.output_dir);
  int iexist = gSystem->mkdir(workdir);
  if (iexist != 0) {
    cerr << "ERROR: directory "<< workdir << " already exists. Please change the output name." << endl;
    exit(200);
  }
  gSystem->cd(workdir);
  cerr << "Writing analysis results to output directory : "<< workdir <<endl;

  // Set it to kTRUE if you do not run interactively
  gROOT->SetBatch(kTRUE); 
  //gROOT->SetBatch(kFALSE); 

  // Initialize Root application
  TRint* app = new TRint("Root Application", &argc, argv);

  MuE::Setup* params = new MuE::Setup(); 
  MuE::Event* event = new MuE::Event();
  
  TChain chain("MuEtree"); // chain event trees
  TChain parchain("MuEsetup"); // chain parameter trees

  cerr << "chaining trees... " << endl;
  cout << "Input root files: " << endl;
  unsigned ifi = 0;
  for (const auto & ifname : ifnames) {
    ++ifi;
    cout<< "File "<< ifi << ". " << ifname << endl;
    parchain.Add(ifname.data());
    chain.Add(ifname.data());
  }
  Long64_t nhdtot = parchain.GetEntries();
  if (nhdtot != nfiles) {
    cerr << "***ERROR: number of headers = "<<nhdtot
	      <<" does not correspond to the number of chained files = "<<nfiles<<endl;
    exit(300);  
  }
  cerr << "Number of input files         = " << nfiles << endl;
  Long64_t nevtot = chain.GetEntries();
  cerr << "Total number of input events  = " << nevtot << endl;

  parchain.SetBranchAddress("MuEparams", &params);
  chain.SetBranchAddress("MuE", &event);

  if (_debug_) cerr << "parameter chain ... " << endl; 
  MuE::MCpara pargen; // parameters as read from the first tree in the chain
  MuE::MCstat sumc;    // sums over the chained trees
  Double_t WmaxFinal(0);

  for (Long64_t ifil=0; ifil < nfiles; ++ifil) {
    Long64_t ientry = parchain.LoadTree(ifil);
    if (ientry < 0) {
      cerr << "***ERROR in chaining MuEsetup, ifil = "<<ifil<<", ientry = "<<ientry<<endl;
      break;
    }
    parchain.GetEntry(ifil);
    if (_debug_) parchain.Show();
    if (ifil==0) {
      cout<<"\n"<<"========================================================================"<< endl;
      cout<<"MuE generator: "<< params->program_version << endl;
      string strwgt = params->UNWGT ? "Unweighted" : "Weighted";
      cout<< strwgt << " events generation" << endl;
      cout<<"Mode : "<< params->Mode << endl;
      cout<<"muon beam energy        = "<< params->Ebeam << " GeV" << endl;
      cout<<"RMS beam energy spread  = "<< params->EbeamRMS << " GeV" << endl;
      cout<<"muon charge             = "<< params->charge_mu << endl;
      cout<<"electron mass           = "<< params->mass_e << " GeV" << endl;
      cout<<"muon mass               = "<< params->mass_mu << " GeV" << endl;
      cout<<"1/alpha                 = "<< params->invalfa0 << endl;
      cout<<"Minimum electron energy = "<< params->Emin_e << " GeV" << endl;
      cout<<"Initial normalization Wnorm = "<< params->Wnorm << " ub" << endl;
      cout<<"========================================================================"<< endl;
      pargen = params->GetMCpara();
    }
    else {
      bool checkOk = CheckParameters(pargen, params->GetMCpara());
      if (!(checkOk)) exit(400);
    }

    mc_inputs.push_back(make_pair(idirnames.at(ifil), params->MCsums.Nwgt));

    // summary printouts for a single file
    cout<<"\n"<<"Input file : "<< ifil+1 << ". ==>> "<< ifnames.at(ifil) <<endl;
    cout<<"requested events     = "<< params->Nevreq << endl;
    cout<<"initial random seeds = "<< params->rnd_ext << " " << params->rnd_int << endl;
    cout<<"photon cutoff k_0    = "<< params->k0cut << " GeV" << endl;
    cout<<"initial Wmax         = "<< params->Wmax <<endl;
    cout<<"................................................................"<< endl;
    cout<<"N generated events        = "<< params->MCsums.Nevgen << endl;
    cout<<"N weights                 = "<< params->MCsums.Nwgt << endl;
    cout<<"N negative weights        = "<< params->MCsums.Nwgt_Negative << endl;
    cout<<"N weights above Wmax      = "<< params->MCsums.Nwgt_OverMax << endl;
    cout<<"True Max weight           = "<< params->MCsums.WmaxTrue <<endl;
    cout<<"Cross section             = "<< params->MCsums.Xsec << " +/- " << params->MCsums.XsecErr << endl;
    cout<<"Cross section (negative)  = "<< params->MCsums.Xsec_Negative << " +/- " << params->MCsums.Xsec_Negative_Err << endl;
    cout<<"Cross section (above max) = "<< params->MCsums.Xsec_OverMax << " +/- " << params->MCsums.Xsec_OverMax_Err << endl;

    // incremental sums
    sumc.Nevgen         += params->MCsums.Nevgen;
    sumc.Nwgt           += params->MCsums.Nwgt;
    sumc.Nwgt_Negative  += params->MCsums.Nwgt_Negative;
    sumc.Swgt           += params->MCsums.Swgt;
    sumc.Swgt_Negative  += params->MCsums.Swgt_Negative;
    sumc.SQwgt          += params->MCsums.SQwgt;
    sumc.SQwgt_Negative += params->MCsums.SQwgt_Negative;
    sumc.Nwgt_OverMax   += params->MCsums.Nwgt_OverMax ;

    // define true Wmax (over all files)
    if (params->MCsums.WmaxTrue > WmaxFinal) WmaxFinal = params->MCsums.WmaxTrue;
  }

  // Final quantities over all the files
  Double_t Avwgt = sumc.Swgt / sumc.Nwgt;
  Double_t ErrAvwgt = sqrt( (sumc.SQwgt/sumc.Nwgt - Avwgt*Avwgt) /sumc.Nwgt );
  
  Double_t Avwgt_Negative(0); 
  Double_t ErrAvwgt_Negative(0); 
  if (sumc.Nwgt_Negative > 0) {
    Avwgt_Negative = sumc.Swgt_Negative / sumc.Nwgt;
    ErrAvwgt_Negative = sqrt( (sumc.SQwgt_Negative/sumc.Nwgt - Avwgt_Negative*Avwgt_Negative) /sumc.Nwgt );
  }

  Double_t sigma0 = pargen.Wnorm;
  sumc.Xsec = sigma0 * Avwgt;
  sumc.XsecErr = sigma0 * ErrAvwgt;
  sumc.Xsec_Negative = sigma0 * Avwgt_Negative;
  sumc.Xsec_Negative_Err = sigma0 * ErrAvwgt_Negative;

  // FINAL printouts
  cout<< endl;
  cout<<"================================================================"<< endl;
  cout<<"==== TOTAL STATISTICS =========================================="<< endl;
  cout<<"================================================================"<< endl;
  cout<<"N generated events        = "<< sumc.Nevgen << endl;
  cout<<"N weights                 = "<< sumc.Nwgt << endl;
  cout<<"N negative weights        = "<< sumc.Nwgt_Negative << endl;
  cout<<"N weights above Wmax      = "<< sumc.Nwgt_OverMax << endl;
  cout<<"True Max weight           = "<< WmaxFinal <<endl;
  cout<<"Cross section             = "<< sumc.Xsec << " +/- " << sumc.XsecErr << endl;
  cout<<"Cross section (negative)  = "<< sumc.Xsec_Negative << " +/- " << sumc.Xsec_Negative_Err << endl;
  //  cout<<"Cross section (above max) = "<< sumc.Xsec_OverMax << " +/- " << sumc.Xsec_OverMax_Err << endl;
  cout<<"Cross section (above max) =  *** TO BE DEFINED *** " << endl;
  cout<<"================================================================"<< endl;

  // Fast Simulation //
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%
  MuE::FastSim fs(pargen,fsi);  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // ANALYSIS //
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  MuE::Analysis analyzer(cfg, pargen, fsi, an_input);
  analyzer.BeginJob();
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ///////////////////////////////////////////////////////////////////////////////
  ////////////////////////////   EVENT  LOOP   //////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  //
  Long64_t n_events = cfg.n_events;
  
  if (n_events == 0) {
    n_events = nevtot;
    if (sumc.Nevgen != nevtot) 
      cerr << "***WARNING: total number of chained events ("<<nevtot
	   <<") is different from the expected number of generated events ("<<sumc.Nevgen<<")"<<endl;
  }

  if (n_events <0) cerr << "\n" << "Kinematic distributions will be loaded from existing results." << endl;
  else cerr << "Number of events to be analyzed = " << n_events << endl;

  if (_debug_) cerr << "event chain ... " << endl;

  // number of events with negligible weight (skipped)
  int zero_wgt_events = 0;    
    
GammaFunctionGenerator* gamma= new GammaFunctionGenerator;
ECALProperties *ecalprop= new ECALProperties();    
EMECALShowerParametrization *myparam = new EMECALShowerParametrization(ecalprop,{100.0,0.1},{1.0,0.1,100.0,1.0},1,1);
ECAL *TheEcal= new ECAL(5,-7.125,7.125,5,-7.125,7.125);    
    

  for (Long64_t iEvent=0; iEvent < n_events; ++iEvent) {
    Long64_t ientry = chain.LoadTree(iEvent);
    if (ientry <0) {
      cerr << "***ERROR in chaining MuE, event = "<<iEvent<<", ientry = "<<ientry<<endl;
      break;
    }
    if (iEvent % 1000000 == 0 ) cerr << "\n processing event : " << iEvent <<"\r"<< flush;
    
    chain.GetEntry(iEvent);
    if (_debug_ && iEvent<10)  chain.Show();

    Double_t evwgt = event->wgt_full;

    if (std::abs(evwgt) > 1e-17) {

      //%%%%%%%%%%%%%%%%%%%%%%%%%%%
      fs.Process(*event,gamma,myparam,TheEcal);
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%
      analyzer.Analyze(*event, fs);
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%
    }
    else {
      zero_wgt_events++;
      fs.RandomNrSync();
    }
  }
    
    
//questo lo devi fare alla fine di tutti gli eventi
  //TheEcal->Print_();
    
    
    
  cout<<endl<< "End reading MuE events. Read "<< n_events << " events." <<endl;
  cout<<"number of zero-weight events = "<< zero_wgt_events <<endl;

  cout<<endl<<"last Random seed = "<<gRandom->GetSeed()<<endl;
  cout<<"last call to Rndm() = "<< setw(20) << setprecision(17) <<gRandom->Rndm() <<endl;
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  analyzer.EndJob(sumc, mc_inputs);
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  delete event;
  delete params;

  cout << "Finished." << endl;

  if (!gROOT->IsBatch()) app->Run();

  return 0;
}
