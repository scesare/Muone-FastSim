#include <string>
#include <iostream>
#include "Analysis.h"

using namespace MuE;
using namespace std;

void Analysis::BeginJob()
{
  output_hist_file = new TFile("results.root","RECREATE");
  
  if (paran.makeTree) 
    {
      output_tree_file = new TFile("outtree.root","RECREATE");
      atree = new TTree("atree","Analysis output tree");
      Int_t splitBranches = 2;
      atree->Branch("event",&myAna,64000,splitBranches);
    }
  
  output_hist_file->cd();
  
  cout<<"\n"<<"Analysis Inputs: thetaMax for preselection and histograms = "<< paran.thetaMax << " mrad"<<endl;

 // histos = new Histos(pargen, paran);
}

void Analysis::Analyze(const MuE::Event & event, const MuE::FastSim & fs)
{
  //const MuE::KineVars & genKin = fs.GetGenKin();
  //const MuE::KineVars & detKin = fs.GetDetKin();
  const MuE::KineVars & detKinBeamRot = fs.GetDetKinBeamRot();
    
  const MuE::Photon & photon = fs.GetPhoton();

  // apply preselection if both the gen-level and det-level electron angle are above the cut (default 30mrad)&&
  //if (genKin.the > paran.thetaMax && detKin.the > paran.thetaMax ) return;
    if (detKinBeamRot.def_angle_e > paran.thetaMax) return;
    //if (detKinBeamRot.tar==0 && detKinBeamRot.def_angle_e >35 ) return;
    //if (detKinBeamRot.tar==1 && detKinBeamRot.def_angle_e >70 ) return;

    
    

  // filling my analysis variables
  myAna.RunNr = event.RunNr;
  myAna.EventNr = event.EventNr;
  myAna.wgt_full = event.wgt_full;
  myAna.wgt_norun = event.wgt_norun;
  myAna.wgt_lep = event.wgt_lep;
  myAna.wgt_LO = event.wgt_LO;
  myAna.E_mu_in = event.E_mu_in;      
  //myAna.genKin = genKin;
  //myAna.detKin = detKin;
  myAna.detKinBeamRot = detKinBeamRot;
  myAna.photon = photon;
  
  // filling my analysis tree
  if (paran.makeTree) atree->Fill();
  
  // filling my histos
  //histos->Fill(event, myAna);
}

void Analysis::EndJob(const MCstat & mcsums) 
{
  Long64_t n_events = mcsums.Nevgen;

  if (n_events >=0) {
    /* cerr << "Final histogram normalizations, ratios, fits, plots "<<endl<<endl;
    histos->SetSums(mcsums);
    histos->Normalize(n_events);

    histos->Do_Ratios();
    histos->Plot();
    if (paran.doTemplates) histos->Plot2D();
    histos->Fit();
    histos->PlotResolutions();*/

    // write out root tree
    if (paran.makeTree) {
      atree->Print();
      output_tree_file->Write();
    }
    
    // write out hist root file
    //output_hist_file->Write();
    
    if (n_events ==0) 
      cerr <<endl<< "Processed all input events, histograms written to file: results.root" << endl;
    else
      cerr <<endl<< "Processed " << n_events << " events, histograms written to file: results.root" << endl;
  }

  else {
    cerr <<"\n"<< "*** WARNING: Analysis::EndJob missing argument mc_inputs with n_events < 0 \n"<<endl;
    exit(999);
  }

}

void  Analysis::EndJob(const MCstat & mcsums, const std::vector<std::pair<std::string,Long64_t> > & mc_inputs)
{
  Long64_t n_events = parmain.n_events;

  if (n_events >=0) {
  /*  cerr <<endl<< "Final histogram normalizations, ratios, fits, plots "<<endl<<endl;
    histos->SetSums(mcsums);
    histos->Normalize(n_events);

    histos->Do_Ratios();
    histos->Plot();
    if (paran.doTemplates) histos->Plot2D();
    histos->Fit();
    histos->PlotResolutions();
*/
    // write out root tree
    if (paran.makeTree) {
      atree->Print();
      output_tree_file->Write();
    }
    
    // write out hist root file
    output_hist_file->Write();
    
 /*   if (n_events ==0) 
      cerr <<endl<< "Processed all input events, histograms written to file: results.root" << endl;
    else
      cerr <<endl<< "Processed " << n_events << " events, histograms written to file: results.root" << endl;
  }

  else {
    // load histos from existing external file
    TFile *fp = new TFile(parmain.histo_ifname.c_str());

    //    output_hist_file->cd();
    
    TFile *projfile = new TFile("projstat.root","RECREATE");
    projfile->cd();
    histos->RatioFinal(n_events, fp, "hn_thmu");
    histos->RatioFinal(n_events, fp, "hn_the");
    projfile->Write();

    cerr <<"\n"<< "*** WARNING: missing code to normalize histograms to be used in analysis !!! \n"<<endl;

    // test numerical precision comparing results from FORTRAN code to ROOT
    //histos->LoadExtHistos(fp);
    //histos->LoadCarlos(mc_inputs);
    //histos->CompareWithCarlos();
  }*/
  
}}

