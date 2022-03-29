#ifndef Analysis_H
#define Analysis_H

//////////////////////////////////////////////////////
// class for Analysis of NLO MC MuE scattering events
//
// G.Abbiendi  4/Dec/2018 
//////////////////////////////////////////////////////

#include "TFile.h"
#include "TTree.h"
#include "MuEtree.h"
#include "MuEana.h"
#include "FastSim.h"
//#include "Histos.h"
#include "Inputs.h"

namespace MuE {
  
  class Analysis {
    
  public:
    Analysis(const One_Input & p, const MCpara & g, const FS_Input & fsi, const AN_Input & a)
      : parmain(p), pargen(g), parsim(fsi), paran(a) {}

    virtual ~Analysis(){};

    void BeginJob();
    void Analyze(const Event & event, const FastSim & fs);
    void EndJob(const MCstat & mcsums);
    void EndJob(const MCstat & mcsums, const std::vector<std::pair<std::string,Long64_t> > & mc_inputs);

  private:
    const One_Input & parmain; 
    const MCpara & pargen;
    const FS_Input & parsim; 
    const AN_Input & paran;

    MuEana myAna;
    TFile * output_hist_file;
    TFile * output_tree_file;
    TTree * atree;
    
   // Histos * histos;

  };
}

#endif
