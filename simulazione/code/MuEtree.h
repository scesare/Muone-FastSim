#ifndef MuEtree_H
#define MuEtree_H

///////////////////////////////////////////////
// Data Formats for MuE MC events
//
// G.Abbiendi  5/Jul/2018 
///////////////////////////////////////////////

#include <vector>
#include <string>
#include <fstream>
#include "TROOT.h"

namespace MuE {

  class Setup;

  class MCpara {
  public:
    std::string program_version;
    Int_t process_ID;
    std::string running_on;
    std::string start_time;
    Long64_t Nevreq;
    Bool_t UNWGT;
    std::string Mode;
    UInt_t rnd_ext, rnd_int;
    Double_t Ebeam, EbeamRMS;
    Double_t charge_mu;
    Double_t mass_mu, mass_e;
    Double_t invalfa0;
    Double_t k0cut; 
    Double_t Emin_e; 
    Double_t Wnorm; 
    Double_t Wmax;

    MCpara(std::ifstream &);
    MCpara(const Setup & s);    

    MCpara():
    process_ID(0),Nevreq(0),UNWGT(false),rnd_ext(0),rnd_int(0),
    Ebeam(0),EbeamRMS(0),charge_mu(0),mass_mu(0),mass_e(0),invalfa0(0),k0cut(0),Emin_e(0),Wnorm(0),Wmax(0)
    {};

    virtual ~MCpara(){};
    
    ClassDef(MCpara,1)
  };

  class MCstat {

  public:
    Long64_t Nevgen;
    Long64_t Nwgt, Nwgt_Negative;
    Double_t Swgt, Swgt_Negative;
    Double_t SQwgt, SQwgt_Negative;
    Long64_t Nwgt_OverMax;
    Double_t WmaxTrue;
    Double_t Xsec, XsecErr;
    Double_t Xsec_Negative, Xsec_Negative_Err;
    Double_t Xsec_OverMax, Xsec_OverMax_Err;
   
    MCstat(std::ifstream &);

    MCstat():
    Nevgen(0), Nwgt(0), Nwgt_Negative(0),
    Swgt(0), Swgt_Negative(0), 
    SQwgt(0), SQwgt_Negative(0),
    Nwgt_OverMax(0),
    WmaxTrue(0),
    Xsec(0), XsecErr(0),
    Xsec_Negative(0), Xsec_Negative_Err(0),
    Xsec_OverMax(0), Xsec_OverMax_Err(0)
    {};

    virtual ~MCstat(){};
    
    ClassDef(MCstat,2)
  };

  class Setup {

  public:
    std::string program_version;
    Int_t process_ID;
    std::string running_on;
    std::string start_time;
    Long64_t Nevreq;
    Bool_t UNWGT;
    std::string Mode;
    UInt_t rnd_ext, rnd_int;
    Double_t Ebeam, EbeamRMS;
    Double_t charge_mu;
    Double_t mass_mu, mass_e;
    Double_t invalfa0;
    Double_t k0cut; 
    Double_t Emin_e; 
    Double_t Wnorm; 
    Double_t Wmax;

    MCstat MCsums;

    Setup(const MCpara & m, const MCstat & s);

    Setup():
    process_ID(0),Nevreq(0),UNWGT(false),rnd_ext(0),rnd_int(0),
    Ebeam(0),EbeamRMS(0),charge_mu(0),mass_mu(0),mass_e(0),invalfa0(0),k0cut(0),Emin_e(0),Wnorm(0),Wmax(0)
    {};

    inline MCpara GetMCpara() const {return MCpara(*this);}

    virtual ~Setup(){};

    ClassDef(Setup,5)
  };

  class P4 {

  public:
    Double_t E;
    Double_t px;
    Double_t py;
    Double_t pz;
    
    P4():
    E(0),px(0),py(0),pz(0)
    {};

    virtual ~P4(){};
    
    ClassDef(P4,2)
  };

  class Event {

  public:    
    UInt_t RunNr;
    Long64_t EventNr;
    Double_t wgt_full, wgt_norun, wgt_lep, wgt_LO;
    Double_t E_mu_in; 
    P4 P_mu_out;
    P4 P_e_out;
    std::vector<P4> photons;
    
    Event():
    RunNr(0),EventNr(0),wgt_full(0),wgt_norun(0),wgt_lep(0),wgt_LO(0),E_mu_in(0)
    {};

    bool Read(std::ifstream &);

    virtual ~Event(){};
    
    ClassDef(Event,5)
  };

}

#endif
