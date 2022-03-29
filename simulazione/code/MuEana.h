#ifndef MuEana_H
#define MuEana_H
#include <vector>
///////////////////////////////////////////////
// Classes defining MuE analysis variables
//
// G.Abbiendi  4/Dec/2018 
///////////////////////////////////////////////

namespace MuE {

  // Analysis variables in the output tree
  
  class KineVars {

  public:
    Double_t t13; // Mandelstam t (muon leg)
    Double_t t24; // Mandelstam t (electron leg)
    Double_t x13; // Feynman x (muon leg)
    Double_t x24; // Feynman x (electron leg)
    Double_t tt_e; // t computed from electron angle with LO formulas
    Double_t xt_e; // x computed from electron angle with LO formulas
    Double_t x_in;
    Double_t y_in;
    Double_t Ee; // electron energy
    Double_t Emu; // muon energy
    Double_t the; // electron theta (in mrad)
    Double_t thmu; // muon theta (in mrad)
    Double_t phe; // electron phi (from -pi to +pi)
    Double_t phmu; // muon phi (from -pi to +pi) 

    Double_t deltaPhi; // acoplanarity (deltaPhi)
    Double_t openingAngle; // opening angle mu-e out in the Lab
    Double_t tripleProduct; // triple product btw normalized vectors i . mu x e
    Double_t cooXe;
    Double_t cooXmu;
    Double_t cooYe;
    Double_t cooYmu;
    Double_t pXmu;
/*Double_t pYmu;
Double_t pZmu;
Double_t pXe;
Double_t pYe;
Double_t pZe;
Double_t pXmu_out;
Double_t pYmu_out;
Double_t pZmu_out;*/
Double_t pXe_out;
Double_t pYe_out;
Double_t pZe_out;
Double_t Pmu_out;
Double_t Pe_out;
Double_t tar;
/*Double_t ThEl_interaction;
Double_t ThMu_interaction;*/
Double_t def_angle_mu;
Double_t def_angle_e; 
/*Double_t n_max_Cell;
Double_t E_clus3x3;
Double_t E_1;*/
Double_t Ecell1;
Double_t Ecell2;
Double_t Ecell3;
Double_t Ecell4;
Double_t Ecell5;
Double_t Ecell6;
Double_t Ecell7;
Double_t Ecell8;
Double_t Ecell9;
Double_t Ecell10;
Double_t Ecell11;
Double_t Ecell12;
Double_t Ecell13;
Double_t Ecell14;
Double_t Ecell15;
Double_t Ecell16;
Double_t Ecell17;
Double_t Ecell18;
Double_t Ecell19;
Double_t Ecell20;
Double_t Ecell21;
Double_t Ecell22;
Double_t Ecell23;
Double_t Ecell24;
Double_t Ecell25;
Double_t n_cell_e;

    KineVars():
      t13(0),t24(0),x13(0),x24(0),tt_e(0),xt_e(0),
      x_in(0),y_in(0),Ee(0),Emu(0),the(0),thmu(0),phe(0),phmu(0),deltaPhi(0),openingAngle(0),tripleProduct(0),cooXe(0),cooXmu(0),cooYe(0),cooYmu(0)/*,pXmu(0),pYmu(0),pZmu(0),pXe(0),pYe(0),pZe(0),pXmu_out(0),pYmu_out(0),pZmu_out(0)*/,pXe_out(0),pYe_out(0),pZe_out(0),Pmu_out(0),Pe_out(0),tar(-1)/*,ThEl_interaction(0),ThMu_interaction(0)*/,def_angle_mu(0),def_angle_e(0)/*,n_max_Cell(0),E_clus3x3(0),E_1(0)*/,Ecell1(0),Ecell2(0),Ecell3(0),Ecell4(0),Ecell5(0),Ecell6(0),Ecell7(0),Ecell8(0),Ecell9(0),Ecell10(0),Ecell11(0),Ecell12(0),Ecell13(0),Ecell14(0),Ecell15(0),Ecell16(0),Ecell17(0),Ecell18(0),Ecell19(0),Ecell20(0),Ecell21(0),Ecell22(0),Ecell23(0),Ecell24(0),Ecell25(0),n_cell_e(0)
    {};
    

      

    virtual ~KineVars(){};

    ClassDef(KineVars,1)
  };

  class Photon {

  public:
    Double_t energy;    // photon energy in the Lab frame
    Double_t theta;     //   "    theta in the Lab frame (in mrad)
    Double_t phi;       //   "    phi in the Lab frame (in rad)
    Double_t energyCoM; // photon energy in the Centre-of-Mass frame
    Double_t def_angle_ph;//   "    theta in the mu_in frame (in mrad)
   Double_t coox;
   Double_t cooy;
    Double_t n_cell_ph;

    
  Photon():
    energy(-1),theta(-1),phi(0),energyCoM(-1),def_angle_ph(-1),coox(-1),cooy(-1),n_cell_ph(0)
      {};
virtual ~Photon(){};

    ClassDef(Photon,1)
  };

  class MuEana {

   public:
    UInt_t RunNr;
    Long64_t EventNr;
    Double_t wgt_full, wgt_norun, wgt_lep, wgt_LO;   // event weights 
    Double_t E_mu_in;  // incoming muon energy
    KineVars genKin;   // kinematic variables at Generator-level for e and mu tracks
   // KineVars detKin; // kinematic variables at Detector-level for e and mu tracks
    KineVars detKinBeamRot; 
     Photon photon;     // photon kinematic variables at Gen-level
    
    MuEana():
     RunNr(0),EventNr(0),wgt_full(0),wgt_norun(0),wgt_lep(0),wgt_LO(0),E_mu_in(0)
       {};
    
     virtual ~MuEana(){};   

     ClassDef(MuEana,1)
  };
}
  

#endif
