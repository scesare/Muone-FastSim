#ifndef FastSim_H
#define FastSim_H

///////////////////////////////////////////////
// Fast Simulation of MuE scattering
//
// G.Abbiendi  4/Sep/2018 
///////////////////////////////////////////////
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/Vector3Dfwd.h"
#include "MuEtree.h"
#include "MuEana.h"
#include "Inputs.h"
#include "TMatrixD.h"
#include "TMatrixFBase.h"
#include <TMatrixFSym.h>
#include "TString.h"
#include "GammaFunctionGenerator.h"
#include "EMECALShowerParametrization.h"
#include "ECAL.h"

namespace MuE {

  class FastSim {

  public:
   FastSim(const MCpara & pargen, const FS_Input & fsi, bool _debug_=false);
    virtual ~FastSim(){};

    void Process(const Event & event,GammaFunctionGenerator* & gamma,EMECALShowerParametrization* const & myParam, ECAL* const & myGrid);

    const KineVars & GetGenKin() const {return genKin;}
    // const KineVars & GetDetKin() const {return detKin;}
    const KineVars & GetDetKinBeamRot() const {return detKinBeamRot;}
    const Photon & GetPhoton() const {return photon;}

    void RandomNrSync();

  private:
    typedef ROOT::Math::PxPyPzEVector PxPyPzEVector;
	typedef ROOT::Math::XYZVector XYZVector;
    Double_t P_2bodies_CoM(Double_t Mass, Double_t mass_mu, Double_t mass_e) const;

    PxPyPzEVector Lorentz_ToCoM(const PxPyPzEVector & plab) const;
    PxPyPzEVector Lorentz_ToLab(const PxPyPzEVector & pcm) const;

    // Double_t ThetaRMS(const PxPyPzEVector & p) const; 
    TMatrixD RotDivIN(const PxPyPzEVector & p) const; 
    PxPyPzEVector RotDiv(const PxPyPzEVector & p,const PxPyPzEVector & o) const; 
    //PxPyPzEVector Smear(const PxPyPzEVector & p) const; 
   // PxPyPzEVector SmearX(const PxPyPzEVector & p) const; 
   // PxPyPzEVector SmearPolar(const PxPyPzEVector & p) const; 
    //TMatrixD coo(const Double_t & a, const Double_t & s,const Double_t & ae, const Double_t & se) const; 
    //TMatrixD coo(const PxPyPzEVector & p,const PxPyPzEVector & q) const; 
    // TMatrixD MCSin(const PxPyPzEVector & k) const; 
    TMatrixD MCSout(const PxPyPzEVector & kin, const PxPyPzEVector & k, const PxPyPzEVector & ke) const; 
    TMatrixD MCSphoton(const PxPyPzEVector & p_gamma_Lab_div,const Double_t & xin,const Double_t & yin) const;
    TMatrixD Def_angle(const PxPyPzEVector & p_mu_in_div,const PxPyPzEVector & p_mu_out_div,const PxPyPzEVector & p_e_out_div) const;
    Double_t Def_angle_ph(const PxPyPzEVector & p_mu_in_div,const PxPyPzEVector & p_gamma_lab_div) const;
    Int_t ECALe(const Double_t & x,const Double_t & y) const;
    Int_t ECALph(const Double_t & x,const Double_t & y) const;
    void LoadKineVars(const PxPyPzEVector & p_mu_in,  const PxPyPzEVector & p_e_in, const PxPyPzEVector & p_mu_out, const PxPyPzEVector & p_e_out,  const TMatrixD & coo, const Double_t & TheINT,const Double_t & ThMuINT, KineVars & kv, ECAL* const & myGrid);
    void LoadPhoton(const Event & event, Photon & photon, const PxPyPzEVector & p_gamma_lab_div,const PxPyPzEVector & p_mu_in,const Double_t & x,const Double_t & y, ECAL* const & myGrid);
    void LoadECAL(KineVars & kv, ECAL* const & myGrid,int j);

    static const Double_t mm_PDG; // PDG muon mass 
    static const Double_t me_PDG; // PDG electron mass

    const Double_t & mm; // muon mass
    const Double_t & me; // electron mass
    const Double_t & Ebeam; // average beam energy
    const Double_t & EbeamRMS; // beam energy spread

    const Int_t & model; // model for detector smearing (gaussian resolution)
    const Int_t & MSopt; // options for multiple scattering (default/only Xplane/only polar)
    const Double_t & thickness; // material thickness (in X0) for model_=0
    const Double_t & intrinsic_resolution; // intrinsic resolution (in mrad) for model_=0

    Double_t sigSI; //sigma distribuzione angolare MCS da PDG nel silicio
    Double_t sigBE; //sigma distribuzione angolare MCS da PDG nel Berillio
    Double_t sigBE2in; //sigma distribuzione angolare MCS da PDG nel Berillio fino dove interagisce
   const Double_t sSin; //m spessore silicio
   const Double_t x0S; // m
   const Double_t sB; //m spessore berillio
   const Double_t x0B; // m
    Int_t tar; // target where mu interacts
    Double_t vertex; // where mu interacts in the target

      
    bool debug;

    PxPyPzEVector p_system; // mu-e centre-of-mass system fourmomentum
    Double_t Minv; // event invariant mass = sqrt(s)
    KineVars genKin; // kinematic variables at Gen-level for e and mu track
   // KineVars detKin; // kinematic variables at Detector-level for e and mu track
    KineVars detKinBeamRot; // kinematic variables at Detector-level for e and mu track with divergence
    Photon photon; // photon variables at Gen-level


  };
}

#endif
