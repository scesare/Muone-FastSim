#include <iostream>
#include <cmath>
#include <cstdio>
#include "dalpha.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////
// FORTRAN Function giving the contribution of lepton or top quark to Dalpha
// mass : mass of the lepton or top (in GeV)
// q2 : invariant momentum transfer (in GeV^2), negative for space-like
// i : 1=electron; 2=muon; 3=tau; 4=TOP quark
extern "C" double summa_(double *mass, double *q2, int *i);
////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FORTRAN subroutine (Fred Jegerlehner's Dalpha_had as a function of the scale)
// de      : energy scale (in GeV) // spacelike region : de = -sqrt(-t)
// dst2    : sin^2(theta_ew)
// dder    : Dalpha_had
// derrder : uncertainty on Dalpha_had
// ddeg and derrdeg regard the weak SU2 coupling
extern "C" void dhadr5n12_(double *de, double *dst2, double *dder, double *derrder, double *ddeg, double *derrdeg);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
// Leptonic Delta(alpha) from C.Carloni
//
double dalep(double t) {
  const double alfa0 = 1./137.03599907430637;
  //  const double pi = acos(-1.);
  const double pi = 3.1415926535897932384626433832795029;
  const double alpi = alfa0/pi;
  // masses: electron, muon, tau, top quark
  //  used in Carlo's code (for reweighting)
  double mass[4] = {0.510998928e-3,0.1056583715,1.77682,175.6};

  double Sum = 0.;
  for (int i=0; i<4; ++i) {
    int ifla = i+1;
    Sum += summa_(&mass[i], &t, &ifla);
  }
  return alpi*Sum;
}
////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Interface to Fred Jegerlehner's Dalpha_had(t)
// t = space-like momentum transfer (t<0)
//
double dahadFred(double t) {
  double tn = -t;
  double Q = -sqrt(tn);
  double st2 = 0.2322;
  double der(0.), errder(0.), deg(0.), errdeg(0.);
  dhadr5n12_(&Q, &st2, &der, &errder, &deg, &errdeg);
  return der;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parametrizations of the Delta(alpha) hadronic
// iparam = 0 : pol2 
// iparam = 1 : Lepton-Like (Carloni's original)
// iparam = 2 : modified Lepton-Like, (par[0])_LLmod = (par[0]/par[1])_LL 
//
double dahadPar(int iparam, double t, double par[2]) {

  double dahad(0);

  if (iparam==0) {
    dahad = par[0]*t+par[1]*t*t;
  }
  else if (iparam==2) {
    double squ = sqrt(1.-4.*par[1]/t);
    dahad = par[0]*par[1]*
      (-5./9.-4./3.*par[1]/t + (4./3.*par[1]*par[1]/t/t + par[1]/3./t -1./6.)*2./squ*log(std::abs((1.-squ)/(1.+squ))));
  }
  else if (iparam==1) {
    double squ = sqrt(1.-4.*par[1]/t);
    dahad = par[0] *
      (-5./9.-4./3.*par[1]/t + (4./3.*par[1]*par[1]/t/t + par[1]/3./t -1./6.)*2./squ*log(std::abs((1.-squ)/(1.+squ))));
  }
  else 
    {cout<<"\n ***ERROR: unknown dahad parametrization: iparam = "<< iparam <<endl;}

  return dahad;
}
