//FAMOS headers
#include "GammaFunctionGenerator.h"
#include "IncGamma.h"

#include <TRandom3.h>

GammaFunctionGenerator::GammaFunctionGenerator() {
  xmax = 30.;
    


  for (unsigned i = 1; i <= 12; ++i) {
    // For all let's put the limit at 2.*(alpha-1)      (alpha-1 is the max of the dist)
    approxLimit.push_back(2 * ((double)i));
    //TF1* gamma=myIncompleteGamma.MyGamma();
    //myIncompleteGamma.Set_a((double)i,gamma);
    //myIncompleteGamma.a().setValue((double)i);
    integralToApproxLimit.push_back(myIncompleteGamma.MyGamma((double)i,approxLimit[i - 1]));
    theGammas.push_back(GammaNumericalGenerator((double)i, 1., 0, approxLimit[i - 1] + 1.));
  }
  coreCoeff.push_back(0.);  // alpha=1 not used
  coreCoeff.push_back(1. / 8.24659e-01);
  coreCoeff.push_back(1. / 7.55976e-01);
  coreCoeff.push_back(1. / 7.12570e-01);
  coreCoeff.push_back(1. / 6.79062e-01);
  coreCoeff.push_back(1. / 6.65496e-01);
  coreCoeff.push_back(1. / 6.48736e-01);
  coreCoeff.push_back(1. / 6.25185e-01);
  coreCoeff.push_back(1. / 6.09188e-01);
  coreCoeff.push_back(1. / 6.06221e-01);
  coreCoeff.push_back(1. / 6.05057e-01);

}

GammaFunctionGenerator::~GammaFunctionGenerator(){}

double GammaFunctionGenerator::shoot() const {
  if (alpha < 0.)
    return -1.;
  if (badRange)
    return xmin / beta;
  if (alpha < 12) {
    if (alpha == na) {
      return gammaInt() / beta;
    } else if (na == 0) {
      return gammaFrac() / beta;
    } else {
      double gi = gammaInt();
      double gf = gammaFrac();
      return (gi + gf) / beta;
    }
  } else {
    // an other generator has to be used in such a case
    return -1.;
  }
}

double GammaFunctionGenerator::gammaFrac() const {
  /* This is exercise 16 from Knuth; see page 135, and the solution is
     on page 551.  */

  double p, q, x, u, v;
  p = M_E / (frac + M_E);
  do {
    u = gRandom->Uniform(0.0,1.0);
    v = gRandom->Uniform(0.0,1.0);

    if (u < p) {
      x = exp((1 / frac) * log(v));
      q = exp(-x);
    } else {
      x = 1 - log(v);
      q = exp((frac - 1) * log(x));
    }
  } while (gRandom->Uniform(0.0,1.0) >= q);

  return x;
}

double GammaFunctionGenerator::gammaInt() const {
  // Exponential distribution : no approximation
  if (na == 1) {
    return xmin - log(gRandom->Uniform(0.0,1.0));
  }

  unsigned gn = na - 1;

  // are we sure to be in the tail
  if (coreProba == 0.)
    return xmin - coreCoeff[gn] * log(gRandom->Uniform(0.0,1.0));

  // core-tail interval
  if (gRandom->Uniform(0.0,1.0) < coreProba) {
    return theGammas[gn].gamma_lin();
  }
  //  std::cout << " Tail called " << std::endl;
  return approxLimit[gn] - coreCoeff[gn] * log(gRandom->Uniform(0.0,1.0));
}

void GammaFunctionGenerator::setParameters(double a, double b, double xm) {

std::cout << "Setting parameters " << std::endl;
alpha = a;
std::cout << "Setting parameters " << std::endl;
beta = b;
xmin = xm * beta;
    
  if (xm > xmax) {
    badRange = true;
    return;
  }
  badRange = false;
  na = 0;
   
  if (alpha > 0. && alpha < 12)
    na = (unsigned)floor(alpha);

  frac = alpha - na;
  // Now calculate the probability to shoot between approxLimit and xmax
  // The Incomplete gamma is normalized to 1
  if (na <= 1)
    return;

  //myIncompleteGamma.a().setValue((double)na);
  //TF1* gamma=myIncompleteGamma.MyGamma();
//myIncompleteGamma.Set_a((double)na,gamma);

  unsigned gn = na - 1;
  //  std::cout << " na " << na << " xm " << xm << " beta " << beta << " xmin " << xmin << " approxLimit " << approxLimit[gn] << std::endl;
  if (xmin > approxLimit[gn]) {
    coreProba = 0.;
  } else {
    double tmp = (xmin != 0.) ? myIncompleteGamma.MyGamma((double)na,xmin) : 0.;
    coreProba = (integralToApproxLimit[gn] - tmp) / (1. - tmp);
    theGammas[gn].setSubInterval(xmin, approxLimit[gn]);
  }
    
  //  std::cout << " Proba " << coreProba << std::endl;
}