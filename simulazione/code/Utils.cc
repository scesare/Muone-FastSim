#include <string>
#include <sstream>
#include <iostream>
#include "Utils.h"

using namespace std;

namespace MuE {

istringstream input_line (ifstream & input_file, bool debug) {
  string line;
  getline(input_file, line);
  istringstream stream(line);
  if (debug) cout << "\t" << stream.str() << endl;
  return stream;
}

bool CheckParameters(const MCpara & pargen_0, const MCpara & pargen) 
{
  bool checkOk = true;
  
  if (pargen.program_version != pargen_0.program_version) {
    cerr << "***ERROR: inconsistent input files. program_version = "
	 <<pargen_0.program_version<<" (ref), "<<pargen.program_version<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.UNWGT != pargen_0.UNWGT) {
    cerr << "***ERROR: inconsistent input files. UNWGT = "
	 <<pargen_0.UNWGT<<" (ref), "<<pargen.UNWGT<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.Mode != pargen_0.Mode) {
    cerr << "***ERROR: inconsistent input files. Mode = "
	 <<pargen_0.Mode<<" (ref), "<<pargen.Mode<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.rnd_ext == pargen_0.rnd_ext) {
    cerr << "***ERROR: input files have the same random seed. rnd_ext, rnd_int = "
	 <<pargen_0.rnd_ext<<", "<<pargen_0.rnd_int <<endl;
    checkOk = false;
  }
  if (pargen.Ebeam != pargen_0.Ebeam) {
    cerr << "***ERROR: inconsistent input files. Ebeam = "
	 <<pargen_0.Ebeam<<" (ref), "<<pargen.Ebeam<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.EbeamRMS != pargen_0.EbeamRMS) {
    cerr << "***ERROR: inconsistent input files. EbeamRMS = "
	 <<pargen_0.EbeamRMS<<" (ref), "<<pargen.EbeamRMS<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.charge_mu != pargen_0.charge_mu) {
    cerr << "***ERROR: inconsistent input files. charge_mu = "
	 <<pargen_0.charge_mu<<" (ref), "<<pargen.charge_mu<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.mass_mu != pargen_0.mass_mu) {
    cerr << "***ERROR: inconsistent input files. mass_mu = "
	 <<pargen_0.mass_mu<<" (ref), "<<pargen.mass_mu<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.mass_e != pargen_0.mass_e) {
    cerr << "***ERROR: inconsistent input files. mass_e = "
	 <<pargen_0.mass_e<<" (ref), "<<pargen.mass_e<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.invalfa0 != pargen_0.invalfa0) {
    cerr << "***ERROR: inconsistent input files. invalfa0 = "
	 <<pargen_0.invalfa0<<" (ref), "<<pargen.invalfa0<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.k0cut != pargen_0.k0cut) {
    cerr << "***Warning: input files have different k0cut: "
	 << pargen_0.k0cut << " (ref), "<< pargen.k0cut << " (new)" <<endl;
    cerr << "\t this may not a problem, k0cut is just a technical parameter, but it has to be small enough not to affect the physics results."<<endl;
  }
  if (pargen.Emin_e != pargen_0.Emin_e) {
    cerr << "***ERROR: inconsistent input files. Emin_e = "
	 <<pargen_0.Emin_e<<" (ref), "<<pargen.Emin_e<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.Wnorm != pargen_0.Wnorm) {
    cerr << "***ERROR: inconsistent input files. Wnorm = "
	 <<pargen_0.Wnorm<<" (ref), "<<pargen.Wnorm<<" (new)" <<endl;
    checkOk = false;
  }
  if (pargen.Wmax != pargen_0.Wmax) {
    cerr << "***Warning: input files have different Wmax: "
	 << pargen_0.Wmax << " (ref), "<< pargen.Wmax << " (new)" <<endl;
    cerr << "\t this may not a problem, Wmax is just a technical parameter, but it has to be greater than the true Wmax known at the end."<<endl;
  }
  
  return checkOk;
}

}
