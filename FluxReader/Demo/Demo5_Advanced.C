// This fifth and final demo introduces many of the other functions
// provided by the FluxReader framework
// It is certainly not exhausted, but should give a good overview
// It discusses creating entirely new Parents and Detectors,
// how to reuse Dk2Nu entries at detectors with smearing,
// how to set the number and subset of files to run over,
// and discusses cross sections a bit more

#ifdef __CINT__
void Demo5_Advanced()
{
  std::cout << "Sorry, you must run in compiled mode." << std::endl;
}
#else

// C/C++ Includes
#include <iostream>
#include <string>

// ROOT Includes
#include "TFile.h"

// Package Includes
#include "Combiner.h"
#include "Detectors.h"
#include "FluxReader.h"
#include "Parameters.h"
#include "ParticleParam.h"
#include "Utilities.h"
#include "Vars.h"

using namespace flxrd;

void Demo5_Advanced()
{
  Parameters p(false);

  p.AddDetector(kNOvA_ND);
  // p.AddDetector(kNOvA_FD);

  // The Parent object is used to define a neutrino... parent
  // This is just a name (std::string) and pdg (int)
  // Thus, the user can generate custom ones
  // Parent kKShort("KShort", 310);
  // p.AddParent(kKShort);

  // The Detector object defines... a detector
  // It consists of a name (std::string), dominant nuclear target (std::string),
  // a position in space (cm in detector coordinates, 3 doubles),
  // A size (cm, 3 doubles), and a "number of uses" (int)
  // The number of uses defines how many times to smear
  // a neutrino ray through the detector
  /*Detector kNOvA_ND_Shift("NOvA-ND-Shift", "CH2",
                           1141.4, -345.6, 99566.5,     // position, z shifted by +100 cm
                           262.14, 393.27, 1424.52698,  // size
                           1); // uses */
  // p.AddDetector(kNOvA_ND_Shift);

  string dk2nu_loc = "/nusoft/data/flux/dk2nu/nova/2010/flugg_mn000z200i_20101117.gpcfgrid_lowth/";
  dk2nu_loc += "*dk2nu.root";
  FluxReader *fr = new FluxReader(dk2nu_loc, 2);

  // The FluxReader Constructor generates the same expanded list of files,
  // everytime it expands a wildcard (assuming there are no name changes, additions, or removals)
  // The user can set the number of files to use, and the number of files to skip
  // Default: use all files
  // FluxReader *fr = new FluxReader(dk2nu_loc);
  // Second input: set the number of files to use (100 in this example)
  // FluxReader *fr = new FluxReader(dk2nu_loc, 100);
  // Third input: how many files to skip (300 in this example)
  // FluxReader *fr = new FluxReader(dk2nu_loc, 100, 300);
  // To use all files after skipping some, leave the second input as 0 (in this example, all files after the first 200):
  // FluxReader *fr = new FluxReader(dk2nu_loc, 0, 200);
  // Note that the code will quit if there are no files after construction
  // Furthermore, if there are less files than what is specified by the second input,
  // the code will run on what is left but output the (smaller) number of files
  // This condition can occur if the user specifies too many files to begin with,
  // or if too many files get skipped

  // At the NOvA ND, it could make sense to smear neutrino rays throughout the volume
  // This functionality exists in the Parameters object with Parameters::SetDetUses
  // The first argument is the detector name (std::string),
  // and the second is the number of times to smear the neutrino rays
  // As a note: each use of the neutrino ray picks a random point inside the detector,
  // and these points are different for each neutrino ray
  // (I.e., the same "random" points are not used for each neutrino ray)
  // It is NOT possible to run the same Detector in different Spectra with different uses--
  // only one of the values will be used
  // It is possible to predict which uses value will be the one picked,
  // but it is safer to set one value and stick with it,
  // especially since this is the only one that will be used anyway
  // p.SetDetUses("NOvA-ND", 10);

  // Add a Spectra object
  fr->AddSpectra(p, "enu", "Energy (GeV)", Bins(100, 0., 10.), kEnergy);

  TFile* out = new TFile("/nova/ana/users/gkafka/FluxReader/demo5.root", "RECREATE");

  fr->ReadFlux(out);
  out->Close();
  delete fr;

  // The cross section class, XSec, can generate cross section plots
  // It gets this information from a file found using the $GENIEXSECPATH environment variable
  // The class is constructed with no arguments
  // XSec* xsec = new XSec();
  // It can generate plots as a TGraph*, TSpline3*, or TH1*, with the functions
  // GetGraph, GetXSec, and GetHist
  // GetGraph and GetXSec take the same inputs:
  // a neutrino PDG (SIGNED int), a nuclear target (std::string), and current (std::string)
  // TGraph* g   = xsec->GetGraph(+14, "CH2", "tot_cc");
  // TSpline3* s = xsec->GetXSec( +14, "CH2", "tot_cc");
  // GetHist takes a TSpline3*, number of bins (int), and either:
  // a minimum and maximum edge (2 doubles) for equally sized bins, or
  // a double* of bin edges (really, a pointer to an array of bin edges)
  // TH1* h      = xsec->GetHist(s, 120, 0., 120.); // Note s from above!
  // Cross section ratios can be generated as well
  // GetGraphRatio and GetXSecRatio take two sets of inputs,
  // the first as numerator, the second as denominator
  // TSpline3* r = xsec->GetXSecRatio(+14, "CH2", "tot_cc", +14, "CH2", "tot_nc");
  // Since GetHist takes a TSpline3* as an input anyway,
  // to generate a ratio as a TH1* does not require a new function
  // Supported neutrino PDG are +/-12, 14 and 16
  // Supported nuclei are H, C, N, O, S, Cl, Ti, Fe, CH2
  // For a list of currents, uncomment the following line
  // xsec->ListIntTypes();
  // delete xsec;
}

#endif
