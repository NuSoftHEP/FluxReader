// This third demo introduces the Var and Weight classes
// It shows how to create a new Var object
// It next shows a different Weight object
// Finally, it demonstrates using weights external to the framework

#ifdef __CINT__
void Demo2_VarWeight()
{
  std::cout << "Sorry, you must run in compiled mode." << std::endl;
}
#else

// C/C++ Includes
#include <iostream>
#include <string>

// ROOT Includes
#include "TFile.h"
#include "TSpline.h"

// Package Includes
#include "Detectors.h"
#include "FluxReader.h"
#include "Parameters.h"
#include "Utilities.h"
#include "Vars.h"
#include "Weight.h"
#include "XSec.h"

using namespace flxrd;

void Demo2_VarWeight()
{
  Parameters p(false);

  // Add a couple of detectors
  p.AddDetector(knova_nd);
  p.AddDetector(knova_fd);

  FluxReader *fr = new FluxReader("/nova/ana/users/rhatcher/dk2nu-data/fullsplit/generic*.root");

  // Add a Spectra object
  fr->AddSpectra(p, "enu", "Energy (GeV)", Bins(100, 0., 10.), kEnergy);

  // kEnergy above is a Var object, and it is defined in Vars.h
  // Its definition is reproduced here with notes found inside parentheses
  // (a)const (b)Var (c)kEnergy((d){"nuray", "nuray.E"},
  //                  (e)[](const bsim::Dk2Nu* nu, const int& (*)i_nuray)
  //                  (f){ return nu->nuray[i_nuray].E; });
  // (a) The const label means that the Var will not change once it is created
  //     i.e., every entry will be calculated the same way
  // (b) This is just specifying the type of object, like int, or double
  // (c) This is the name of the object, like 'x' in "int x = 3;"
  // (d) This is the list of branch names needed for the variable
  // (e) This slightly opaque line is the type of function that will be used.
  //     to calculate the value for a given entry
  //     (*)If the Var does not access the "nuray" branch, i_nuray can be omitted
  // (f) The block of code inside the curly braces ("{", "}")
  //     determines what is returned, and how it is calculated
  //     kEnergy is simple enough to have this done in one line,
  //     but other more complicated variables can span multiple lines
  //     Simply separate lines by semicolons as in regular code

  // Let's make a new variable, of the absolute value of neutrino flavor, and call it kAbsFlavor
  // % const Var kAbsFlavor...
  // Looking at Dk2Nu documentation, this is the decay.ntype branch
  // % const Var kAbsFlavor({"decay", "decay.ntype"},...
  // Since this does not access the nuray branch, we can omit i_nuray
  // % const Var kAbsFlavor({"decay", "decay.ntype"},
  // %                      [](const bsim::Dk2Nu* nu, const int&)...
  // We could probably return the quantity of interest in one line, but let's be more explicit
  const Var kAbsFlavor({"decay", "decay.ntype"},
                       [](const bsim::Dk2Nu* nu, const int&)
                       { int nuflav = nu->decay.ntype;
                         return abs(nuflav); });

  // Now let's plot this variable
  // The binning will have to be different than before, since all entries should be 12 or 14
  fr->AddSpectra(p, "nuflav", "PDG", Bins(20, 0., 20.), kAbsFlavor);
  // This still corresponds to an event rate of sorts...

  // These won't be needed in the rest of the script
  p.RemoveXSec("tot_cc");
  p.RemoveXSec("tot_nc");

  // We could check how often the simulation has various mesons decay
  // to muon or electron neutrinos, but we would want to weight each neutrino as 1
  // There is a default Weight, kDefaultW, and also an object that applies no weight,
  // kNoWeight, both defined in Weight.h
  // There is also a file Weights.h for other commonly used weights
  // Adding a weight to a Spectra is done in the FluxReader::AddSpectra function,
  // as an optional argument after the Var
  fr->AddSpectra(p, "dkrate", "PDG", Bins(20, 0., 20.), kAbsFlavor, kNoWeight);
  // This Spectra will show the number of neutrinos of a given flavor, from a given parent

  // The Weight object is defined nearly identically to the Var object
  // The main difference comes in the 'type of function' line, note (e) above
  // For a Weight, this line becomes
  // (const double& w, const bsim::Dk2Nu* nu, const int& i_nuray, const TObject* extW)
  // where w is a default weight, and extW is external weights

  // Let's try applying a cross section ourselves
  // (Cross sections are described in more detail in a later demo)
  // We will assume that the cross section will be in a TSpline3* object, which is continuous
  // We need only the neutrino energy, found in nuray.E,
  // which of course accesses the nuray branch
  // % const Weight kAppXSec({"nuray", "nuray.E"},
  // %               [](const double& w, const bsim::Dk2Nu* nu, const int& i_nuray, const TObject* extW)...
  // Now we have to write the return function
  // It is important that we tell the function what kind of class the external weights are in,
  // otherwise we can ONLY access TObject functions, and we need the TSpline::Eval function
  // % const Weight kAppXSec({"nuray", "nuray.E"},
  // %               [](const double& w, const bsim::Dk2Nu* nu, const int& i_nuray, const TObject* extW)
  // %               { TSpline3* xsecspline = (TSpline3*)extW; ...
  // Now, we simply pull out the energy, evaluate the cross section here,
  // and return the cross section times the default weight, w
  const Weight kAppXSec({"nuray", "nuray.E"},
                        [](const double& w, const bsim::Dk2Nu* nu, const int& i_nuray, const TObject* extW)
                        { TSpline3* xsecspline = (TSpline3*)extW;
                          double energy = nu->nuray[i_nuray].E;
                          double xsec = xsecspline->Eval(energy);
                          if(xsec < 0.) {
                            xsec = 0.; // Just in case!
                          }
                          return w*xsec; });

  // Now we need to actual make the cross section spline
  XSec* xsec = new XSec();
  TSpline3* spline = xsec->GetXSec(14, "CH2", "tot_cc");

  // Since this is JUST muon neutrinos, let's remove other flavors
  p.RemoveNuFlav(-14);
  p.RemoveNuFlav(+12);
  p.RemoveNuFlav(-12);

  // Now let's make the Spectra
  // Note the last argument (which is optional):
  // we have to give the Spectra the external weights object!
  fr->AddSpectra(p, "xsec", "Energy (GeV)", Bins(100, 0., 10.), kEnergy, kAppXSec, spline);

  TFile* out = new TFile("/nova/ana/users/gkafka/FluxReader/demo2.root", "RECREATE");

  fr->ReadFlux(out);
  out->Close();
  delete fr;
  delete xsec;

  // If this was a lot of information, don't worry--
  // Most Vars and Weights should be predefined in Vars.h, Weights.h, and Weight.h
  // If you need to create a new one, try to follow this demo closely,
  // check other documentaiton, ask questions, and keep trying!
}

#endif
