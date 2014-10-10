// This second demo takes a much deeper look into the Parameters object
// It introduces adding and removing parameters,
// and the intertwining of this functionality with FluxReader

#ifdef __CINT__
void Demo1_Parameters()
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
#include "Detectors.h"
#include "FluxReader.h"
#include "Parameters.h"
#include "Utilities.h"
#include "Vars.h"

using namespace flxrd;

void Demo1_Parameters()
{
  // Recall the default parameters (ignoring parent sign):
  // The default flavors are nue, anti-nue, numu, and anti-numu
  // The default parents are muons, pions, kaons, and K-Long
  // The default cross sections are none, CC, and NC
  Parameters p(false);

  // Add a couple of detectors
  p.AddDetector(kNOvA_ND);
  p.AddDetector(kNOvA_FD);

  // What if we don't care about neutrinos from muons?
  // These can be easily removed in 3 different ways
  // By PDG
  p.RemoveParent(13);
  // By name
  p.RemoveParent("muon");
  // By Parent object
  p.RemoveParent(Parent::kMuon);

  // But let's add it back in, which requires a Parent object
  p.AddParent(Parent::kMuon);

  // We can remove neutrino flavors in the same way as Parents
  // But to "add" them, we have to reset them
  p.ResetNuFlavs();
  // This put nutaus back in the mix; let's get rid of them
  p.RemoveNuFlav(16);
  p.RemoveNuFlav("anutau");

  // Cross sections can be added and removed
  // The associated functions takes a string
  // The string should match something in XSec::ListIntTypes()
  p.AddXSec("tot_cc_p");
  p.RemoveXSec("tot_cc_p");
  // Note that order matters--if the above lines were reversed,
  // Nothing would have been removed since tot_cc_p is not a default
  // But then the cross section would have been added!

  string dk2nu_loc = "/nusoft/data/flux/blackbird-numix/flugg_mn000z200i_rp11_lowth_pnut_f11f093bbird/dk2nu/";
  dk2nu_loc += "*dk2nu.root";
  FluxReader *fr = new FluxReader(dk2nu_loc, 2);

  // Add a Spectra object
  fr->AddSpectra(p, "enu1", "Energy (GeV)", Bins(100, 0., 10.), kEnergy);

  // Order also matters to Spectra
  // The current configuration of the Parameters object is what will get created
  // If we removed a detector now:
  p.RemoveDetector("NOvA-FD");
  // Then added a new Spectra (which will be the same, otherwise, for comparison purposes)
  fr->AddSpectra(p, "enu2", "Energy (GeV)", Bins(100, 0., 10.), kEnergy);
  // This Spectra will have the same flavors, parents and cross sections,
  // but it will NOT include the NOvA FD

  TFile* out = new TFile("/nova/ana/users/gkafka/FluxReader/demo1.root", "RECREATE");

  fr->ReadFlux(out);
  out->Close();
  delete fr;

  // This has not shown all of the Add/Remove functions for Parameters;
  // check other documentation for this
  // However, the functions do all have similar functionality,
  // so none are wildly different from the others
}

#endif
