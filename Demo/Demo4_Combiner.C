// This fourth demo introduces the Combiner class
// It details how to use the class, and its limitations

#ifdef __CINT__
void Demo4_Combiner()
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
#include "Utilities.h"
#include "Vars.h"

using namespace flxrd;

void Demo4_Combiner()
{
  Parameters p(false);

  // Add a couple of detectors
  p.AddDetector(kNOvA_ND);
  p.AddDetector(kNOvA_FD);

  string dk2nu_loc = "/nusoft/data/flux/blackbird-numix/flugg_mn000z200i_rp11_lowth_pnut_f11f093bbird/dk2nu/";
  dk2nu_loc += "*dk2nu.root";
  FluxReader *fr = new FluxReader(dk2nu_loc, 2);

  // Add a Spectra object
  fr->AddSpectra(p, "enu", "Energy (GeV)", Bins(100, 0., 10.), kEnergy);

  // Add this "detector correlated" Spectra again
  fr->AddSpectra(p, "bmmat", "NOvA-ND", "NOvA-FD", "Energy(GeV)", Bins(100, 0., 10.), kEnergy);

  TFile* out = new TFile("/nova/ana/users/gkafka/FluxReader/demo4.root", "RECREATE");

  fr->ReadFlux(out);
  out->Close();
  delete fr;

  // So far, there has been nothing new
  // Now let's introduce the Combiner class, which can add together histograms in an automated way
  // When it is constructed, it is given a file to open
  Combiner* c = new Combiner("/nova/ana/users/gkafka/FluxReader/demo4.root");

  // This can combine all histograms with like neutrino flavor or like neutrino parent
  // It always requires the cross section and detector to be the same
  // To combine all flavors and parents, use the CombineAll function
  c->CombineAll();

  // The other two functions that exist are CombineNuFlavs() and CombineParents()
  // These are called in the same fashion, i.e., c->CombineNuFlavs();
  // All three functions write their contents out in the appropriate Spectra and detector directory
  // CombineAll() calls CombineNuFlavs() and CombineParents() as part of its work,
  // so simply calling this function should be good enough
  delete c;

  // In more detail, CombineNuFlavs() looks at all plots with like parent, cross section, and detector,
  // and it adds them all together
  // This gives all neutrinos that decay from a particular meson parent
  // Likewise CombineParents() gives spectra for all neutrinos of a specific flavor,
  // and CombineAll() gives all neutrinos (with the same cross section and detector)

  // Combiner only works on 1D, 2D, and 3D Spectra, but not on "detector correlated" Spectra
  // The detector correlated spectra have normalizations applied,
  // so they cannot be simply added using the TH1::Add function
  // Consequently, the combined plots are always made for this type of Spectra,
  // but a Combiner must be used for the other types of Spectra
  // This design decision allows the user to decide later to combine Spectra,
  // something that would have been a possibility for detector correlated Spectra--
  // but only if the aforementioned normalizations were saved to file as well
}

#endif
