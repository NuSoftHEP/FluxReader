// This first demo script shows the basics for running FluxReader
// This can be considered the bare essentials necessary to run FluxReader
// It introduces the Parameters class and FluxReader class,
// and the essential functions necessary for proper running conditions

#ifdef __CINT__
void Demo0_FluxReader()
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

void Demo0_FluxReader()
{
  // First, we need to set up the parameters to run over
  // Just by constructing the object, defaults are set for
  // neutrino flavors, parents, and cross sections.
  // The default flavors are nue, anti-nue, numu, and anti-numu
  // The default parents are muons, pions, kaons, and K-Long
  // The boolean input determines whether the parent sign is considered (true), or ignored (false)
  // The default cross sections are none, CC, and NC
  Parameters p(false);

  // Add a couple of detectors
  // These are predefined in Detectors.h
  p.AddDetector(kNOvA_ND);
  p.AddDetector(kNOvA_FD);

  // Next, we'll create a FluxReader object
  // In the most basic constructor, we only provide a wildcard path name to the input files
  // The constructor will expand this into all files matching this string
  // Optional second and third arguments can be added--
  // the second sets the number of files to use;
  // the third sets the number of files to skip over
  string dk2nu_loc = "/nusoft/data/flux/dk2nu/nova/2010/flugg_mn000z200i_20101117.gpcfgrid_lowth/";
  dk2nu_loc += "*dk2nu.root";
  FluxReader *fr = new FluxReader(dk2nu_loc, 2);

  // The FluxReader needs to generate something!
  // This is the most basic construction of a Spectra object,
  // and will generate a set of 1D histograms
  // The first input is the Parameters object, and this determines the "set" of histograms
  // The second input is the Spectra label; it will be the directory label in the output file,
  // and it will also be part of the title of each histogram
  // The third input is the x axis label
  // The fourth input is a vector of bin edges
  // Bins is a function defined in Utilities.h and returns a vector of bin edges,
  // so this standard method creates a set of equally sized bins
  // The fifth input is a Var object, and determines what is filled in the histograms
  // kEnergy is one of the predefined Vars in Vars.h
  fr->AddSpectra(p, "enu", "Energy (GeV)", Bins(100, 0., 10.), kEnergy);

  // The last thing to do before running is to set up an output file
  TFile* out = new TFile("/nova/ana/users/gkafka/FluxReader/demo0.root", "RECREATE");

  // This function loops over the files, and fills all the histograms!
  fr->ReadFlux(out);
  out->Close(); // Close output file
  delete fr; // Clean up
}

#endif
