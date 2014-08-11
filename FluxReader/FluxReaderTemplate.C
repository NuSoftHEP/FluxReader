// This template script includes all of, and nothing more than,
// the basic necessities for running FluxReader

#ifdef __CINT__
void FluxReaderTemplate()
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

void FluxReaderTemplate()
{
  // First, we need to set up a Parameters object
  Parameters p(false);

  // Add at least one detector
  p.AddDetector(knova_fd);

  // Next, we'll create a FluxReader object
  string dk2nu_loc = "/nusoft/data/flux/dk2nu/nova/2010/flugg_mn000z200i_20101117.gpcfgrid_lowth/";
  dk2nu_loc += "*dk2nu.root";
  FluxReader *fr = new FluxReader(dk2nu_loc, 2);

  // The FluxReader needs to generate something!
  fr->AddSpectra(p, "enu", "Energy (GeV)", Bins(100, 0., 10.), kEnergy);

  // The last thing to do before running is to set up an output file
  TFile* out = new TFile("/nova/ana/users/gkafka/FluxReader/HelloWorld.root", "RECREATE");

  // This function loops over the files, and fills all the histograms!
  fr->ReadFlux(out);
  out->Close(); // Close output file
  delete fr; // Clean up
}

#endif
