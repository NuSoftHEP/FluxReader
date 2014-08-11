// This fourth demo explores the AddSpectra functions in more detail,
// showing an example of creating a 2D Spectra
// It also shows some flexibility with binning, and variably sized bins
// Lastly, it demonstrates "detector correlated" Spectra,
// great for generating beam matrices

#ifdef __CINT__
void Demo3_Spectra()
{
  std::cout << "Sorry, you must run in compiled mode." << std::endl;
}
#else

// C/C++ Includes
#include <iostream>
#include <string>
#include <vector>

// ROOT Includes
#include "TFile.h"

// Package Includes
#include "Detectors.h"
#include "FluxReader.h"
#include "Parameters.h"
#include "Utilities.h"
#include "Vars.h"

using namespace flxrd;

void Demo3_Spectra()
{
  Parameters p(false);

  // Add a couple of detectors
  p.AddDetector(knova_nd);
  p.AddDetector(knova_fd);

  string dk2nu_loc = "/nusoft/data/flux/dk2nu/nova/2010/flugg_mn000z200i_20101117.gpcfgrid_lowth/";
  dk2nu_loc += "*dk2nu.root";
  FluxReader *fr = new FluxReader(dk2nu_loc, 2);

  // Add a Spectra object
  fr->AddSpectra(p, "enu1", "Energy (GeV)", Bins(100, 0., 10.), kEnergy);

  // So far we've only dealt with 1D Spectra, but there are others
  // FluxReader can make 2D, 3D, and "detector correlated" Spectra
  // The 2D and 3D Spectra simply take extra arguments for their y and z axes
  // (The optional arguments shown in the last demo ALWAYS come at the end)
  // Here is a pT pz histogram as an example of a 2D Spectra
  // Note that kpz and kpT are both defined in Vars.h
  fr->AddSpectra(p, "pTpz", "p_{z} (GeV)", Bins(120, 0., 120.), kpz,
                            "p_{T} (GeV)", Bins(40,  0., 4.),   kpT);

  // These examples have all been using the Bins function, but this is not required
  // The particular input only requires a std::vector<double>
  // Consequently, FluxReader can handle varied bins
  // The following code snippet by Chris Backhouse,
  // from novasoft's CAFAna/Binning.cxx, and generates a vector of bin edges
  // ***** //
  const int kNumTrueEnergyBins = 100;

  // N+1 bin low edges
  std::vector<double> edges(kNumTrueEnergyBins+1);

  const double Emin = .5; // 500 MeV: there's really no events below there

  // How many edges to generate. Allow room for 0-Emin bin
  const double N = kNumTrueEnergyBins-1;
  const double A = N*Emin;

  edges[0] = 0;

  for(int i = 1; i <= N; ++i){
    edges[kNumTrueEnergyBins-i] = A/i;
  }

  edges[kNumTrueEnergyBins] = 120; // Replace the infinity that would be here
  // ***** //

  // Check out the difference yourself!
  fr->AddSpectra(p, "enu2", "Energy (GeV)", edges, kEnergy);

  // Back to Spectra for the "detector correlated" example
  // This Spectra plots the same variable on the x and y axes, but at different detectors
  // After the Spectra title is given (2nd input),
  // the 3rd and 4th inputs are strings for the detector names
  // Both of the detectors MUST be in the Parameters object given, and all others are ignored
  // After these inputs, the remaining inputs be have the same way as a 1D Spectra
  // This particular example is also called a beam matrix
  // When examining the output, you may notice some extra plots with this type of Spectra--
  // those will be the subject of the next tutorial
  fr->AddSpectra(p, "bmmat", "nova_nd", "nova_fd", "Energy(GeV)", Bins(100, 0., 10.), kEnergy);

  TFile* out = new TFile("/nova/ana/users/gkafka/FluxReader/demo3.root", "RECREATE");

  fr->ReadFlux(out);
  out->Close();
  delete fr;
}

#endif
