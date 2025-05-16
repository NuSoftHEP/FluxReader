#pragma once

// C/C++ Includes
#include <map>
#include <set>
#include <string>

// Package Includes
#include "Detector.h"

// Forward Class Definitions
class TFile;
class TH1;
class TH2D;

namespace flxrd
{
  class Parameters;

  /// Combiner is a class which reads a FluxReader output file,
  /// reads its parameters, and can combine its contents
  class Combiner
  {
  public:
    /// \param out A string that is the path to a FluxReader output file
    Combiner(std::string out);

    ~Combiner();

    /// Combine histograms with common detector, cross section and decay parent
    /// The result is stored in the same detector and histogram folder,
    /// and the original histograms are retained unchanged
    void CombineNuFlavs();

    /// Combine histograms with common detector, cross section and neutrino flavor
    void CombineParents();

    /// Calls CombineNuFlavs() and CombineParents(),
    /// then takes the combined parent histograms and combines the neutrino flavors,
    /// i.e., the result combines all neutrinos 
    void CombineAll();

  private:
    /// Returns true if any histgorams in a FluxReader output file contain
    /// the string "search" in their name
    bool CombineAlreadyCalled(std::string search);

    /// Outputs the parameters and histogram types found in the FluxReader file
    /// that was input to the Combiner constructor 
    void InitialMessage();

    /// Sets up fParams to match the parameters found in the FluxReader file
    /// that was input to the Combiner constructor
    void SetupParameters(Parameters* params,
                         std::set<std::string>& dets,
                         std::set<std::string>& nuflavs,
                         std::set<std::string>& parents,
                         std::set<std::string>& xsecs);

    /// Map pointing from a Spectra name to its associated Parameters
    std::map<std::string, Parameters> fParamsMap;

    /// List of Spectra in the file given to the constructor
    std::set<std::string> fSpectra;

    TFile* fOut;
  };
}
