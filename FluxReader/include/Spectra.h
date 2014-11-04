#pragma once

// C/C++ Includes
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

// Package Includes
#include "Detector.h"
#include "Parameters.h"
#include "Var.h"
#include "Weight.h"

// Forward Class Definitions
class TDirectory;
class TH1;
class TObject;
class TSpline3;

namespace flxrd
{
  /// This abstract class sets up some common elements for a dimensional Spectra
  /// The version with a set dimension will include a vector of histograms
  class Spectra
  {
  public:
    friend class FluxReader;

    ~Spectra();
 
    /// Access one of the histograms
    virtual TH1* GetHist(int i_hist) = 0;

    /// Access the common piece of the titles
    std::string GetTitle() const { return fTitle; }

  protected:
    /// These inputs will be common to any dimensional Spectra
    Spectra(Parameters params, std::string title,
            const Var& varx, const Weight& wei,
            TObject* extWeights = nullptr);

    /// These are the branches necessary for the Var and Weight
    std::set<std::string> BranchesToAdd() const { return fBranches; }

    /// Returns all detectors needed for the Spectra
    std::set<Detector> Detectors() const;

    /// Fill one of the histograms with an entry
    /// The correct histogram will be determined from fParams
    virtual void Fill(bsim::Dk2Nu* nu, std::map<std::string, int> nurayIndices) = 0;

    /// Fills the cross section spline map
    void SetupXSec();

    /// Write all of the histograms in the input directory
    virtual void WriteHists(TDirectory* dir) = 0;

    /// Create a cross section label to identifty specific splines
    std::string XSecName();

    std::set<std::string> fBranches; ///< List of flux file branches needed to be activated

    const double fDefaultWeightCorrection = 1./(10000. * M_PI);

    TObject* fExtWeights; ///< External weights to be applied to histogram entries

    Parameters fParams; ///< The Parameters to be used

    std::string fTitle; ///< Label to prefix all of the histograms

    Var    fVarX; ///< Variable to fill the x axis
    Weight fWei; ///< How to weight each entry

    std::map<std::string, TSpline3*> fXSecSplines; ///< Map of cross section splines
  };
}
