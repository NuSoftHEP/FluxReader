#pragma once

// Package Includes
#include "Spectra.h"

// Forward Class Definitions
class TH1D;

namespace flxrd
{
  /// One dimensional implementation of the abstract Spectra class
  class Spectra1D: public Spectra
  {
  public:
    friend class FluxReader;

    ~Spectra1D();

    /// Access one of the histograms
    TH1* GetHist(int i_hist);

  protected:
    /// Fill one of the histograms with an entry
    /// The correct histogram will be determined from fParams,
    /// which was declared in the abstract base Spectra class
    void Fill(bsim::Dk2Nu* nu, std::map<std::string, int> nurayIndices);

    void WriteHists(TDirectory* out);

  private:
    /// Spectra1D specific constructor
    Spectra1D(Parameters params, std::string title,
              std::string labelx, std::vector<double> binsx, const Var& varx,
              const Weight& wei, TObject* extWeights = nullptr);

    /// Creates the histograms
    /// Called inside the constructor
    void CreateHists(std::string labelx, std::vector<double> binsx);

    std::vector<TH1D*> fHists; ///< Vector of 1D histograms
  };

}
