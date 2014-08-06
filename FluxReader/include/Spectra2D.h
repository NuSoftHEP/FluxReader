#pragma once

// Package Includes
#include "Spectra.h"

// Forward Class Definitions
class TH2D;

namespace flxrd
{
  /// Two dimensional implementation of the abstract Spectra class
  /// See the documentation for the abstract or 1D implementation for more details
  class Spectra2D: public Spectra
  {
  public:
    friend class FluxReader;

    ~Spectra2D();

    TH1* GetHist(int i_hist);

  protected:
    void Fill(bsim::Dk2Nu* nu, std::map<std::string, int> nurayIndices);

    void WriteHists(TDirectory* out);

    Var fVarY;

  private:
    /// Spectra2D specific constructor
    Spectra2D(Parameters params, std::string title,
              std::string labelx, std::vector<double> binsx, const Var& varx,
              std::string labely, std::vector<double> binsy, const Var& vary,
              const Weight& wei, TObject* extWeights = nullptr);

    void CreateHists(std::string labelx, std::vector<double> binsx,
                     std::string labely, std::vector<double> binsy);

    std::vector<TH2D*> fHists; ///< Vector of 2D histograms
  };
}
