#pragma once

// Package Includes
#include "Spectra.h"

// Forward Class Definitions
class TH3D;

namespace flxrd
{
  /// Three dimensional implementation of the abstract Spectra class
  /// See the documentation for the abstract or 1D implementation for more details
  class Spectra3D: public Spectra
  {
  public:
    friend class FluxReader;

    ~Spectra3D();

    TH1* GetHist(int i_hist);

  protected:
    void Fill(bsim::Dk2Nu* nu, std::map<std::string, int> nurayIndices);

    void WriteHists(TDirectory* out);

    Var fVarY;
    Var fVarZ;

  private:
    /// Spectra3D specific constructor
    Spectra3D(Parameters params, std::string title,
              std::string labelx, std::vector<double> binsx, const Var& varx,
              std::string labely, std::vector<double> binsy, const Var& vary,
              std::string labelz, std::vector<double> binsz, const Var& varz,
              const Weight& wei, TObject* extWeights = nullptr);

    void CreateHists(std::string labelx, std::vector<double> binsx,
                     std::string labely, std::vector<double> binsy,
                     std::string labelz, std::vector<double> binsz);

    std::vector<TH3D*> fHists; ///< Vector of 3D histograms
  };
}
