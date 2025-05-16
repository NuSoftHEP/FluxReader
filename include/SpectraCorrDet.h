#pragma once

// Package Includes
#include "Spectra.h"

// Forward Class Definitions
class TH1D;
class TH2D;

namespace flxrd
{
  /// Implementation of the abstract Spectra class correlating two detectors
  /// See the documentation for the abstract or 1D implementation for more details
  class SpectraCorrDet: public Spectra
  {
  public:
    friend class FluxReader;

    ~SpectraCorrDet();

    TH1* GetHist(int i_hist);

  protected:
    /// Fill the full 2D histograms and associated normalization histograms
    /// fHists gets filled using the weight from the detY neutrino ray
    /// fNorms gets filled using the weight from the detX neutrino ray
    void Fill(bsim::Dk2Nu* nu, std::map<std::string, int> nurayIndices);

    void WriteHists(TDirectory* out);

  private:
    SpectraCorrDet(Parameters params, std::string title,
                   std::string detX, std::string detY,
                   std::string labelx, std::vector<double> binsx, const Var& varx,
                   const Weight& wei, TObject* extWeights = nullptr);

    /// These functions combine histograms, much like a flxrd::Combiner
    /// However, the normalization histograms do not get written to file,
    /// so combining these plots is not (correctly) possible with a Combiner
    /// These functions combine the histograms correctly, before normalization
    void CombineNuFlavs(std::vector<TH2D*>& newHists, std::vector<TH1D*>& newNorms);
    void CombineParents(std::vector<TH2D*>& newHists, std::vector<TH1D*>& newNorms);
    void CombineAll();
/// FIX THE NAMING ISSUE
    void CreateHists(std::string detX, std::string detY,
                     std::string labelx, std::vector<double> binsx);

    /// For each histogram in fHists,
    /// normalize each column in fHists (variable value at detX)
    /// by the corresponding fNorms histogram bin (total event weight at detX)
    void Normalize();

    std::vector<TH2D*> fHists; ///< Vector of 2D histograms of detX vs detY
    std::vector<TH1D*> fNorms; ///< Vector of 1D histograms of events at detX

    int i_detX; ///< Index of the x axis detector in the internal Parameters object
    int i_detY; ///< Index of the y axis detector in the internal Parameters object

    bool fIsNormalized; ///< Helper to determine whether fHists have been normalized by fNorms yet
    bool fAlreadyCombined; ///< Helper to determine whether fHists have been combined yet
  };

}
