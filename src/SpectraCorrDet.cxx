#include "SpectraCorrDet.h"

// C/C++ Includes
#include <cassert>
#include <iostream>

// Root Includes
#include "TDirectory.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TObject.h"
#include "TSpline.h"

namespace flxrd
{
  //---------------------------------------------------------------------------
  SpectraCorrDet::SpectraCorrDet(Parameters params, std::string title,
                                 std::string detX, std::string detY,
                                 std::string labelx, std::vector<double> binsx, const Var& varx,
                                 const Weight& wei, TObject* extWeights)
    : Spectra(params, title, varx, wei, extWeights), fIsNormalized(false), fAlreadyCombined(false)
  {
    assert(params.NDet() >= 2); // There need to be at least two detectors in the Parameters object

    // Initialize the two indices
    i_detX = -1;
    i_detY = -1;

    // Loop through the input parameters object,
    // look for detectors that match the names of detX and detY inputs
    for(unsigned int i_det = 0, n_det = params.NDet(); i_det < n_det; ++i_det) {
      if(!detX.compare(params.GetDetName(i_det))) {
        i_detX = i_det; // Store the index that matched detX
      }
      if(!detY.compare(params.GetDetName(i_det))) {
        i_detY = i_det; // Store the index that matched detY
      }
    }

    assert(i_detX != -1 && i_detY != -1); // Make sure both detectors were found

    fParams.SetCurrentDet(i_detY); // Set the current detector to detY

    CreateHists(detX, detY, labelx, binsx);
  }

  //---------------------------------------------------------------------------
  TH1* SpectraCorrDet::GetHist(int i_hist)
  {
    // Normalize the histograms by fNorms if they have not already been
    if(!fIsNormalized) {
      Normalize();
    }

    const int n_hist = fHists.size();
    if(i_hist < 0 || i_hist >= n_hist) {
      std::cout << "Input histogram index is out of range." << std::endl;
      assert(false);
    }

    return fHists[i_hist];
  }

  //---------------------------------------------------------------------------
  void SpectraCorrDet::Fill(bsim::Dk2Nu* nu, std::map<std::string, int> nurayIndices)
  {
    int nuPDG = nu->decay.ntype;

    if(!fParams.SetCurrentNuFlav(nuPDG)) {
      return;
    }

    int parPDG = ( (fParams.IsSignSensitive()) ?     GetAncestorPDG(nu)
                                               : abs(GetAncestorPDG(nu)) );

    if(!fParams.SetCurrentParent(parPDG)) {
      return;
    }

    // Get the first and last indices in the NuRay vector corresponding to the x axis detector
    int first_nuray_x = nurayIndices[fParams.GetDetName(i_detX)];
    int last_nuray_x  = first_nuray_x + fParams.GetDetector(i_detX).GetUses();
    if(last_nuray_x == first_nuray_x) {
      ++last_nuray_x;
    }

    // Get the first and last indices in the NuRay vector corresponding to the y axis detector
    int first_nuray_y = nurayIndices[fParams.GetDetName(i_detY)];
    int last_nuray_y  = first_nuray_y + fParams.GetDetector(i_detY).GetUses();
      if(last_nuray_y == first_nuray_y) {
        ++last_nuray_y;
      }

    fParams.SetCurrentDet(i_detY); // Set the current detector to detY

    for(int i_xsec = 0, n_xsec = fParams.NXSec(); i_xsec < n_xsec; ++i_xsec) {
      fParams.SetCurrentXSec(i_xsec);

      int i_hist = fParams.GetCurrentMaster() - fParams.MaxMaster(i_detY - 1); // Get the correct histogram index

      for(int i_nuray_x = first_nuray_x; i_nuray_x < last_nuray_x; ++i_nuray_x) {
        // Calculate the standard weight at the x axis detector
        double weight_x =   nu->decay.nimpwt * nu->nuray[i_nuray_x].wgt
                          * fXSecSplines[XSecName()]->Eval(nu->nuray[i_nuray_x].E)
                          * fDefaultWeightCorrection;

        for(int i_nuray_y = first_nuray_y; i_nuray_y < last_nuray_y; ++i_nuray_y) {
          // Calculate the standard weight at the y axis detector
          double weight_y =   nu->decay.nimpwt * nu->nuray[i_nuray_y].wgt
                            * fXSecSplines[XSecName()]->Eval(nu->nuray[i_nuray_y].E)
                            * fDefaultWeightCorrection;

          // Evaluate the variables and weights, fill the histograms
          // Both axes variables evaluate fVarX, but the x axis is evaluated at detX, and the y axis at detY
          // The weight applied is the weight at detY
          fHists[i_hist]->Fill(fVarX(nu, i_nuray_x), fVarX(nu, i_nuray_y), fWei(weight_y, nu, i_nuray_y, fExtWeights));

          // This evaluates the weight at detX
          fNorms[i_hist]->Fill(fVarX(nu, i_nuray_x),                       fWei(weight_x, nu, i_nuray_x, fExtWeights));
        } // Loop over y detector uses
      } // Loop over x detector uses
    } // Loop over cross sections

    return;
  }

  //---------------------------------------------------------------------------
  void SpectraCorrDet::WriteHists(TDirectory* out)
  {
    // Normalize the histograms by fNorms if they have not already been
    if(!fIsNormalized) {
      Normalize();
    }

    TDirectory* temp = gDirectory;

    out->cd();

    // These histograms do not need detector directories, so just write them
    for(unsigned int i_hist = 0, n_hist = fHists.size(); i_hist < n_hist; ++i_hist) {
      gDirectory->WriteTObject(fHists[i_hist]);
    }

    temp->cd();
    return;
  }

  //---------------------------------------------------------------------------
  void SpectraCorrDet::CombineNuFlavs(std::vector<TH2D*>& newHists,
                                      std::vector<TH1D*>& newNorms)
  {
    std::string rep_str = "allnu"; // This string will replace the neutrino flavor name

    // Store number of each parameter in this Spectra
    const unsigned int n_flav = fParams.NFlav();
    const unsigned int n_par  = fParams.NPar();
    const unsigned int n_xsec = fParams.NXSec();

    for(unsigned int i_xsec = 0; i_xsec < n_xsec; ++i_xsec) {
      for(unsigned int i_par  = 0; i_par  < n_par;  ++i_par) {
        // This corresponds to the way Parameters does indexing.
        // NuFlav index 0 is used by not including "+ i_flav" at the end.
        int index =   n_flav*n_par*n_xsec*i_detY
                    + n_flav*n_par*i_xsec
                    + n_flav*i_par;

        // CreateHists does not loop over detectors,
        // So offset the fHists index to compensate
        int i_hist = index - fParams.MaxMaster(i_detY - 1);

        // Find/create  a stored histogram and copy it into a new histogram for combining
        TH2D* hHist = (TH2D*)fHists[i_hist]->Clone();
        TH1D* hNorm = (TH1D*)fNorms[i_hist]->Clone();

        for(unsigned int i_flav = 1; i_flav < n_flav; ++i_flav) {
          ++i_hist; // This corresponds to an increment of the NuFlav index

          // Add in the next histograms
          hHist->Add(fHists[i_hist]);
          hNorm->Add(fNorms[i_hist]);
        }

        // Replace neutrino flavor name by the replacement string
        // The histogram name has format title_nuflav_par_xsec_det
        std::string hName = hHist->GetName(); // Get the name
        int firstPos = hName.find('_');
        int secndPos = hName.find('_', firstPos+1); // Search only after position specified in second argument
        hName.replace(firstPos+1, secndPos-firstPos-1, rep_str);
        hHist->SetName(hName.c_str()); // Give the combined histogram the correct name for writing to file

        // Store the final added copies
        newHists.push_back(hHist);
        newNorms.push_back(hNorm);
      } // Loop over parents
    } // Loop over cross sections

    return;
  }

  //---------------------------------------------------------------------------
  void SpectraCorrDet::CombineParents(std::vector<TH2D*>& newHists,
                                      std::vector<TH1D*>& newNorms)
  {
    // This function is similar to CombineNuFlav; see its comments for more details.

    std::string rep_str = "allpar"; // This string will replace the parent name

    const unsigned int n_flav = fParams.NFlav();
    const unsigned int n_par  = fParams.NPar();
    const unsigned int n_xsec = fParams.NXSec();

    for(unsigned int i_xsec = 0; i_xsec < n_xsec; ++i_xsec) {
      for(unsigned int i_flav = 0; i_flav < n_flav; ++i_flav) {
        // This corresponds to the way Parameters does indexing.
        // Parent index 0 is used by not including the line "n_flav*i_par".
        int index =   n_flav*n_par*n_xsec*i_detY
                    + n_flav*n_par*i_xsec
                    + i_flav;

        int i_hist = index - fParams.MaxMaster(i_detY - 1);

        TH2D* hHist = (TH2D*)fHists[i_hist]->Clone();
        TH1D* hNorm = (TH1D*)fNorms[i_hist]->Clone();

        for(unsigned int i_par  = 1; i_par  < n_par;  ++i_par) {
          i_hist += n_flav; // This corresponds to an increment of the Parent index

          hHist->Add(fHists[i_hist]);
          hNorm->Add(fNorms[i_hist]);
        }

        std::string hName = hHist->GetName();
        int firstPos = hName.find('_');
        firstPos = hName.find('_', firstPos+1); // Update to the next occurence of '_'
        int secndPos = hName.find('_', firstPos+1);
        hName.replace(firstPos+1, secndPos-firstPos-1, rep_str);
        hHist->SetName(hName.c_str());

        newHists.push_back(hHist);
        newNorms.push_back(hNorm);
      } // Loop over flavors
    } // Loop over cross sections

    return;
  }

  //---------------------------------------------------------------------------
  void SpectraCorrDet::CombineAll()
  {
    // This function is similar to CombineNuFlav; see its comments for more details.

    std::string rep_str  = "allnu"; // Replacement string for neutrino flavor name

    std::vector<TH2D*> vecCombinedNuFlavHists;
    std::vector<TH1D*> vecCombinedNuFlavNorms;
    std::vector<TH2D*> vecCombinedParentHists;
    std::vector<TH1D*> vecCombinedParentNorms;

    CombineNuFlavs(vecCombinedNuFlavHists, vecCombinedNuFlavNorms); // Combine neutrino flavors
    CombineParents(vecCombinedParentHists, vecCombinedParentNorms); // Combine neutrino parents

    // Store the combined flavor and parent histograms in main histogram vectors
    for(unsigned int i = 0, n = vecCombinedNuFlavHists.size(); i < n; ++i) {
      fHists.push_back(vecCombinedNuFlavHists[i]);
      fNorms.push_back(vecCombinedNuFlavNorms[i]);
    }
    for(unsigned int i = 0, n = vecCombinedParentHists.size(); i < n; ++i) {
      fHists.push_back(vecCombinedParentHists[i]);
      fNorms.push_back(vecCombinedParentNorms[i]);
    }

    const unsigned int n_flav = fParams.NFlav();
    const unsigned int n_xsec = fParams.NXSec();

    // This block loops through the combined parent vector, and combines flavors
    // The combined parent vector should have n_xsec*n_flav entries
    // All histograms of like cross section will be consecutive; see CombineParents
    for(unsigned int i_xsec = 0; i_xsec < n_xsec; ++i_xsec) {
      int index = i_xsec*n_flav; // The first index of cross section i_xsec

      // Get the combined parent histograms
      TH2D* hHist = (TH2D*)vecCombinedParentHists[index]->Clone();
      TH1D* hNorm = (TH1D*)vecCombinedParentNorms[index]->Clone();

      for(unsigned int i_flav = 1; i_flav < n_flav; ++i_flav) {
        ++index; // Increment the flavor index

        hHist->Add(vecCombinedParentHists[index]);
        hNorm->Add(vecCombinedParentNorms[index]);
      }

      // The parent name is already replaced by pulling the combined parent histograms
      // Now replace the neutrino flavor name
      std::string hName = hHist->GetName();
      int firstPos = hName.find('_');
      int secndPos = hName.find('_', firstPos+1);
      hName.replace(firstPos+1, secndPos-firstPos-1, rep_str);
      hHist->SetName(hName.c_str());

      fHists.push_back(hHist);
      fNorms.push_back(hNorm);
    } // Loop over cross sections

    return;
  }

  //---------------------------------------------------------------------------
  void SpectraCorrDet::CreateHists(std::string detX, std::string detY, 
                                   std::string labelx, std::vector<double> binsx)
  {
    std::string hist_title = ""; // This will become the title for writing to file
    std::string both_det_str = detX + "_" + detY;
    std::string axis_label = ";" + detX + " " + labelx + ";" + detY + " " + labelx; // This is the x and y axes labels
    const int nBinsX = binsx.size() - 1;

    for(int i = fParams.MaxMaster(i_detY - 1), n = fParams.MaxMaster(i_detY); i < n; ++i) {
      hist_title = fTitle + "_" + fParams.NameTag(i); // Create the full title for writing to file
      hist_title.replace(hist_title.length() - detY.length(), detY.length(), both_det_str);

      // Create the 2D histogram of detX vs detY
      TH2D* h2 = new TH2D(hist_title.c_str(), axis_label.c_str(),
                          nBinsX, &binsx[0], nBinsX, &binsx[0]);
      fHists.push_back(h2);

      // Create the 1D histogram of detX events
      TH1D* h1 = new TH1D("", "", nBinsX, &binsx[0]);
      fNorms.push_back(h1);
    }

    return;
  }

  //---------------------------------------------------------------------------
  void SpectraCorrDet::Normalize()
  {
    // The histograms must be combined BEFORE normalizing,
    // but this should only be done once
    if(!fAlreadyCombined) {
      CombineAll();
    }
    fAlreadyCombined = true;

    // Store the number of x axis and y axis bins to loop over
    const int X = fHists[0]->GetNbinsX();
    const int Y = fHists[0]->GetNbinsY();

    for(int i_hist = 0, n_hist = fHists.size(); i_hist < n_hist; ++i_hist) {
      for(int i = 0; i <= X+1; ++i) {
        for(int j = 0; j <= Y+1; ++j) {
          // Get the normalization value
          double detX_norm_value = fNorms[i_hist]->GetBinContent(i);

          // Divide by the normalization value if it is greater than 0
          if(detX_norm_value > 0) {
            fHists[i_hist]->SetBinContent(i,j, fHists[i_hist]->GetBinContent(i,j)
                                              /detX_norm_value);
          }
          else { // Otherwise set the bin to 0
            fHists[i_hist]->SetBinContent(i,j, 0);
          }
        } // Loop over y axis
      } // Loop over x axis
    } // Loop over histograms

    fIsNormalized = true;

    return;
  }
}
