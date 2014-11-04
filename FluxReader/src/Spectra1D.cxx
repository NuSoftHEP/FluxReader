#include "Spectra1D.h"

// C/C++ Includes
#include <cassert>
#include <iostream>

// Root Includes
#include "TDirectory.h"
#include "TH1.h"
#include "TH1D.h"
#include "TSpline.h"

namespace flxrd
{
  //---------------------------------------------------------------------------
  Spectra1D::Spectra1D(Parameters params, std::string title,
                       std::string labelx, std::vector<double> binsx, const Var& varx,
                       const Weight& wei, TObject* extWeights)
    : Spectra(params, title, varx, wei, extWeights)
  {
    CreateHists(labelx, binsx); // Set up the histograms
  }

  //---------------------------------------------------------------------------
  TH1* Spectra1D::GetHist(int i_hist)
  {
    // Check that there actually is a histogram to return
    const int n_hist = fHists.size();
    if(i_hist < 0 || i_hist >= n_hist) {
      std::cout << "Input histogram index is out of range." << std::endl;
      assert(false);
    }

    return fHists[i_hist];
  }

  //---------------------------------------------------------------------------
  void Spectra1D::Fill(bsim::Dk2Nu* nu, std::map<std::string, int> nurayIndices)
  {
    int nuPDG = nu->decay.ntype; // Get the neutrino flavor from the flux object

    if(!fParams.SetCurrentNuFlav(nuPDG)) {
      return;
    }

    // Get the neutrino parent PDG from the flux object (absolute value if applicable)
    int parPDG = ( (fParams.IsSignSensitive()) ?     nu->decay.ptype
                                               : abs(nu->decay.ptype) );

    if(!fParams.SetCurrentParent(parPDG)) {
      return;
    }

    for(int i_det = 0, n_det = fParams.NDet(); i_det < n_det; ++i_det) {
      fParams.SetCurrentDet(i_det);

      // Get the first and last indices in the NuRay vector corresponding to the current detector
      int first_nuray = nurayIndices[fParams.GetDetName(i_det)];
      int last_nuray  = first_nuray + fParams.GetDetector(i_det).GetUses();

      for(int i_xsec = 0, n_xsec = fParams.NXSec(); i_xsec < n_xsec; ++i_xsec) {
        fParams.SetCurrentXSec(i_xsec);

        int i_hist = fParams.GetCurrentMaster(); // Get the correct histogram index

        for(int i_nuray = first_nuray; i_nuray < last_nuray; ++i_nuray) {
          // Calculate the standard weight
          double weight =   nu->decay.nimpwt * nu->nuray[i_nuray].wgt
                          * fXSecSplines[XSecName()]->Eval(nu->nuray[i_nuray].E)
                          * fDefaultWeightCorrection;

          // Evaluate the variable and weight, fill the histogram
          fHists[i_hist]->Fill(fVarX(nu, i_nuray), fWei(weight, nu, i_nuray, fExtWeights));
        } // Loop over uses
      } // Loop over cross sections
    } // Loop over detectors

    return;
  }

  //---------------------------------------------------------------------------
  void Spectra1D::WriteHists(TDirectory* out)
  {
    TDirectory* temp = gDirectory; // Store the current directory to come back to later

    std::string det_name = ""; // Used to compare the current detector to the previous one

    for(const auto& index : fParams) { // Loop over all Parameters indices
      fParams.SetIndices(index); // Set the current master

      // Check if the detector has changed since the last Parameters master index
      if(det_name.compare(fParams.GetDetName(fParams.GetCurrentDet()))) {
        // Set the new detector name as the new "previous" string
        det_name = fParams.GetDetName(fParams.GetCurrentDet());
        out->cd(); // Go to the top level of the input directory

        // Check if this detector already has a directory made for it
        // Make the directory if not
        if( !(out->GetListOfKeys()->FindObject(det_name.c_str())) ) {
          out->mkdir(fParams.GetDetName(fParams.GetCurrentDet()).c_str());
        }

        // Go into the detector directory
        out->cd(fParams.GetDetName(fParams.GetCurrentDet()).c_str());
      }

      // Write the current histogram
      gDirectory->WriteTObject(fHists[index]);
    }

    temp->cd(); // Go back to the original directory
    return;
  }

  //---------------------------------------------------------------------------
  void Spectra1D::CreateHists(std::string labelx, std::vector<double> binsx)
  {
    std::string hist_title = ""; // This will become the title for writing to file
    std::string axis_label = ";"+labelx+";"; // This is the x axis label
    const int nBinsX = binsx.size() - 1;

    // Create a histogram for each Parameters master index
    for(int i = 0, n = fParams.MaxMaster(); i < n; ++i) {
      hist_title = fTitle + "_" + fParams.NameTag(i);
      TH1D* h1 = new TH1D(hist_title.c_str(), axis_label.c_str(), nBinsX, &binsx[0]);
      fHists.push_back(h1);
    }

    return;
  }
}
