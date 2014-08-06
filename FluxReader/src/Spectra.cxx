#include "Spectra.h"

// Root Includes
#include "TF1.h"
#include "TObject.h"
#include "TSpline.h"

// Package Includes
#include "ParticleParam.h"
#include "Utilities.h"
#include "XSec.h"

// Other External Includes
#include "dk2nu.h"

namespace flxrd
{
  //---------------------------------------------------------------------------
  Spectra::Spectra(Parameters params, std::string title,
                   const Var& varx, const Weight& wei, TObject* extWeights)
    : fParams(params), fTitle(title), fVarX(varx), fWei(wei)
  {
    if(extWeights) {
      fExtWeights = extWeights;
    }
    else {
      fExtWeights = nullptr;
    }

    // These default branches are needed to split between neutrino flavor and parent type,
    // in order to fill or write the correct histogram, and to evaluate cross sections
    fBranches = {"nuray", "nuray.E", "nuray.wgt", "decay", "decay.ntype", "decay.ptype", "decay.nimpwt"};

    // Add variables needed by x axis variable and weight to the list of branches
    fBranches.insert(fVarX.Branches().begin(), fVarX.Branches().end());
    fBranches.insert(fWei .Branches().begin(), fWei .Branches().end());

    SetupXSec(); // Create and store the necessary cross section splines
  }

  //---------------------------------------------------------------------------
  std::set<Detector> Spectra::Detectors() const
  {
    std::set<Detector> ret;

    for(int i_det = 0, n_det = fParams.NDet(); i_det < n_det; ++i_det) {
      ret.insert(fParams.GetDetector(i_det));
    }

    return ret;
  }

  //---------------------------------------------------------------------------
  void Spectra::SetupXSec()
  {
    XSec* xsec = new XSec(); // Create XSec object

    std::string xsecname = "";

    for(const auto& index : fParams) { // Loop over all Parameters indices
      fParams.SetIndices(index); // Set the master index

      xsecname = XSecName(); // Generate a cross section label

      // Check if the label is a new one
      if(fXSecSplines.find(xsecname) == fXSecSplines.end()) {
        // Get the parameters to create a cross section spline: PDG, target nucleus, interaction current
        int pdg = fParams.GetNuFlavPDG(fParams.GetCurrentNuFlav());
        std::string tar = fParams.GetDetector(fParams.GetCurrentDet()).GetTarget();
        std::string curr = fParams.GetXSecName(fParams.GetCurrentXSec());

        if(curr.compare("NoXSec")) { // The comparison returns 0 if the two strings match...
          fXSecSplines[xsecname] = xsec->GetXSec(pdg, tar, curr); // Create and store the spline
        }
        else {
          // For NoXSec, or no cross section, create a line at 1
          TF1* f = new TF1("f", "1", 0., 120.);
          TSpline3* s = new TSpline3("", 0., 120., f, 120);
          fXSecSplines[xsecname] = s;
        }
      }
    }

    delete xsec; // Clean up

    return;
  }

  //---------------------------------------------------------------------------
  std::string Spectra::XSecName()
  {
    // Create a label by conglomerating the three strings
    std::string ret =   fParams.GetNuFlav(  fParams.GetCurrentNuFlav()).GetName()
                      + fParams.GetXSecName(fParams.GetCurrentXSec())
                      + fParams.GetDetName( fParams.GetCurrentDet());

    return ret;
  }
}
