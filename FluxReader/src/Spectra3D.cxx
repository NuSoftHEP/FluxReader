#include "Spectra3D.h"

// C/C++ Includes
#include <cassert>
#include <iostream>

// Root Includes
#include "TDirectory.h"
#include "TH1.h"
#include "TH3D.h"
#include "TSpline.h"

namespace flxrd
{
  //---------------------------------------------------------------------------
  Spectra3D::Spectra3D(Parameters params, std::string title,
                       std::string labelx, std::vector<double> binsx, const Var& varx,
                       std::string labely, std::vector<double> binsy, const Var& vary,
                       std::string labelz, std::vector<double> binsz, const Var& varz,
                       const Weight& wei, TObject* extWeights)
    : Spectra(params, title, varx, wei, extWeights), fVarY(vary), fVarZ(varz)
  {
    // Variables required by x axis variable and weight are set by Spectra constructor above
    // Add variables required by y and z axes variable to the list of branches
    fBranches.insert(fVarY.Branches().begin(), fVarY.Branches().end());
    fBranches.insert(fVarZ.Branches().begin(), fVarZ.Branches().end());

    CreateHists(labelx, binsx, labely, binsy, labelz, binsz);
  }

  //---------------------------------------------------------------------------
  TH1* Spectra3D::GetHist(int i_hist)
  {
    const int n_hist = fHists.size();
    if(i_hist < 0 || i_hist >= n_hist) {
      std::cout << "Input histogram index is out of range." << std::endl;
      assert(false);
    }

    return fHists[i_hist];
  }

  //---------------------------------------------------------------------------
  void Spectra3D::Fill(bsim::Dk2Nu* nu, std::map<std::string, int> nurayIndices)
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

    for(int i_det = 0, n_det = fParams.NDet(); i_det < n_det; ++i_det) {
      fParams.SetCurrentDet(i_det);

      int first_nuray = nurayIndices[fParams.GetDetName(i_det)];
      int last_nuray  = first_nuray + fParams.GetDetector(i_det).GetUses();
      if(last_nuray == first_nuray) {
        ++last_nuray;
      }

      for(int i_xsec = 0, n_xsec = fParams.NXSec(); i_xsec < n_xsec; ++i_xsec) {
        fParams.SetCurrentXSec(i_xsec);

        int i_hist = fParams.GetCurrentMaster();

        for(int i_nuray = first_nuray; i_nuray < last_nuray; ++i_nuray) {
          double weight =   nu->decay.nimpwt * nu->nuray[i_nuray].wgt
                          * fXSecSplines[XSecName()]->Eval(nu->nuray[i_nuray].E)
                          * fDefaultWeightCorrection;

          fHists[i_hist]->Fill(fVarX(nu, i_nuray), fVarY(nu, i_nuray), fVarZ(nu, i_nuray), fWei(weight, nu, i_nuray, fExtWeights));
        } // Loop over uses
      } // Loop over cross sections
    } // Loop over detectors

    return;
  }

  //---------------------------------------------------------------------------
  void Spectra3D::WriteHists(TDirectory* out)
  {
    TDirectory* temp = gDirectory;

    std::string det_name = "";

    for(const auto& index : fParams) {
      fParams.SetIndices(index);

      if(det_name.compare(fParams.GetDetName(fParams.GetCurrentDet()))) {
        det_name = fParams.GetDetName(fParams.GetCurrentDet());
        out->cd();

        if( !(out->GetListOfKeys()->FindObject(det_name.c_str())) ) {
          out->mkdir(fParams.GetDetName(fParams.GetCurrentDet()).c_str());
        }

        out->cd(fParams.GetDetName(fParams.GetCurrentDet()).c_str());
      }

      gDirectory->WriteTObject(fHists[index]);
    }

    temp->cd();
    return;
  }

  //---------------------------------------------------------------------------
  void Spectra3D::CreateHists(std::string labelx, std::vector<double> binsx,
                              std::string labely, std::vector<double> binsy,
                              std::string labelz, std::vector<double> binsz)
  {
    std::string hist_title = "";
    std::string axis_label = ";"+labelx+";"+labely+";"+labelz; // This is all axes labels
    const int nBinsX = binsx.size() - 1;
    const int nBinsY = binsy.size() - 1;
    const int nBinsZ = binsz.size() - 1;

    for(int i = 0, n = fParams.MaxMaster(); i < n; ++i)
    {
      hist_title = fTitle + "_" + fParams.NameTag(i);
      TH3D* h3 = new TH3D(hist_title.c_str(), axis_label.c_str(),
                          nBinsX, &binsx[0], nBinsY, &binsy[0], nBinsZ, &binsz[0]);
      fHists.push_back(h3);
    }

    return;
  }
}
