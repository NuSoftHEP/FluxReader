#pragma once

// Root Includes
#include "TH1.h"
#include "TMath.h"

// Package Includes
#include "Weight.h"

namespace flxrd
{
  // This file contains common weights

  // The Weight object closely resembles the Var object.
  // Refer to the documentation in Vars.h for a more detailed description
  // of the format for a Var/Weight.

  // Weight by the standard input, w,
  // multiplied by an external weight identified by the neutrino parent pT and pz
  const Weight kExtWeightBypTpz({"ancestor","ancestor.stoppx","ancestor.stoppy","ancestor.stoppz"},
                 [](const double& w, const bsim::Dk2Nu* nu, const int& i_nuray, const TObject* extW)
                 { TH1* hExtW = (TH1*)extW;
                   double px = nu->ancestor[0].stoppx;
                   double py = nu->ancestor[0].stoppy;
                   double pz = nu->ancestor[0].stoppz;
                   int bin = hExtW->FindFixBin(sqrt(px*px + py*py), pz);
                   if(bin == -1)
                     return 0.;
                   return w*hExtW->GetBinContent(bin); });
}
