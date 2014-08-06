#include "Detector.h"

namespace flxrd
{
  //---------------------------------------------------------------------------
  Detector::Detector(const std::string& det_name, const std::string& target,
                     const double& coordx, const double& coordy, const double& coordz,
                     const double& sizex,  const double& sizey,  const double& sizez,
                     const int& nuses)
    : fDetName(det_name), fTarget(target), fUses(nuses)
  {
    // Make sure vectors are empty
    if(!fCoord.empty()) {
      fCoord.clear();
    }
    if(!fSize.empty()) {
      fSize.clear();
    }

    // Populate coordinate vector
    fCoord.push_back(coordx);
    fCoord.push_back(coordy);
    fCoord.push_back(coordz);

    // Populate size vector
    fSize.push_back(sizex);
    fSize.push_back(sizey);
    fSize.push_back(sizez);
  }

  //---------------------------------------------------------------------------
  TVector3 Detector::GetTCoords() const
  {
    return TVector3(fCoord[0], fCoord[1], fCoord[2]);
  }
}
