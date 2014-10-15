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
  Detector::Detector(const std::string& det_name, const std::string& target,
                     const std::vector<double>& coords,
                     const std::vector<double>& sizes,
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
    // Make sure the input vector has the right size
    if(coords.size() == 3) {
      fCoord.push_back(coords[0]);
      fCoord.push_back(coords[1]);
      fCoord.push_back(coords[2]);
    }
    else {
      fCoord.push_back(0.);
      fCoord.push_back(0.);
      fCoord.push_back(0.);
    }

    // Populate size vector
    // Make sure the input vector has the right size
    if(sizes.size() == 3) {
      fSize.push_back(sizes[0]);
      fSize.push_back(sizes[1]);
      fSize.push_back(sizes[2]);
    }
    else {
      fSize.push_back(0.);
      fSize.push_back(0.);
      fSize.push_back(0.);
    }
  }

  //---------------------------------------------------------------------------
  TVector3 Detector::GetTCoords() const
  {
    return TVector3(fCoord[0], fCoord[1], fCoord[2]);
  }
}
