#pragma once

// C/C++ Includes
#include <string>
#include <vector>

// Root Includes
#include "TVector3.h"

namespace flxrd
{
  /// Class encoding information about detectors
  class Detector
  {
  public:
    /// \param det_name The name of the detector; this will be what is saved to file
    Detector(const std::string& det_name, const std::string& target,
             const double& coordx, const double& coordy, const double& coordz,
             const double& sizex,  const double& sizey,  const double& sizez,
             const int& nuses);

    Detector(const std::string& det_name, const std::string& target,
             const std::vector<double>& coords,
             const std::vector<double>& sizes,
             const int& nuses);

    /// Get detector name
    std::string GetDetName() const { return fDetName; }

    /// Get the detector target nucleus type
    std::string GetTarget() const { return fTarget; }

    /// Get detector coordinate(s)
    double GetCoordX() const { return fCoord[0]; }
    double GetCoordY() const { return fCoord[1]; }
    double GetCoordZ() const { return fCoord[2]; }
    std::vector<double> GetCoords() const { return fCoord; }

    /// Get detector coordinates as a TVector3
    TVector3 GetTCoords() const;

    /// Get detector size(s)
    double GetSizeX() const { return fSize[0]; }
    double GetSizeY() const { return fSize[1]; }
    double GetSizeZ() const { return fSize[2]; }
    std::vector<double> GetSizes() const { return fSize; }

    /// Preset functions for returning half of the detector size
    double GetHalfSizeX() const { return fSize[0]/2.; }
    double GetHalfSizeY() const { return fSize[1]/2.; }
    double GetHalfSizeZ() const { return fSize[2]/2.; }

    /// Get the number of times to use a neutrino ray in the detector
    int GetUses() const { return fUses; }

    /// Set the number fo times to use a neutrino ray in the detector
    void SetUses(int nuses) { fUses = nuses; }

    /// Simple method for comparing detectors, based on name comparison
    bool operator<(const Detector& otherDet) const
    {
      return (fDetName.compare(otherDet.fDetName) < 0);
    }

  private:
    std::string fDetName;       ///< Detector name
    std::string fTarget;        ///< Detector target nucleus
    std::vector<double> fCoord; ///< Detector coordinates (cm)
    std::vector<double> fSize;  ///< Detector size (cm)

    int fUses; ///< Number of times to use a neutrino ray in the detector
  };
}
