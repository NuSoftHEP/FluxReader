#pragma once

// C/C++ Includes
#include <map>
#include <set>
#include <string>
#include <vector>

// Forward Class Definitions
namespace bsim { class Dk2Nu; }

namespace flxrd
{
  /// Make a vector of equally spaced bin edges from a min value, max value, and number of bins
  std::vector<double> Bins(int nbins, double min, double max);

  /// Create a map of Dk2Nu branch name labels pointing to corresponding locations in a Dk2Nu object
  std::map<std::string, void*> OverrideAddresses(bsim::Dk2Nu* nu);

  /// Expand an input wildcard string into a vector of all file names matching the pattern
  std::vector<std::string> Wildcard(std::string fileWildcard);
}
