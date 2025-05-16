#pragma once

// C/C++ Includes
#include <functional>
#include <set>
#include <string>

// Other External Includes
#include "dk2nu.h"

// Forward Class Definitions
class TObject;

namespace flxrd
{
  /// \brief A class which represents a weight applied to neutrino events
  ///
  /// A weight takes a list of variables that need to be read from a flux file,
  /// and a function which determines how the weight is calculated
  class Weight
  {
  public:
    /// This is the standard format for the Weight function
    /// The double input represents a simple weight commonly applied to all events
    /// The Dk2Nu object stores all necessary values for a given entry
    /// The i_nuray integer points to the necessary index in the Dk2Nu NuRay vector
    /// The TObject* pointer is a set of weights calculated externally from the FluxReader package
    typedef double (WeiFunc_t)(const double& w, const bsim::Dk2Nu* nu, const int& i_nuray, const TObject* extW);

    Weight(const std::set<std::string>& branches,
           const std::function<WeiFunc_t>& func)
      : fBranches(branches), fFunc(func) {}

    /// Copy constructor
    Weight(const Weight& copy) : fBranches(copy.fBranches), fFunc(copy.fFunc) {}

    /// Return the list of branches needed for the Weight
    const std::set<std::string>& Branches() const { return fBranches; }

    /// Allow the Weight to be called as a function, i.e., wei(w, nu, i_nuray, extW)
    double operator()(const double& w, const bsim::Dk2Nu* nu, const int& i_nuray, const TObject* extW)
    {
      return fFunc(w, nu, i_nuray, extW);
    }

  protected:
    std::set<std::string> fBranches; ///< List of branch names needed from the input flux file
    std::function<WeiFunc_t> fFunc; ///< The function to calculate the weight
  };

  /// All entries get weighted by 'importance weight'*'propagation weight'*'cross section'
  /// Instead of having each Weight need to include these branches in their lists,
  /// this default weight takes the input double as the weight.
  /// This allows the product above to be passed as an input
  const Weight kDefaultW({}, [](const double& w, const bsim::Dk2Nu*, const int&, const TObject*)
                         { return w; });

  /// All entries have weight 1
  const Weight kNoWeight({}, [](const double&, const bsim::Dk2Nu*, const int&, const TObject*)
                         { return 1.; });

  /// All entries have weight c
  const Weight kConstant(double c);
}
