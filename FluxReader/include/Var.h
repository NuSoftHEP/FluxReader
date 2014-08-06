#pragma once

// C/C++ Includes
#include <functional>
#include <set>
#include <string>

// Other External Includes
#include "dk2nu.h"

namespace flxrd
{
  /// \brief A class which represents a variable used to bin events
  ///
  /// A Var takes a list of variables that need to be read from a flux file,
  /// and a function which determines how the variable is calculated
  /// See Vars.h for common variables
  class Var
  {
  public:
    /// This is the standard format for the Var function
    /// The Dk2Nu object stores all necessary values for a given entry
    /// The i_nuray integer points to the necessary index in the Dk2Nu NuRay vector
    typedef double (VarFunc_t)(const bsim::Dk2Nu* nu, const int& i_nuray);

    Var(const std::set<std::string>& branches,
        const std::function<VarFunc_t>& func)
      : fBranches(branches), fFunc(func) {}

    /// Copy constructor
    Var(const Var& copy) : fBranches(copy.fBranches), fFunc(copy.fFunc) {}

    /// Return the list of branches needed for the Var
    const std::set<std::string>& Branches() const { return fBranches; }

    /// Allow the Var to be called as a function, i.e., var(nu, i_nuray)
    double operator()(const bsim::Dk2Nu* nu, const int& i_nuray)
    {
      return fFunc(nu, i_nuray);
    }

  protected:
    std::set<std::string> fBranches; ///< List of variable names needed from the input flux file
    std::function<VarFunc_t> fFunc; ///< The function to calculate the variable
  };
}
