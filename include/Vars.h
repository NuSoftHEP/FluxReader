#pragma once

// Root Includes
#include "TMath.h"

// Package Includes
#include "Var.h"

namespace flxrd
{
  // This file should contain common variables

  // Note the format of each Var:
  // (a)const (b)Var (c)kEnergy((d){"nuray", "nuray.E"},
  //                  (e)[](const bsim::Dk2Nu* nu, const int& (*)i_nuray)
  //                  (f){ return nu->nuray[i_nuray].E; });
  // (a) The const label means that the Var will not change once it is created
  //     i.e., every entry will be calculated the same way.
  // (b) This is just specifying the type of object, like int, or double.
  // (c) This is the name of the object, like 'x' in "int x = 3;".
  // (d) This is the list of branch names needed for the variable.
  // (e) This slightly opaque line is the type of function that will be used.
  //     to calculate the value for a given entry.
  //     (*)If the Var does not access the "nuray" branch, i_nuray can be omitted,
  //        see kpT or kpz for an example of this.
  // (f) The block of code inside the curly braces ("{", "}")
  //     determines what is returned, and how it is calculated.
  //     kEnergy is simple enough to have this done in one line,
  //     but other more complicated variables can span multiple lines.
  //     Simply separate lines by semicolons as in regular code,
  //     and see kpT or kpz for examples of blocks that use more than one line.

  // Neutrino energy
  const Var kEnergy({"nuray", "nuray.E"},
                    [](const bsim::Dk2Nu* nu, const int& i_nuray)
                    { return nu->nuray[i_nuray].E; });

  // These are probably better incantations of pT and pz,
  // but the ancestor branch may not be filled.
  // Momentum of neutrino parent transverse to beam direction
  /*const Var kpT({"ancestor", "ancestor.stoppx", "ancestor.stoppy"},
                [](const bsim::Dk2Nu* nu, const int&)
                { double px = nu->ancestor[0].stoppx;
                  double py = nu->ancestor[0].stoppy;
                  return sqrt(px*px + py*py); });

  // Momentum of neutrino parent along beam direction
  const Var kpz({"ancestor", "ancestor.stoppz"},
                [](const bsim::Dk2Nu* nu, const int&)
                { return nu->ancestor[0].stoppz; });*/

  // Momentum of neutrino parent transverse to beam direction
  const Var kpT({"decay", "decay.pdpx", "decay.pdpy"},
                [](const bsim::Dk2Nu* nu, const int&)
                { double px = nu->decay.pdpx;
                  double py = nu->decay.pdpy;
                  return sqrt(px*px + py*py); });

  // Momentum of neutrino parent along beam direction
  const Var kpz({"decay", "decay.pdpz"},
                [](const bsim::Dk2Nu* nu, const int&)
                { return nu->decay.pdpz; });

  // Momentum of neutrino ancestor (not necessarily parent) transverse to beam direction,
  // as it leaves the NuMI Target
  const Var kTargetExitpT({"tgtexit", "tgtexit.tpx", "tgtexit.tpy"},
                          [](const bsim::Dk2Nu* nu, const int&)
                          { double px = nu->tgtexit.tpx;
                          double py = nu->tgtexit.tpy;
                          return sqrt(px*px + py*py); });

  // Momentum of neutrino ancestor (not necessarily parent) along beam direction,
  // as it leaves the NuMI Target
  const Var kTargetExitpz({"tgtexit", "tgtexit.tpz"},
                          [](const bsim::Dk2Nu* nu, const int&)
                          { return nu->tgtexit.tpz; });
}
