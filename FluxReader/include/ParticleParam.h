#pragma once

// C/C++ Includes
#include <string>
#include <vector>

namespace flxrd
{
  /// Class containing a particle name and PDG
  class ParticleParam
  {
  public:
    ParticleParam(std::string name, int pdg);

    /// Get particle name
    std::string GetName() const { return fName; }

    /// Get particle PDG
    int         GetPDG()  const { return fPDG; }

  protected:
    std::string fName; ///< Particle name
    int         fPDG;  ///< Particle PDG
  };

  /// Operator to compare two ParticleParam objects
  /// Compares name and PDG
  bool operator==(const ParticleParam& pA, const ParticleParam& pB);

  /// A ParticleParam that is specificly a neutrino flavor
  class NuFlav: public ParticleParam
  {
  public:
    NuFlav(std::string name, int pdg) : ParticleParam(name, pdg) {}

    /// Preset NuFlavs
    static const NuFlav kNue;
    static const NuFlav kANue;
    static const NuFlav kNumu;
    static const NuFlav kANumu;
    static const NuFlav kNutau;
    static const NuFlav kANutau;

    /// Returns a vector of all neutrino flavors
    static std::vector<NuFlav> AllNuFlavs(bool SignSensitive = true);

    /// Remove a NuFlav from a vector of NuFlavs by PDG
    static void RemoveNuFlav(std::vector<NuFlav> &nuflavs, int rmpdg);

    /// Remove a NuFlav from a vector of NuFlavs by name
    static void RemoveNuFlav(std::vector<NuFlav> &nuflavs, std::string rmname);

    /// Remove a NuFlav from a vector of NuFlavs by NuFlav object
    static void RemoveNuFlav(std::vector<NuFlav> &nuflavs, const NuFlav& rmflav);
  };

  /// A ParticleParam that is specificly a neutrino parent
  class Parent: public ParticleParam
  {
  public:
    Parent(std::string name, int pdg) : ParticleParam(name, pdg) {}

    /// Present Parents
    static const Parent kMuPlus;
    static const Parent kMuMinus;
    static const Parent kPiPlus;
    static const Parent kPiMinus;
    static const Parent kKPlus;
    static const Parent kKMinus;

    static const Parent kMuon;
    static const Parent kPion;
    static const Parent kKaon;
    static const Parent kKLong;

    /// Returns a vector of all neutrino parents
    static std::vector<Parent> AllParents(bool SignSensitive = true);

    /// Remove a Parent from a vector of Parents by PDG
    static void RemoveParent(std::vector<Parent> &parents, int rmpdg);

    /// Remove a Parent from a vector of Parents by name
    static void RemoveParent(std::vector<Parent> &parents, std::string rmname);

    /// Remove a Parent from a vector of Parents by Parent object
    static void RemoveParent(std::vector<Parent> &parents, const Parent& rmpar);
  };
}
