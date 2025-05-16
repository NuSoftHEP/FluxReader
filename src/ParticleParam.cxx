#include "ParticleParam.h"

namespace flxrd
{
  //---------------------------------------------------------------------------
  // Preset neutrino flavors
  const NuFlav NuFlav::kNue    = NuFlav("nue",    +12);
  const NuFlav NuFlav::kANue   = NuFlav("anue",   -12);
  const NuFlav NuFlav::kNumu   = NuFlav("numu",   +14);
  const NuFlav NuFlav::kANumu  = NuFlav("anumu",  -14);
  const NuFlav NuFlav::kNutau  = NuFlav("nutau",  +16);
  const NuFlav NuFlav::kANutau = NuFlav("anutau", -16);

  // Present neutrino parents
  const Parent Parent::kMuPlus  = Parent("muplus",  -13);
  const Parent Parent::kMuMinus = Parent("muminus", +13);
  const Parent Parent::kPiPlus  = Parent("piplus",  +211);
  const Parent Parent::kPiMinus = Parent("piminus", -211);
  const Parent Parent::kKPlus   = Parent("Kplus",   +321);
  const Parent Parent::kKMinus  = Parent("Kminus",  -321);

  // Preset parents used when PDG sign is ignored
  const Parent Parent::kMuon  = Parent("mu", 13);
  const Parent Parent::kPion  = Parent("pi", 211);
  const Parent Parent::kKaon  = Parent("K",  321);
  const Parent Parent::kKLong = Parent("KL", 130);

  //---------------------------------------------------------------------------
  ParticleParam::ParticleParam(std::string name, int pdg)
    : fName(name), fPDG(pdg) // Populate name and PDG parameters
  {}

  //---------------------------------------------------------------------------
  bool operator==(const ParticleParam& pA, const ParticleParam& pB)
  {
    return((pA.GetPDG() == pB.GetPDG()) &&
           !pA.GetName().compare(pB.GetName()));
  }

  //---------------------------------------------------------------------------
  std::vector<NuFlav> NuFlav::AllNuFlavs(bool SignSensitive)
  {
    std::vector<NuFlav> ret;

    // Add neutrinos to vector
    ret.push_back(kNue);
    if(SignSensitive) ret.push_back(kANue);
    ret.push_back(kNumu);
    if(SignSensitive) ret.push_back(kANumu);
    ret.push_back(kNutau);
    if(SignSensitive) ret.push_back(kANutau);

    return ret;
  }

  //---------------------------------------------------------------------------
  std::vector<Parent> Parent::AllParents(bool SignSensitive)
  {
    std::vector<Parent> ret;

    if(SignSensitive) {
      // Add parents to vector, discerning PDG sign
      ret.push_back(kMuPlus);
      ret.push_back(kMuMinus);
      ret.push_back(kPiPlus);
      ret.push_back(kPiMinus);
      ret.push_back(kKPlus);
      ret.push_back(kKMinus);
    }
    else {
      // Add parents to vector, ignoring PDG sign
      ret.push_back(kMuon);
      ret.push_back(kPion);
      ret.push_back(kKaon);
    }

    // Add K_L to vector
    ret.push_back(kKLong);

    return ret;
  }

  //---------------------------------------------------------------------------
  void NuFlav::RemoveNuFlav(std::vector<NuFlav> &nuflavs, int rmpdg)
  {
    unsigned int n_flav = nuflavs.size();
    for(unsigned int i_flav = 0; i_flav < n_flav; ) { // Iterator is incremented inside the loop
      if(nuflavs[i_flav].GetPDG() == rmpdg) {
        nuflavs.erase(nuflavs.begin() + i_flav); // Remove the neutrino flavor
        --n_flav; // Lower the number of neutrino flavors accordingly
      }
      else {
        ++i_flav; // Move to the next NuFlav entry
      }
    }

    return;
  }

  //---------------------------------------------------------------------------
  void NuFlav::RemoveNuFlav(std::vector<NuFlav> &nuflavs, std::string rmname)
  {
    // Similar to the call with an integer input
    unsigned int n_flav = nuflavs.size();
    for(unsigned int i_flav = 0; i_flav < n_flav; ) {
      if(!nuflavs[i_flav].GetName().compare(rmname)) {
        nuflavs.erase(nuflavs.begin() + i_flav);
        --n_flav;
      }
      else {
        ++i_flav;
      }
    }

    return;
  }

  //---------------------------------------------------------------------------
  void NuFlav::RemoveNuFlav(std::vector<NuFlav> &nuflavs, const NuFlav& rmflav)
  {
    // Similar to the call with an integer input
    unsigned int n_flav = nuflavs.size();
    for(unsigned int i_flav = 0; i_flav < n_flav; ) {
      if(nuflavs[i_flav] == rmflav) {
        nuflavs.erase(nuflavs.begin() + i_flav);
        --n_flav;
      }
      else {
        ++i_flav;
      }
    }

    return;
  }

  //---------------------------------------------------------------------------
  void Parent::RemoveParent(std::vector<Parent> &parents, int rmpdg)
  {
    // Similar to the RemoveNuFlav version
    unsigned int n_par = parents.size();
    for(unsigned int i_par = 0; i_par < n_par; ) {
      if(parents[i_par].GetPDG() == rmpdg) {
        parents.erase(parents.begin() + i_par);
        --n_par;
      }
      else {
        ++i_par;
      }
    }

    return;
  }

  //---------------------------------------------------------------------------
  void Parent::RemoveParent(std::vector<Parent> &parents, std::string rmname)
  {
    // Similar to the RemoveNuFlav version
    unsigned int n_par = parents.size();
    for(unsigned int i_par = 0; i_par < n_par; ) {
      if(!parents[i_par].GetName().compare(rmname)) {
        parents.erase(parents.begin() + i_par);
        --n_par;
      }
      else {
        ++i_par;
      }
    }

    return;
  }

  //---------------------------------------------------------------------------
  void Parent::RemoveParent(std::vector<Parent> &parents, const Parent& rmpar)
  {
    // Similar to the RemoveNuFlav version
    unsigned int n_par = parents.size();
    for(unsigned int i_par = 0; i_par < n_par; ) {
      if(parents[i_par] == rmpar) {
        parents.erase(parents.begin() + i_par);
        --n_par;
      }
      else {
        ++i_par;
      }
    }

    return;
  }
}
