#include "Parameters.h"

// C/C++ Includes
#include <iostream>

// Package Includes
#include "XSec.h"

namespace flxrd
{
  //---------------------------------------------------------------------------
  int Indices::GetCurrentMaster() const
  {
    // Standard convention "definition"
    int master =   iFlav
                 + iPar *nFlav
                 + iXSec*nPar *nFlav
                 + iDet *nXSec*nPar *nFlav;

    return master;
  }

  //---------------------------------------------------------------------------
  bool Indices::operator!=(const Indices& indices) const
  {
    return (indices.GetCurrentMaster() != GetCurrentMaster());
  }

  //---------------------------------------------------------------------------
  int Indices::operator*() const
  {
    return GetCurrentMaster();
  }

  //---------------------------------------------------------------------------
  const Indices& Indices::operator++()
  {
    // Do nothing if the current master is already at its max
    if(iDet >= nDet) {
      return *this;
    }

    ++iFlav; // Start by trying to increment iFlav

    // Set iFlav back to 0 if it is "out of bounds",
    // then try incrementing iPar
    if(iFlav >= nFlav) {
      iFlav = 0;
      ++iPar;
    }

    // And so on...
    if(iPar >= nPar) {
      iPar = 0;
      ++iXSec;
    }

    if(iXSec >= nXSec) {
      iXSec = 0;
      ++iDet;
    }

    return *this;
  }

  //---------------------------------------------------------------------------
  Parameters::Parameters(bool SignSensitive, bool verbosity)
    : fAncestorPar(true),
      fSignSensitive(SignSensitive),
      fVerbosity(verbosity)
  {
    SetDefaults(SignSensitive); // Set sensible default flavors, parents, and cross sections

    RemoveNuTaus(); // By default, remove tau neutrinos
  }

  //---------------------------------------------------------------------------
  Parameters::Parameters(const Parameters& params)
    : fAncestorPar(params.fAncestorPar),
      fSignSensitive(params.IsSignSensitive()),
      fVerbosity(params.fVerbosity),
      fNuFlav(params.fNuFlav),
      fParent(params.fParent),
      fXSec  (params.fXSec),
      fDet   (params.fDet),
      fIndices(params.fIndices)
  {
    fIndices.iFlav = 0;
    fIndices.iPar  = 0;
    fIndices.iXSec = 0;
    fIndices.iDet  = 0;
  }

  //---------------------------------------------------------------------------
  void Parameters::AddDetector(const Detector& det)
  {
    // Make sure the detector does not already exist
    for(unsigned int i_det = 0, n_det = fDet.size(); i_det < n_det; ++i_det) {
      if(!fDet[i_det].GetDetName().compare(det.GetDetName())) {
        std::cout << "This detector name already exists for another detector." << std::endl;
        return;
      }
    }

    fDet.push_back(det); // Add the detector
    UpdateIndices(); // Make sure the Parameters' Indices object is aware of this change
    return;
  }

  //---------------------------------------------------------------------------
  void Parameters::AddParent(Parent parent)
  {
    // Make sure the parent is not being duplicated
    for(unsigned int i_par = 0, n_par = fParent.size(); i_par < n_par; ++i_par) {
      if(parent.GetPDG() == fParent[i_par].GetPDG()) {
        std::cout << "This parent is already included as " << parent.GetName() << "." << std::endl;
        return;
      }
    }

    fParent.push_back(parent); // Add parent to vector
    UpdateIndices(); // Make sure the Parameters' Indices object is aware of this change
    return;
  }

  //---------------------------------------------------------------------------
  void Parameters::AddXSec(std::string xsec)
  {
    // Check that the cross section does not already exist
    for(unsigned int i_xsec = 0, n_xsec = fXSec.size(); i_xsec < n_xsec; ++i_xsec) {
      if(!xsec.compare(fXSec[i_xsec])) {
        std::cout << "This cross section is already included." << std::endl;
        return;
      }
    }

    // Create XSec object to get a list of valid cross sections
    XSec* xs = new XSec();

    // If the cross section was invalid, show the user a list of valid inputs, clean up, and return
    if(!xs->IsValidProcess(xsec) && xsec.compare("NoXSec")) {
      std::cout << "The input cross section is not valid. The following are valid:" << std::endl;
      xs->ListIntTypes();
      delete xs;
      return;
    }

    fXSec.push_back(xsec); // Add cross section to vector

    delete xs; // Clean up
    UpdateIndices(); // Make sure the Parameters' Indices object is aware of this change
    return;
  }

  //---------------------------------------------------------------------------
  int Parameters::GetCurrentDet() const
  {
    return fIndices.iDet;
  }

  //---------------------------------------------------------------------------
  int Parameters::GetCurrentNuFlav() const
  {
    return fIndices.iFlav;
  }

  //---------------------------------------------------------------------------
  int Parameters::GetCurrentParent() const
  {
    return fIndices.iPar;
  }

  //---------------------------------------------------------------------------
  int Parameters::GetCurrentXSec() const
  {
    return fIndices.iXSec;
  }

  //---------------------------------------------------------------------------
  int Parameters::GetCurrentMaster() const
  {
    return fIndices.GetCurrentMaster();
  }

  //---------------------------------------------------------------------------
  Detector Parameters::GetDetector(int i_det) const
  {
    if(i_det < 0 || i_det >= NDet()) {
      std::cout << "Input index is out of range." << std::endl; 
      std::cout << "Returning first (0th) detector." << std::endl;
      return fDet[0];
    }

    return fDet[i_det];
  }

  //---------------------------------------------------------------------------
  NuFlav Parameters::GetNuFlav(int i_flav) const
  {
    if(i_flav < 0 || i_flav >= NFlav()) {
      std::cout << "Input index is out of range." << std::endl; 
      std::cout << "Returning first (0th) NuFlav." << std::endl;
      return fNuFlav[0];
    }

    return fNuFlav[i_flav];
  }

  //---------------------------------------------------------------------------
  std::string Parameters::GetDetName(int i_det) const
  {
    if(i_det < 0 || i_det >= NDet()) {
      std::cout << "Input index to GetDetName is out of range." << std::endl;
      return "";
    }

    return fDet[i_det].GetDetName();
  }

  //---------------------------------------------------------------------------
  int Parameters::GetNuFlavPDG(int i_flav) const
  {
    if(i_flav < 0 || i_flav >= NFlav()) {
      std::cout << "Input index to GetNuFlavPDG is out of range." << std::endl;
      return -1;
    }

    return fNuFlav[i_flav].GetPDG();
  }

  //---------------------------------------------------------------------------
  int Parameters::GetParentPDG(int i_par) const
  {
    if(i_par < 0 || i_par >= NPar()) {
      std::cout << "Input index to GetParentPDG is out of range." << std::endl;
      return -1;
    }

    return fParent[i_par].GetPDG();
  }

  //---------------------------------------------------------------------------
  std::string Parameters::GetXSecName(int i_xsec) const
  {
    if(i_xsec < 0 || i_xsec >= NXSec()) {
      std::cout << "Input index to GetXSecName is out of range." << std::endl;
      return "";
    }

    return fXSec[i_xsec];
  }

  //---------------------------------------------------------------------------
  int Parameters::MaxMaster() const
  {
    return NFlav()*NPar()*NXSec()*NDet();
  }

  //---------------------------------------------------------------------------
  int Parameters::MaxMaster(int i_det) const
  {
    return NFlav()*NPar()*NXSec()*(i_det+1);
  }

  //---------------------------------------------------------------------------
  void Parameters::RemoveDetector(std::string rmname)
  {
    // This code is similar to removal functions defined in ParticleParam
    unsigned int n_det = fDet.size();
    for(unsigned int i_det = 0; i_det < n_det; ) {
      if(!fDet[i_det].GetDetName().compare(rmname)) {
        fDet.erase(fDet.begin() + i_det);
        --n_det;
      }
      else {
        ++i_det;
      }
    }

    UpdateIndices(); // Make sure the Parameters' Indices object is aware of this change
    return;
  }

  //---------------------------------------------------------------------------
  void Parameters::RemoveNuFlav(int rmpdg)
  {
    NuFlav::RemoveNuFlav(fNuFlav, rmpdg); // Remove the actual NuFlav object
    UpdateIndices(); // Make sure the Parameters' Indices object is aware of this change
    return;
  }

  //---------------------------------------------------------------------------
  void Parameters::RemoveNuFlav(std::string rmname)
  {
    NuFlav::RemoveNuFlav(fNuFlav, rmname);
    UpdateIndices(); // Make sure the Parameters' Indices object is aware of this change
    return;
  }

  //---------------------------------------------------------------------------
  void Parameters::RemoveNuFlav(const NuFlav& rmflav)
  {
    NuFlav::RemoveNuFlav(fNuFlav, rmflav);
    UpdateIndices(); // Make sure the Parameters' Indices object is aware of this change
    return;
  }

  //---------------------------------------------------------------------------
  void Parameters::RemoveNuTaus()
  {
    // Specific calls to RemoveNuFlav
    RemoveNuFlav(+16);
    RemoveNuFlav(-16);

    UpdateIndices(); // Make sure the Parameters' Indices object is aware of these changes
    return;
  }

  //---------------------------------------------------------------------------
  void Parameters::RemoveParent(int rmpdg)
  {
    Parent::RemoveParent(fParent, rmpdg); // Remove the actual Parent object
    UpdateIndices(); // Make sure the Parameters' Indices object is aware of this change
    return;
  }

  //---------------------------------------------------------------------------
  void Parameters::RemoveParent(std::string rmname)
  {
    Parent::RemoveParent(fParent, rmname);
    UpdateIndices(); // Make sure the Parameters' Indices object is aware of this change
    return;
  }

  //---------------------------------------------------------------------------
  void Parameters::RemoveParent(const Parent& rmpar)
  {
    Parent::RemoveParent(fParent, rmpar);
    UpdateIndices(); // Make sure the Parameters' Indices object is aware of this change
    return;
  }

  //---------------------------------------------------------------------------
  void Parameters::RemoveXSec(std::string rmxsec)
  {
    // This code is similar to removal functions defined in ParticleParam
    unsigned int n_xsec = fXSec.size();
    for(unsigned int i_xsec = 0; i_xsec < n_xsec; ) {
      if(!fXSec[i_xsec].compare(rmxsec)) {
        fXSec.erase(fXSec.begin() + i_xsec);
        --n_xsec;
      }
      else {
        ++i_xsec;
      }
    }

    UpdateIndices(); // Make sure the Parameters' Indices object is aware of this change
    return;
  }

  //---------------------------------------------------------------------------
  void Parameters::ResetNuFlavs()
  {
    // Make sure vector is cleared
    if(!fNuFlav.empty()) {
      fNuFlav.clear();
    }

    fNuFlav = NuFlav::AllNuFlavs(); // Populate the vector with all flavors
    UpdateIndices(); // Make sure the Parameters' Indices object is aware of this change
    return;
  }

  //---------------------------------------------------------------------------
  void Parameters::SetAncestorPar()
  {
    fAncestorPar = true;
    return;
  }

  //---------------------------------------------------------------------------
  void Parameters::SetAncestorTgt()
  {
    fAncestorPar = false;
    return;
  }

  //---------------------------------------------------------------------------
  void Parameters::SetDefaults(bool SignSensitive)
  {
    fNuFlav = NuFlav::AllNuFlavs(); // All neutrino flavors
    fParent = Parent::AllParents(SignSensitive); // All parents

    if(!fXSec.empty()) {
      fXSec.clear();
    }
    // No cross section, CC, and NC
    fXSec.push_back("NoXSec");
    fXSec.push_back("tot_cc");
    fXSec.push_back("tot_nc");

    UpdateIndices(); // Make sure the Parameters' Indices object is aware of these changes
    return;
  }

  //---------------------------------------------------------------------------
  void Parameters::SetDetUses(std::string detname, int nuses)
  {
    for(unsigned int i_det = 0, n_det = fDet.size(); i_det < n_det; ++i_det) {
      if(!detname.compare(fDet[i_det].GetDetName())) { // Find the corresponding detector
        fDet[i_det].SetUses(nuses); // Set the number of uses
      }
    }

    return;
  }

  //---------------------------------------------------------------------------
  void Parameters::ClearAll()
  {
    fNuFlav.clear();
    fParent.clear();
    fXSec.clear();
    fDet.clear();

    fIndices.iFlav = 0;
    fIndices.iPar  = 0;
    fIndices.iXSec = 0;
    fIndices.iDet  = 0;

    return;
  }

  //---------------------------------------------------------------------------
  std::string Parameters::NameTag(int master)
  {
    if(!SetIndices(master)) { // Check that the index is within its proper bounds
      return "";
    }

    // Conglomerate all the strings, separated by underscores
    std::string ret = fNuFlav[fIndices.iFlav].GetName();
    ret += "_";
    ret += fParent[fIndices.iPar].GetName();
    ret += "_";
    ret += fXSec[fIndices.iXSec];
    ret += "_";
    ret += fDet[fIndices.iDet].GetDetName();

    return ret;
  }

  //---------------------------------------------------------------------------
  bool Parameters::SetCurrentDet(int i_det)
  {
    if(i_det < 0 || i_det >= NDet()) {
      std::cout << "Input detector index is out of range." << std::endl;
      return false;
    }

    fIndices.iDet = i_det;
    return true;
  }

  //---------------------------------------------------------------------------
  bool Parameters::SetCurrentNuFlav(int PDG)
  {
    for(int i_flav = 0; i_flav < NFlav(); ++i_flav) {
      if(PDG == fNuFlav[i_flav].GetPDG()) {
        fIndices.iFlav = i_flav;
        return true;
      }
    }

    if(fVerbosity) {
      std::cout << "Could not find " << PDG << " in flavor vector." << std::endl;
    }
    return false;
  }

  //---------------------------------------------------------------------------
  bool Parameters::SetCurrentParent(int PDG)
  {
    for(int i_par = 0; i_par < NPar(); ++i_par) {
      if(PDG == fParent[i_par].GetPDG()) {
        fIndices.iPar = i_par;
        return true;
      }
    }

    if(fVerbosity) {
      std::cout << "Could not find " << PDG << " in parent vector." << std::endl;
    }
    return false;
  }

  //---------------------------------------------------------------------------
  bool Parameters::SetCurrentXSec(int i_xsec)
  {
    if(i_xsec < 0 || i_xsec >= NXSec()) {
      std::cout << "Input cross section index is out of range." << std::endl;
      return false;
    }

    fIndices.iXSec = i_xsec;
    return true;
  }

  //---------------------------------------------------------------------------
  bool Parameters::SetIndices(int master)
  {
    // Since master is NOT passed by pointer or reference,
    // this function makes a copy of master

    if(master < 0 || master >= MaxMaster()) {
      std::cout << "Input master index is out of range." << std::endl;
      return false;
    }

    // Set the detector index and subtract the relevant value from the copy of master
    fIndices.iDet = (int)(master/(NFlav()*NPar()*NXSec()));
    master -= fIndices.iDet*NFlav()*NPar()*NXSec();

    // Set the cross section index and subtract the relevant value from the copy of master
    fIndices.iXSec = (int)(master/(NFlav()*NPar()));
    master -= fIndices.iXSec*NFlav()*NPar();

    // Set the parent index and subtract the relevant value from the copy of master
    fIndices.iPar  = (int)(master/NFlav());
    master -= fIndices.iPar*NFlav();

    // What remains of the copy of master corresponds to just the flavor index
    fIndices.iFlav = master;

    return true;
  }

  //---------------------------------------------------------------------------
  void Parameters::UpdateIndices()
  {
    fIndices.nFlav = fNuFlav.size();
    fIndices.nPar  = fParent.size();
    fIndices.nXSec = fXSec.size();
    fIndices.nDet  = fDet.size();
  }

  //---------------------------------------------------------------------------
  Indices Parameters::begin()
  {
    fIndices.iFlav = 0;
    fIndices.iPar  = 0;
    fIndices.iXSec = 0;
    fIndices.iDet  = 0;

    return fIndices;
  }

  //---------------------------------------------------------------------------
  Indices Parameters::end() const
  {
    Indices i;

    i.nFlav = fNuFlav.size();
    i.nPar  = fParent.size();
    i.nXSec = fXSec.size();
    i.nDet  = fDet.size();

    i.iFlav = 0;
    i.iPar  = 0;
    i.iXSec = 0;
    i.iDet  = fDet.size();

    return i;
  }
}
