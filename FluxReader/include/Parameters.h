#pragma once

// C/C++ Includes
#include <string>
#include <vector>

// Package Includes
#include "Detector.h"
#include "ParticleParam.h"

namespace flxrd
{
  class Parameters; // Forward definition for reference in Indices

  /// A class that points to a specific set of parameters and can increment through them
  /// It stores index for neutrino flavor, parent, cross section, and detector
  /// These indices can be combined to form a "master" index
  /// By analogy to numbers written in powers of ten,
  /// think of iFlav as the 'units digit,' iPar as the 'tens digit,'
  /// iXSec as the 'hundreds digit,' and iDet as the 'thousands digit'
  /// Whenever iFlav reaches nFlav, it is reset to 0 and iPar is incremented, and so on
  /// With zero-indexing, the maximum value is when every index satisfies i<Param> = n<Param>-1
  class Indices
  {
  public:
    friend class Parameters; // Allow Parameters to change the private fields

    Indices() : iFlav(0), iPar(0), iXSec(0), iDet(0),
                nFlav(0), nPar(0), nXSec(0), nDet(0) {}

    /// Calculate and return the current master index
    int GetCurrentMaster() const;

    /// Compare a different Indices object to 'this' one by comparing current masters
    bool operator!=(const Indices& indices) const;

    /// Dereference operator: return current master
    int operator*() const;

    /// Increment operator
    /// This operator MUST PRECEED the object; i.e., ++IndicesObj, not IndicesObj++
    /// This increments the stored indices as appropriate, if possible
    const Indices& operator++();

  private:
    int iFlav;
    int iPar;
    int iXSec;
    int iDet;

    int nFlav;
    int nPar;
    int nXSec;
    int nDet;
  };

  /// A class that stores the parameters to apply to a FluxReader output file
  /// Includes neutrino flavors, neutrino parents, cross sections, and detectors
  class Parameters
  {
  public:
    friend class Combiner;
    friend class FluxReader;
    friend class Spectra;
    friend class Spectra1D;
    friend class Spectra2D;
    friend class Spectra3D;
    friend class SpectraCorrDet;

    /// This is the default constructor for users
    Parameters(bool SignSensitive = true, bool verbosity = true);

    /// Add a detector to run over
    void AddDetector(const Detector& det);

    /// Add a parent to run over
    void AddParent(Parent parent);

    /// Add a cross section to apply
    void AddXSec(std::string xsec);

    /// Access the internal index pointing to the current parameter tag
    int GetCurrentDet   () const;
    int GetCurrentNuFlav() const;
    int GetCurrentParent() const;
    int GetCurrentXSec  () const;

    /// Access the current master index
    int GetCurrentMaster() const;

    /// Get whether to split neutrino ancestor by parent or ancestor exiting target
    bool GetAncestorPar() const { return fAncestorPar; }

    /// Pull a stored detector to call its class functions
    Detector GetDetector(int i_det) const;

    /// Pull a stored NuFlav to call its class functions
    NuFlav GetNuFlav(  int i_flav) const;

    /// Shortcut to access necessary fields
    std::string GetDetName  (int i_det)  const;
    int         GetNuFlavPDG(int i_flav) const;
    int         GetParentPDG(int i_par)  const;
    std::string GetXSecName (int i_xsec) const;

    /// Get whether the sign of the neutrino parent is considered (true) or ignored (false)
    bool IsSignSensitive() const { return fSignSensitive; }

    /// Maximum master index over all parameters
    int MaxMaster() const;

    /// Maximum master index with fixed detector i_det
    int MaxMaster(int i_det) const;

    /// Number of each of the four parameter tags
    int NFlav() const { return fNuFlav.size(); }
    int NPar()  const { return fParent.size(); }
    int NXSec() const { return fXSec  .size(); }
    int NDet()  const { return fDet   .size(); }

    /// Remove detector to be used
    void RemoveDetector(std::string det);

    /// Remove a neutrino flavor to run over
    void RemoveNuFlav(int rmpdg);
    void RemoveNuFlav(std::string rmname);
    void RemoveNuFlav(const NuFlav& rmflav);

    /// Remove tau neutrinos. This is two pre-defined calls to RemoveNuFlav
    void RemoveNuTaus();

    /// Remove neutrino parent to run over
    void RemoveParent(int rmpdg);
    void RemoveParent(std::string rmname);
    void RemoveParent(const Parent& rmpar);

    /// Remove cross section to be applied
    void RemoveXSec(std::string rmxsec);

    /// Reset fNuFlavs to include all neutrino flavors
    void ResetNuFlavs();

    /// Set Parameters to split by neutrino parent
    void SetAncestorPar();

    /// Set Parameters to split by ancestor leaving target
    void SetAncestorTgt();

    /// Set parameter vectors to some sensible defaults
    /// This does not add any detectors
    void SetDefaults(bool SignSensitive = true);

    /// Set the number of uses for a specific detector
    void SetDetUses(std::string detname, int nuses);

  private:
    /// Copy constructor
    /// This is only used by Spectra during its own construction
    /// It copies the current Parameters into a new object stored by the Spectra object,
    /// allowing for the original one to be configured further
    Parameters(const Parameters& params);

    /// Remove all parameters
    void ClearAll();

    /// Creates a string corresponding to the first 4 index labels,
    /// also sets indices to the input master index
    std::string NameTag(int master);

    /// Set an internal index to the value provided
    bool SetCurrentDet   (int i_det);
    bool SetCurrentNuFlav(int PDG);
    bool SetCurrentParent(int PDG);
    bool SetCurrentXSec  (int i_xsec);

    /// Given a master index, set the internal indices to match
    bool SetIndices(int master);

    /// Make sure the Parameters fIndices object is aware of any parameter additions/removals
    void UpdateIndices();

    /// For looping using Indices
    /// This returns an Indices object with a master index of 0
    Indices begin();

    /// For looping using Indices
    /// This returns an Indices object with a master index of "maximum + 1"
    /// (See the Indices class documentation for the definition of "maximum")
    /// The value "maximum + 1" is analogous to many size or length functions,
    /// which can be thought of as returning "index of last element + 1"
    Indices end() const;

    bool fAncestorPar;   ///< Store whether neutrino plots are split by parent or target exit ancestor
    bool fSignSensitive; ///< Store whether neutrino sign is considered (true) or ignored (false)
    bool fVerbosity;     ///< Determine how much output to print

    /// Parameter tag vectors
    std::vector<NuFlav>      fNuFlav;
    std::vector<Parent>      fParent;
    std::vector<std::string> fXSec;
    std::vector<Detector>    fDet;

    Indices fIndices; ///< These are the actual internal indices for the Parameters object
  };
}
