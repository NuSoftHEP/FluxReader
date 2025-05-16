#pragma once

// C/C++ Includes
#include <map>
#include <set>
#include <string>
#include <vector>

// Root Includes
#include "TVector3.h"

// Package Includes
#include "Parameters.h"
#include "Var.h"
#include "Weight.h"

// Forward Class Definitions
class TBranch;
class TDirectory;
class TH1;
class TTree;

namespace bsim { class Dk2Nu;  }
namespace bsim { class DkMeta; }

namespace flxrd
{
  class Spectra;

  /// A class which reads flux files and outputs user defined histograms
  class FluxReader
  {
  public:
    /// \param fileWildcard A string (which can contain wildcard characters) that
    ///                     is a path to the input flux files
    /// \param numFiles The MAXIMUM number of files to run over
    ///                 (If the input wildcard returns less than numFiles, all will be used)
    /// \param skipFiles The number of files to skip
    FluxReader(std::string fileWildcard,
               unsigned int numFiles = 0,
               unsigned int skipFiles = 0);

    ~FluxReader();

    /// Loops through input files, populates histograms and writes them to file
    void ReadFlux(TDirectory* out); // FIX

    /// Add a Spectra(N)D object to populate
    void AddSpectra(Parameters params, std::string title,
                    std::string labelx, std::vector<double> binsx, const Var& varx,
                    const Weight& wei = kDefaultW, TObject* extWeights = nullptr);
    void AddSpectra(Parameters params, std::string title,
                    std::string labelx, std::vector<double> binsx, const Var& varx,
                    std::string labely, std::vector<double> binsy, const Var& vary,
                    const Weight& wei = kDefaultW, TObject* extWeights = nullptr);
    void AddSpectra(Parameters params, std::string title,
                    std::string labelx, std::vector<double> binsx, const Var& varx,
                    std::string labely, std::vector<double> binsy, const Var& vary,
                    std::string labelz, std::vector<double> binsz, const Var& varz,
                    const Weight& wei = kDefaultW, TObject* extWeights = nullptr);
    void AddSpectra(Parameters params, std::string title,
                    std::string detX, std::string detY,
                    std::string labelx, std::vector<double> binsx, const Var& varx,
                    const Weight& wei = kDefaultW, TObject* extWeights = nullptr);

    /// Allow FluxReader to read files that are not standard Dk2Nu files
    void OverrideTreeName(std::string treepath);
    void OverridePOTPath(std::string metapath, std::string potpath);
    void OverrideDefaultVarName(std::string oldname, std::string newname);

  private:
    /// Add branch(es) to the master list of branches to turn on
    void AddBranch(std::string branchName);
    void AddBranches(std::set<std::string> branchNames);

    /// Add the pre-defined list of branches to the master list
    void AddDefaultBranches();

    /// Notify the user of the parameters to be run over
    void InitialMessage();

    /// Check if any standard Dk2Nu variable names have been overriden
    bool IsStandardDk2Nu();

    /// Set necessary addresses for entries in the flux and metadata trees 
    void SetBranches(TTree* fluxTree, TTree* metaTree);

    /// Set up the map which points a detector name to its first index in the Dk2Nu object's NuRay vector
    void SetNuRayIndices();

    /// Randomly pick a point somewhere in the detector
    /// \param rr For a non square detector, pick a point such that
    ///           x*x + y*y is less than \a rr
    TVector3 Smear(const Detector& det, double rr = -1); // ISSUE

    /// Convert from detector coordinates to beam coordinates
    void ToBeamCoords(const Detector& det, TVector3& xyz);

    std::vector<TBranch*> fBranches; ///< List of TBranches to activate
    std::set<std::string> fBranchNames; ///< List of branch names that will be activated

    std::map<std::string, std::string> fBranchOverrides; ///< Point default Dk2Nu branch name to non-standard branch name

    std::set<Detector> fDetectors; ///< List of detectors to point neutrino rays toward

    std::vector<std::string> fInputFiles; ///< List of input files to run over

    bsim::Dk2Nu*  fNu;   ///< Dk2Nu object that will store values for each input file entry
    bsim::DkMeta* fMeta; ///< DkMeta object that will store metadata about the Dk2Nu tree

    std::map<std::string, int> fNuRayIndex; ///< Map pointing from a detector name to its first index in the Dk2Nu object's NuRay vector

    bool fReweightNuRay; ///< Helper to determine whether neutrino rays need to be reweighted

    /// Spectra vector
    /// All relevant functions are declared in the abstract Spectra class,
    /// so this vector can handle any dimensional Spectra object pointer
    std::vector<Spectra*> fSpectra;

    std::string fTreePath; ///< Path that points to the tree in each input file
    std::string fMetaPath; ///< Path that points to the metadata tree of an input file
    std::string fPOTPath;  ///< Path that points to the POT variable in the metadata tree
  };
}
