#include "FluxReader.h"

// C/C++ Includes
#include <cassert>
#include <iostream>
#include <utility>

// Root Includes
#include "TBranch.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TObject.h"
#include "TRandom.h"
#include "TSpline.h"
#include "TTree.h"

// Package Includes
#include "Detector.h"
#include "Spectra.h"
#include "Spectra1D.h"
#include "Spectra2D.h"
#include "Spectra3D.h"
#include "SpectraCorrDet.h"
#include "Utilities.h"

// Other External Includes
#include "dk2nu.h"
#include "dkmeta.h"
#include "calcLocationWeights.h"

namespace flxrd
{
  //---------------------------------------------------------------------------
  FluxReader::FluxReader(std::string fileWildcard,
                         unsigned int numFiles,
                         unsigned int skipFiles)
  {
    // Clear the input file vector
    if(fInputFiles.size() != 0) {
      fInputFiles.clear();
    }

    // Get the full list of files matching the input string
    fInputFiles = Wildcard(fileWildcard);

    // Remove first number of files specified by skipFiles
    if(skipFiles > fInputFiles.size()) {
      std::cout << "Warning: the number of files to skip is larger than the number of files found." << std::endl;
      std::cout << "There are currently no input files to run over." << std::endl;
      fInputFiles.clear();
    }
    else {
      fInputFiles.erase(fInputFiles.begin(), fInputFiles.begin() + skipFiles);
    }

    // Remove last number of files so that the number of remaining files is equal to numFiles
    if(numFiles > fInputFiles.size()) {
      std::cout << "numFiles input is larger than the number of files found." << std::endl;
      std::cout << "No files will be trimmed." << std::endl;
    }
    if(numFiles < fInputFiles.size() && numFiles != 0) {
      fInputFiles.erase(fInputFiles.begin() + numFiles, fInputFiles.end());
    }

    std::cout << fInputFiles.size() << " files were found matching the input criteria." << std::endl;

    if(fInputFiles.size() == 0) {
      std::cout << "Error: there are no files to run over. Asserting 0." << std::endl;
      assert(0);
    }

    fReweightNuRay = false; // By default, turn this off for speed

    fTreePath = "dk2nuTree";  // This is the default tree name in Dk2Nu files
    fMetaPath = "dkmetaTree"; // This is the default metadata tree name in Dk2Nu files
    fPOTPath  = "pots";       // This is the default POT variable name in Dk2Nu files
  }

  //---------------------------------------------------------------------------
  FluxReader::~FluxReader()
  {
  }

  //---------------------------------------------------------------------------
  void FluxReader::ReadFlux(TDirectory* out)
  {
    AddDefaultBranches(); // Add default branches to list of branches to turn on

    InitialMessage(); // Output the number and parameter types to be run over

    SetNuRayIndices(); // Setup the NuRay map so the detector name points to the first NuRay index for this detector

    // Make a TChain for all of the files
    TChain* fluxChain = new TChain(fTreePath.c_str());
    TChain* metaChain = new TChain(fMetaPath.c_str());
    for(unsigned int i_file = 0, n_file = fInputFiles.size(); i_file < n_file; ++i_file) {
      fluxChain->Add(fInputFiles[i_file].c_str());
      metaChain->Add(fInputFiles[i_file].c_str());
    }

    std::cout << "Looping over " << fluxChain->GetNtrees() << " trees." << std::endl;

    SetBranches(fluxChain, metaChain); // Turn on the necessary branches

    std::cout << "BEGIN!" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl << std::endl;

    int totEntries = 0;  // Total entries over all input files
    double totPOT  = 0.; // Sum of POT found in each file (an int is too small to store this number)
    int treeNumber = -1; // Store the tree number corresponding to the previous entry

    unsigned int i_entry = 0;
    while(metaChain->GetEntry(i_entry)) {
      ++i_entry;

      totPOT += fMeta->pots;
    }

    i_entry = 0; // Reset the entry number to 0
    while(fluxChain->GetEntry(i_entry)) {
      ++i_entry;

      // Let the user know where things stand periodically
      ++totEntries;
      if(totEntries % 250000 == 0) {
        std::cout << "On entry " << totEntries << "." << std::endl;
      }

      // Let the user know when moving to a new tree, i.e., a new file
      if(treeNumber != fluxChain->GetTreeNumber()) {
        treeNumber = fluxChain->GetTreeNumber();
        std::cout << "Moving to tree number " << treeNumber << "." << std::endl;
      }

      // Only the NuRay energy and weight change by detector,
      // so only execute this block if those variables are needed
      if(fReweightNuRay) {
        for(const Detector& det : fDetectors) {
          int index = fNuRayIndex[det.GetDetName()]; // Get the NuRay index for this detector

          double energy = 0., propwt = 0.;

          if(det.GetUses() == 1) {
            TVector3 xyz(0., 0., 0.);
            ToBeamCoords(det, xyz); // Convert coordinates to beam coordinates
            bsim::calcEnuWgt(fNu, xyz, energy, propwt); // Perform the reweight calculation
            fNu->nuray[index].E = energy; // Store the new energy
            fNu->nuray[index].wgt = propwt; // Store the new weight
          }
          else { // Same as above, but perform the calculation for each use at once
            for(int i_use = 0, n_use = det.GetUses(); i_use < n_use; ++i_use) {
              TVector3 xyz = Smear(det); // Smear the ray through the detector
              ToBeamCoords(det, xyz);
              bsim::calcEnuWgt(fNu, xyz, energy, propwt);
              fNu->nuray[index + i_use].E = energy;
              fNu->nuray[index + i_use].wgt = propwt;
            }
          } // end of conditionals if detector uses is 1
        } // end of loop over detectors
      } // end of conditional if NuRay needs to be reweighted

      // Fill histograms with values read from the entry
      for(unsigned int i_spec = 0, n_spec = fSpectra.size(); i_spec < n_spec; ++i_spec) {
        fSpectra[i_spec]->Fill(fNu, fNuRayIndex);
      }

    } // end of loop over flux tree entries

    TDirectory* temp = gDirectory; // Store the current directory to go back to this after running/writing is complete
    out->cd();

    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Total POT: " << totPOT << std::endl;
    std::cout << "Number of entries: " << totEntries << std::endl;

    // Create total POT histogram
    TH1D* hPOT = new TH1D("TotalPOT", ";;POT", 1, 0., 1.);
    hPOT->SetBinContent(1, totPOT);

    // Write histograms to output file
    gDirectory->WriteTObject(hPOT); // Start by recording POT information

    for(unsigned int i_spec = 0, n_spec = fSpectra.size(); i_spec < n_spec; ++i_spec) {
      out->mkdir(fSpectra[i_spec]->GetTitle().c_str()); // Create directory in output file
      out->cd(fSpectra[i_spec]->GetTitle().c_str()); // Go to the new directory
      fSpectra[i_spec]->WriteHists(gDirectory); // Have Spectra object write out its contents
    }

    temp->cd(); // Return to the original directory
    return;
  }

  //---------------------------------------------------------------------------
  void FluxReader::AddSpectra(Parameters params, std::string title,
                              std::string labelx, std::vector<double> binsx, const Var& varx,
                              const Weight& wei, TObject* extWeights)
  {
    // Create the new Spectra1D object
    Spectra1D* s = new Spectra1D(params, title,
                                 labelx, binsx, varx,
                                 wei, extWeights);
    fSpectra.push_back(s); // Add it to the vector of Spectra

    AddBranches(s->BranchesToAdd()); // Add necessary branches to master list

    return;
  }

  //---------------------------------------------------------------------------
  void FluxReader::AddSpectra(Parameters params, std::string title,
                              std::string labelx, std::vector<double> binsx, const Var& varx,
                              std::string labely, std::vector<double> binsy, const Var& vary,
                              const Weight& wei, TObject* extWeights)
  {
    // Create the new Spectra2D object
    Spectra2D* s = new Spectra2D(params, title,
                                 labelx, binsx, varx, 
                                 labely, binsy, vary,
                                 wei, extWeights);
    fSpectra.push_back(s);

    AddBranches(s->BranchesToAdd()); // Add necessary branches to master list

    return;
  }

  //---------------------------------------------------------------------------
  void FluxReader::AddSpectra(Parameters params, std::string title,
                              std::string labelx, std::vector<double> binsx, const Var& varx,
                              std::string labely, std::vector<double> binsy, const Var& vary,
                              std::string labelz, std::vector<double> binsz, const Var& varz,
                              const Weight& wei, TObject* extWeights)
  {
    // Create the new Spectra3D object
    Spectra3D* s = new Spectra3D(params, title,
                                 labelx, binsx, varx,
                                 labely, binsy, vary,
                                 labelz, binsz, varz,
                                 wei, extWeights);
    fSpectra.push_back(s);

    AddBranches(s->BranchesToAdd()); // Add necessary branches to master list

    return;
  }

  //---------------------------------------------------------------------------
  void FluxReader::AddSpectra(Parameters params, std::string title,
                              std::string detX, std::string detY,
                              std::string labelx, std::vector<double> binsx, const Var& varx,
                              const Weight& wei, TObject* extWeights)
  {
    // Create the new SpectraCorrDet object
    SpectraCorrDet* s = new SpectraCorrDet(params, title,
                                           detX, detY,
                                           labelx, binsx, varx,
                                           wei, extWeights);
    fSpectra.push_back(s);

    AddBranches(s->BranchesToAdd()); // Add necessary branches to master list

    return;
  }

  //---------------------------------------------------------------------------
  void FluxReader::OverrideTreeName(std::string treepath)
  {
    fTreePath = treepath; // This will be the new tree path
    return;
  }

  //---------------------------------------------------------------------------
  void FluxReader::OverridePOTPath(std::string metapath, std::string potpath)
  {
    fMetaPath = metapath; // This will be the new metadata tree path
    fPOTPath =  potpath;  // This will be the new path to the POT variable
    return;
  }

  //---------------------------------------------------------------------------
  void FluxReader::OverrideDefaultVarName(std::string oldname, std::string newname)
  {
    // Don't add a new branch if there is nothing to replace
    if(fBranchNames.find(oldname) == fBranchNames.end()) {
      std::cout << oldname << " is not a default branch." << std::endl;
      return;
    }

    fBranchOverrides[oldname] = newname; // Point the old branch name to the new name
    return;
  }

  //---------------------------------------------------------------------------
  void FluxReader::AddBranch(std::string branchName)
  {
    fBranchNames.insert(branchName);
    return;
  }

  //---------------------------------------------------------------------------
  void FluxReader::AddBranches(std::set<std::string> branchNames)
  {
    fBranchNames.insert(branchNames.begin(), branchNames.end());
    return;
  }

  //---------------------------------------------------------------------------
  void FluxReader::AddDefaultBranches()
  {
    // If a ray needs to be reweighted, the calculation needs all of these values
    if(fBranchNames.find("nuray.E")   != fBranchNames.end() ||
       fBranchNames.find("nuray.wgt") != fBranchNames.end()) {
      AddBranch("nuray");
      AddBranch("nuray.E");
      AddBranch("nuray.wgt");
      AddBranch("decay");
      AddBranch("decay.ntype");
      AddBranch("decay.vx");
      AddBranch("decay.vy");
      AddBranch("decay.vz");
      AddBranch("decay.pdpx");
      AddBranch("decay.pdpy");
      AddBranch("decay.pdpz");
      AddBranch("decay.ppdxdz");
      AddBranch("decay.ppdydz");
      AddBranch("decay.pppz");
      AddBranch("decay.ppenergy");
      AddBranch("decay.ptype");
      AddBranch("decay.muparpx");
      AddBranch("decay.muparpy");
      AddBranch("decay.muparpz");
      AddBranch("decay.mupare");
      AddBranch("decay.necm");

      fReweightNuRay = true;
    }

    return;
  }

  //---------------------------------------------------------------------------
  void FluxReader::InitialMessage()
  {
    const int num_per_line = 8;

    std::cout << "Looping over flux files." << std::endl;

    const int n_spec = fSpectra.size();
    std::cout << n_spec << " histogram types will be created:" << std::endl;
    for(int i_spec = 0; i_spec < n_spec; ++i_spec) {
      std::cout << fSpectra[i_spec]->GetTitle();
      if((i_spec < n_spec-1) && ((i_spec+1) % num_per_line != 0)) {
        std::cout << ", ";
      }
      else {
        std::cout << std::endl;
      }
    }

    std::cout << std::endl;
    return;
  }

  //---------------------------------------------------------------------------
  bool FluxReader::IsStandardDk2Nu()
  {
    // Standard tree/variable names
    std::string standardTree = "dk2nuTree";
    std::string standardMeta = "dkmetaTree";
    std::string standardPOT  = "pots";

    if(fBranchOverrides.size() != 0    ||
       fTreePath.compare(standardTree) ||
       fMetaPath.compare(standardMeta) ||
       fPOTPath .compare(standardPOT)) {
      return false;
    }

    return true;
  }

  //---------------------------------------------------------------------------
  void FluxReader::SetBranches(TTree* fluxTree, TTree* metaTree)
  {
    // Start with all branches off
    fluxTree->SetBranchStatus("*", 0);
    metaTree->SetBranchStatus("*", 0);

    // This is the default block for using a normal Dk2Nu file
    if(IsStandardDk2Nu()) {
      fBranches.reserve(fBranchNames.size());

      for(const std::string& branch : fBranchNames) {
        // std::cout << "Turning on branch " << branch << std::endl;
        fluxTree->SetBranchStatus(branch.c_str(), 1); // Turn on the branch

        // Add the actual TBranch to the list of TBranches,
        // and abort if the branch does not exist to avoid a seg fault
        // std::cout << "Adding branch " << branch << std::endl;
        fBranches.push_back(fluxTree->GetBranch(branch.c_str()));
        if(!fBranches.back()) {
          std::cerr << "Tree has no branch \"" << branch
                    << "\". Asserting 0." << std::endl;
          assert(0);
        }

        fluxTree->AddBranchToCache(fBranches.back()); // Theoretically this speed things up...
      } // Loop over branch names

      // Turn on and add the branch for POT
      metaTree->SetBranchStatus(fPOTPath.c_str(), 1);
      fBranches.push_back(metaTree->GetBranch(fPOTPath.c_str()));
      if(!fBranches.back()) {
        std::cerr << "Tree has no branch \"" << fPOTPath
                  << "\". Asserting 0." << std::endl;
        assert(0);
      }
      metaTree->AddBranchToCache(fBranches.back());

      // Make sure these pointers are not pointing to something else
      fNu = 0;
      fMeta = 0;

      // Point the Dk2Nu object in the tree to the class Dk2Nu object, fNu 
      std::string fullTreePath = "dk2nu";
      fluxTree->SetBranchAddress(fullTreePath.c_str(), &fNu);

      // Point the DkMeta object in the tree to the class DkMeta object, fMeta 
      std::string fullMetaPath = "dkmeta";
      metaTree->SetBranchAddress(fullMetaPath.c_str(), &fMeta);
    }
    else {
      // Make sure these pointers are not pointing to something else
      fNu = 0;
      fMeta = 0;

      // Create a map with default Dk2Nu branch names pointing to the actual values in the Dk2Nu object, fNu
      std::map<std::string, void*> m = OverrideAddresses(fNu);

      for(const std::string& branch : fBranchNames) {
        std::string branchPath = "";

        if(m.find(branch) == m.end()) { // Check if the branch is a default Dk2Nu branch
          branchPath = branch; // Use the default name
        }
        else { // Get the name of the branch in the (non-Dk2Nu) tree from the corresponding Dk2Nu branch name
          branchPath = fBranchOverrides[branch];
        }

        fluxTree->SetBranchStatus(branchPath.c_str(), 1); // Turn on the branch
        fluxTree->SetBranchAddress(branchPath.c_str(), m[branch]); // Point the branch into the class Dk2Nu object 
      }

      // Turn on and set the branch for the POTs
      metaTree->SetBranchStatus(fPOTPath.c_str(), 1);
      metaTree->SetBranchAddress(fPOTPath.c_str(), &fMeta->pots);
    }

    const int num_per_line = 8;
    std::cout << "The following branches are active:" << std::endl;
    unsigned int i_branch = 0;
    unsigned int n_branch = fBranchNames.size();
    for(const std::string& branch : fBranchNames) {
      std::cout << branch;
      if((i_branch < n_branch-1) && ((i_branch+1) % num_per_line != 0)) {
        std::cout << ", ";
      }
      else {
        std::cout << std::endl;
      }

      ++i_branch;
    }
    std::cout << std::endl;

    // Make sure the NuRay vector in the class Dk2Nu object is large enough
    // to have an entry/index for all the detectors (and each detector usage)
    if(fReweightNuRay) {
      if(fNu->nuray.size() < (unsigned int)fNuRayIndex["znull"]) {
        for(unsigned int i = fNu->nuray.size(), n = fNuRayIndex["znull"]; i < n; ++i) {
          fNu->nuray.push_back(bsim::NuRay());
        }
      }
    }

    return;
  }

  //---------------------------------------------------------------------------
  void FluxReader::SetNuRayIndices()
  {
    for(unsigned int i_spec = 0, n_spec = fSpectra.size(); i_spec < n_spec; ++i_spec) {
      std::set<Detector> sdets = fSpectra[i_spec]->Detectors(); // Get the detectors to be used in this Spectra
      fDetectors.insert(sdets.begin(), sdets.end()); // Add the Spectra's list to the full FluxReader list
    }

    int index = 0;
    for(const Detector& det : fDetectors) {
      // This will be the first index in the NuRay vector in the Dk2Nu object corresponding to the current detector
      fNuRayIndex[det.GetDetName()] = index;

      index += det.GetUses(); // Increment index by the number of times to use the current detector
    }
    fNuRayIndex["znull"] = index; // This will signal the last NuRay index

    return;
  }

  //---------------------------------------------------------------------------
  TVector3 FluxReader::Smear(const Detector& det, double rr)
  {
    // ISSUE: Add z offset? Fiducial volume cut?

    // Get detector size
    double xrange = det.GetHalfSizeX();
    double yrange = det.GetHalfSizeY();
    double zrange = det.GetHalfSizeZ();

    // Randomly choose point in detector
    double x = gRandom->Uniform(-1.*xrange, xrange);
    double y = gRandom->Uniform(-1.*yrange, yrange);
    double z = gRandom->Uniform(-1.*zrange, zrange);

    // If the detector is not square, make sure point lies inside radius sqrt(rr)
    while((rr < (x*x + y*y)) && (rr > 0.)) {
      x = gRandom->Uniform(-1.*xrange, xrange);
      y = gRandom->Uniform(-1.*yrange, yrange);
    }

    return TVector3(x, y, z);
  }

  //---------------------------------------------------------------------------
  void FluxReader::ToBeamCoords(const Detector& det, TVector3& xyz)
  {
    // Copy old coordinates for readability
    double oldx = xyz.X();
    double oldy = xyz.Y();
    double oldz = xyz.Z();

    // Get the detector coordinates
    double detx = det.GetCoordX();
    double dety = det.GetCoordY();
    double detz = det.GetCoordZ();

    // Update coordinates
    xyz.SetX(oldx + detx);
    xyz.SetY(dety + oldy*TMath::Cos(3.323155*TMath::DegToRad()) +
                    oldz*TMath::Sin(3.323155*TMath::DegToRad()));
    xyz.SetZ(detz + oldz*TMath::Cos(3.323155*TMath::DegToRad()) -
                    oldy*TMath::Sin(3.323155*TMath::DegToRad()));

    return;
  }
}
