#include "Combiner.h"

// C/C++ Includes
#include <cassert>
#include <iostream>

// Root Includes
#include "TAxis.h"
#include "TClass.h"
#include "TCollection.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"
#include "TKey.h"

// Package Includes
#include "Parameters.h"
#include "ParticleParam.h"

namespace flxrd
{
  //---------------------------------------------------------------------------
  Combiner::Combiner(std::string out)
  {
    // Store current directory to come back to later
    TDirectory* temp = gDirectory;

    fOut = new TFile(out.c_str(), "UPDATE"); // Open the input file
    assert(fOut->IsOpen()); // Break if the file is unopened

    // Loop through the file, and find all spectra
    TIter iterSpec(fOut->GetListOfKeys());
    TKey* keySpec;
    while((keySpec = (TKey*)iterSpec())) {
      // This eliminates the TotalPOT histogram from being added
      if( keySpec->ReadObj()->IsA()->InheritsFrom(TDirectory::Class()) ) {
        fSpectra.insert(keySpec->GetName());
      }
    }

    delete keySpec;

    // For each Spectra object, build a Parameters object matching what was used to create the Spectra
    // Keep track of SpectraCorrDet, since these don't get combined
    std::vector<std::string> corrDetSpec;
    for(const std::string& spec: fSpectra) {
      fOut->cd(spec.c_str()); // Go into the Spectra directory

      // List of each parameter (aside from histogram titles, which is stored by the class)
      std::set<std::string> dets;
      std::set<std::string> nuflavs;
      std::set<std::string> parents;
      std::set<std::string> xsecs;

      // Loop through the Spectra folder, and find all detectors
      TIter iterDet(gDirectory->GetListOfKeys());
      TKey* keyDet;
      while((keyDet = (TKey*)iterDet())) {
        // Make sure only directories are added
        if( keyDet->ReadObj()->IsA()->InheritsFrom(TDirectory::Class()) ) {
          dets.insert(keyDet->GetName());
        }
      }

      delete keyDet;

      if(dets.size() == 0) { // This is a SpectraCorrDet
        corrDetSpec.push_back(spec);
      }
      else { // This is a Spectra1D/2D/3D
        // Store, and go into the first detector directory
        std::string firstDet = *(dets.begin());
        gDirectory->cd(firstDet.c_str());

        // Loop through the detector folder, over all histograms
        TIter iterHist(gDirectory->GetListOfKeys());
        TKey* keyHist;
        while((keyHist = (TKey*)iterHist())) {
          std::string histTitle = keyHist->GetName(); // Get a histogram name

          // Remove the histogram title from the beginning, which has form "title_"
          histTitle.erase(histTitle.begin(), histTitle.begin() + spec.length() + 1);
          // Remove detector name from the end, which has the form "_detector"
          histTitle.erase(histTitle.end() - firstDet.length() - 1, histTitle.end());
          // The histogram title will have the remaining form flav_par_xsec

          std::string nuflav = histTitle.substr(0, histTitle.find('_')); // Get the neutrino flavor
          nuflavs.insert(nuflav); // Insert it into the neutrino flavors list
          // Remove it from the histgoram title, leaving par_xsec
          histTitle.erase(histTitle.begin(), histTitle.begin() + nuflav.length() + 1);

          // Repeart for parent name, leaving just xsec
          std::string parent = histTitle.substr(0, histTitle.find('_'));
          parents.insert(parent);
          histTitle.erase(histTitle.begin(), histTitle.begin() + parent.length() + 1);

          xsecs.insert(histTitle); // Add xsec to list
        }

        // Create and set up a Parameters object with the parameters used to make the Spectra
        Parameters* p = new Parameters();
        SetupParameters(p, dets, nuflavs, parents, xsecs);
        fParamsMap[spec] = *p; // Add the Parameters object to the map

        delete keyHist;
      }
    }

    // Get rid of all SpectraCorrDet from the fSpectra list
    for(const std::string& spec: corrDetSpec) {
      fSpectra.erase(spec);
    }

    InitialMessage(); // Output information to user

    temp->cd(); // Go back to original directory
  }

  //---------------------------------------------------------------------------
  Combiner::~Combiner()
  {
    if(fOut) {
      fOut->Close();
    }
  }

  //---------------------------------------------------------------------------
  void Combiner::CombineNuFlavs()
  {
    std::string rep_str = "allnu"; // This string will replace the neutrino flavor name

    // Do nothing if neutrino flavors have already been combined
    if(CombineAlreadyCalled(rep_str)) {
      std::cout << "Neutrino flavors have already been combined." << std::endl;
      return;
    }

    TDirectory* temp = gDirectory;

    // Loop through each Spectra
    for(const std::string& spec: fSpectra) {
      fOut->cd(spec.c_str());

      // Store number of each parameter in this Spectra
      const unsigned int n_flav = fParamsMap[spec].NFlav();
      const unsigned int n_par  = fParamsMap[spec].NPar();
      const unsigned int n_xsec = fParamsMap[spec].NXSec();
      const unsigned int n_det  = fParamsMap[spec].NDet();

      for(unsigned int i_det = 0; i_det < n_det; ++i_det) {
        gDirectory->cd(fParamsMap[spec].GetDetName(i_det).c_str());

        for(unsigned int i_xsec = 0; i_xsec < n_xsec; ++i_xsec) {
          for(unsigned int i_par  = 0; i_par  < n_par;  ++i_par) {
            // This corresponds to the way Parameters does indexing.
            // NuFlav index 0 is used by not including "+ i_flav" at the end.
            int index =   n_flav*n_par*n_xsec*i_det
                        + n_flav*n_par*i_xsec
                        + n_flav*i_par;

            // Find/create the name of a stored histogram and copy it into a new histogram for combining
            std::string hName = spec + "_" + fParamsMap[spec].NameTag(index);
            TH1* h = (TH1*)((TH1*)gDirectory->Get(hName.c_str()))->Clone();

            for(unsigned int i_flav = 1; i_flav < n_flav; ++i_flav) {
              ++index; // This corresponds to an increment of the NuFlav index
              hName = spec + "_" + fParamsMap[spec].NameTag(index); // Next histogram name

              h->Add((TH1*)gDirectory->Get(hName.c_str())); // Add it into the combined histogram
            }

            // Replace neutrino flavor name by the replacement string
            // The histogram name has format title_nuflav_par_xsec_det
            int firstPos = hName.find('_');
            int secndPos = hName.find('_', firstPos+1); // Search only after position specified in second argument
            hName.replace(firstPos+1, secndPos-firstPos-1, rep_str);
            h->SetName(hName.c_str()); // Give the combined histogram the correct name for writing to file

            gDirectory->WriteTObject(h); // Save the histogram
          } // Loop over parents
        } // Loop over cross sections

        fOut->cd(spec.c_str());
      } // Loop over detectors

      fOut->cd();
    } // Loop over Spectra

    temp->cd();
    return;
  }

  //---------------------------------------------------------------------------
  void Combiner::CombineParents()
  {
    // This function is similar to CombineNuFlav; see its comments for more details.

    std::string rep_str = "allpar"; // This string will replace the parent name

    if(CombineAlreadyCalled(rep_str)) {
      std::cout << "Parents have already been combined." << std::endl;
      return;
    }

    TDirectory* temp = gDirectory;

    for(const std::string& spec: fSpectra) {
      fOut->cd(spec.c_str());

      const unsigned int n_flav = fParamsMap[spec].NFlav();
      const unsigned int n_par  = fParamsMap[spec].NPar();
      const unsigned int n_xsec = fParamsMap[spec].NXSec();
      const unsigned int n_det  = fParamsMap[spec].NDet();

      for(unsigned int i_det = 0; i_det < n_det; ++i_det) {
        gDirectory->cd(fParamsMap[spec].GetDetName(i_det).c_str());

        for(unsigned int i_xsec = 0; i_xsec < n_xsec; ++i_xsec) {
          for(unsigned int i_flav = 0; i_flav < n_flav; ++i_flav) {
            // This corresponds to the way Parameters does indexing.
            // Parent index 0 is used by not including the line "n_flav*i_par".
            int index =   n_flav*n_par*n_xsec*i_det
                        + n_flav*n_par*i_xsec
                        + i_flav;

            std::string hName = spec + "_" + fParamsMap[spec].NameTag(index);
            TH1* h = (TH1*)((TH1*)gDirectory->Get(hName.c_str()))->Clone();

            for(unsigned int i_par  = 1; i_par  < n_par;  ++i_par) {
              index += n_par; // This corresponds to an increment of the Parent index
              hName = spec + "_" + fParamsMap[spec].NameTag(index);

              h->Add((TH1*)gDirectory->Get(hName.c_str()));
            }

            // Replace neutrino flavor name by the replacement string
            // The histogram name has format title_nuflav_par_xsec_det
            int firstPos = hName.find('_');
            firstPos = hName.find('_', firstPos+1); // Update to the next occurence of '_'
            int secndPos = hName.find('_', firstPos+1);
            hName.replace(firstPos+1, secndPos-firstPos-1, rep_str);
            h->SetName(hName.c_str());

            gDirectory->WriteTObject(h);
          } // Loop over flavors
        } // Loop over cross sections

        fOut->cd(spec.c_str());
      } // Loop over detectors

      fOut->cd();
    } // Loop over Spectra

    temp->cd();
    return;
  }

  //---------------------------------------------------------------------------
  void Combiner::CombineAll()
  {
    // This function is similar to CombineNuFlav; see its comments for more details.

    std::string nu_str  = "allnu"; // Replacement string for neutrino flavor name
    std::string par_str = "allpar"; // Replacement string for parent name
    std::string rep_str = nu_str + "_" + par_str; // Replacement string for both

    if(CombineAlreadyCalled(rep_str)) {
      std::cout << "All possible plots have already been combined." << std::endl;
      return;
    }

    CombineNuFlavs(); // Combine neutrino flavors
    CombineParents(); // Combine neutrino parents

    TDirectory* temp = gDirectory;

    for(const std::string& spec: fSpectra) {
      fOut->cd(spec.c_str());

      const unsigned int n_flav = fParamsMap[spec].NFlav();
      const unsigned int n_par  = fParamsMap[spec].NPar();
      const unsigned int n_xsec = fParamsMap[spec].NXSec();
      const unsigned int n_det  = fParamsMap[spec].NDet();

      for(unsigned int i_det = 0; i_det < n_det; ++i_det) {
        gDirectory->cd(fParamsMap[spec].GetDetName(i_det).c_str());

        for(unsigned int i_xsec = 0; i_xsec < n_xsec; ++i_xsec) {
          int index =   n_flav*n_par*n_xsec*i_det
                      + n_flav*n_par*i_xsec;

          std::string hName = spec + "_" + fParamsMap[spec].NameTag(index);

          // Create string for one of the combined Parent histogram names
          int firstPos = hName.find('_');
          firstPos = hName.find('_', firstPos+1);
          int secndPos = hName.find('_', firstPos+1);
          hName.replace(firstPos+1, secndPos-firstPos-1, par_str);

          // Get the combined parent histogram
          TH1* h = (TH1*)((TH1*)gDirectory->Get(hName.c_str()))->Clone();

          for(unsigned int i_flav = 1; i_flav < n_flav; ++i_flav) {
            ++index; // Increment the flavor index
            hName = spec + "_" + fParamsMap[spec].NameTag(index);

            // Make sure the combined parent histogram is the one pulled from the file
            firstPos = hName.find('_');
            firstPos = hName.find('_', firstPos+1);
            secndPos = hName.find('_', firstPos+1);
            hName.replace(firstPos+1, secndPos-firstPos-1, par_str);

            h->Add((TH1*)gDirectory->Get(hName.c_str()));
          }

          // The parent name is already replaced by pulling the combined parent histograms
          // Now replace the neutrino flavor name
          firstPos = hName.find('_');
          secndPos = hName.find('_', firstPos+1);
          hName.replace(firstPos+1, secndPos-firstPos-1, nu_str);
          h->SetName(hName.c_str());

          gDirectory->WriteTObject(h);
        } // Loop over cross sections

        fOut->cd(spec.c_str());
      } // Loop over detectors

      fOut->cd();
    } // Loop over Spectra

    temp->cd();
    return;
  }

  //---------------------------------------------------------------------------
  bool Combiner::CombineAlreadyCalled(std::string search)
  {
    TDirectory* temp = gDirectory;

    std::string firstSpec = *(fSpectra.begin());
    fOut->cd(firstSpec.c_str()); // Go to the first spectra directory
    gDirectory->cd(fParamsMap[firstSpec].GetDetName(0).c_str()); // Go to the first detector directory

    // Loop over histograms in the directory
    TIter iter(gDirectory->GetListOfKeys());
    TKey* key;
    while((key = (TKey*)iter())) {
      std::string histTitle = key->GetName(); // Get a histogram name

      // Search for the "search" string in the histogram name
      if(histTitle.find(search) != std::string::npos) {
        temp->cd();
        return true;
      }
    }

    temp->cd();
    return false;
  }

  //---------------------------------------------------------------------------
  void Combiner::InitialMessage()
  {
    const int num_per_line = 8;
    const int n_spec = fSpectra.size();

    std::cout << "Found " << n_spec << " Spectra:" << std::endl;

    int i_spec = 0;
    for(const std::string& spec: fSpectra) {
      std::cout << spec;
      if((i_spec < n_spec-1) && ((i_spec+1) % num_per_line != 0)) {
        std::cout << ", ";
      }
      else {
        std::cout << std::endl;
      }

      ++i_spec;
    }
    std::cout << std::endl;

    for(const std::string& spec: fSpectra) {
      std::cout << "In Spectra " << spec << ":" << std::endl;

      std::cout << "Found " << fParamsMap[spec].NDet()  << " detectors," << std::endl;
      std::cout << "Found " << fParamsMap[spec].NXSec() << " cross sections," << std::endl;
      std::cout << "Found " << fParamsMap[spec].NPar()  << " parents," << std::endl;
      std::cout << "Found " << fParamsMap[spec].NFlav() << " flavors." << std::endl;
    }

    std::cout << std::endl;

    return;
  }

  //---------------------------------------------------------------------------
  void Combiner::SetupParameters(Parameters* params,
                                 std::set<std::string>& dets,
                                 std::set<std::string>& nuflavs,
                                 std::set<std::string>& parents,
                                 std::set<std::string>& xsecs)
  {
    params->ClearAll(); // Get rid of all defaults
    params->ResetNuFlavs(); // Add in all neutrino flavors

    // Neutrino flavor removal
    for(unsigned int i_flav = 0, n_flav = params->NFlav(); i_flav < n_flav; ) {
      std::string nuflav = params->fNuFlav[i_flav].GetName(); // Current flavor name

      // This condition means the flavor was not found in the user input file
      if(nuflavs.find(nuflav) == nuflavs.end()) {
        params->RemoveNuFlav(nuflav); // Remove the flavor
        --n_flav; // Update number of flavors
      }
      else {
        ++i_flav; // Move to next flavor
      }
    }

    // Parent name additions
    int par_num = 0; // dummy value used for pdg
    for(std::set<std::string>::iterator i_par = parents.begin(); i_par != parents.end(); ++i_par) {
      Parent p(*i_par, par_num); // Set up a parent object with a dummy pdg
      params->AddParent(p); // Add it to the parameters
      ++par_num; // Move to next dummy value (necessary to avoid error of adding same pdg twice)
    }

    // Cross Section name additions
    for(std::set<std::string>::iterator i_xsec = xsecs.begin(); i_xsec != xsecs.end(); ++i_xsec) {
      params->AddXSec(*i_xsec); // Nothing special here, just add cross section
    }

    // Detector additions
    for(std::set<std::string>::iterator i_det = dets.begin(); i_det != dets.end(); ++i_det) {
      Detector d(*i_det, "", 0.,0.,0., 0.,0.,0., 0); // Set up a detector object with dummy size and position
      params->AddDetector(d); // Add it to the parameters
    }

    return;
  }
}
