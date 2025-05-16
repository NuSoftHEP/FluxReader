#pragma once

// C/C++ Includes
#include <map>
#include <set>
#include <string>
#include <vector>

// Forward class declarations
class TDirectory;
class TFile;
class TGraph;
class TSpline3;
class TH1;

namespace flxrd
{
  /// Interface to cross-section information
  class XSec
  {
    public: 
      XSec();
      ~XSec();

      /// Pulls the cross section for a given target, neutrino, interaction combinaton as a TGraph
      TGraph* GetGraph(int pdg, std::string tar, std::string type, bool eventRate = false);

      /// Generates a cross section ratio as a TGraph
      /// Inputs with "1" are numerator; "2" are denominator
      TGraph* GetGraphRatio(int pdg1, std::string tar1, std::string type1,
                            int pdg2, std::string tar2, std::string type2,
                            bool eventRate = false);

      /// Generates a cross section plot as a TSpline3*
      TSpline3* GetXSec(int pdg, std::string tar, std::string type, bool eventRate = false,
                        const char* opt = "", double begin_val = 0);

      /// Generates a cross section ratio as a TSpline3*
      /// Inputs with "1" are numerator; "2" are denominator
      TSpline3* GetXSecRatio(int pdg1, std::string tar1, std::string type1,
                             int pdg2, std::string tar2, std::string type2,
                             bool eventRate = false, const char* opt = "", double begin_val = 0);

      /// Generates a cross section plot as a TH1* with equally spaced bins
      TH1* GetHist(TSpline3* s, int nbins, double min, double max);
      TH1* GetHist(TSpline3* s, int nbins, const double* edges);

      /// Get the string pointing to a directory containing specific cross section information
      std::string GetXSecGenStr() const { return fXSecGenStr; }

      /// Check whether the input process is a valid one found in fIntType
      bool IsValidProcess(std::string type) const;

      /// List all the base atomic targets
      void ListBaseTargets() const;

      /// List all the cross section interaction types; for example, "tot_cc".
      void ListIntTypes() const;

      /// List all the molar masses that are stored
      void ListMolarMasses() const;

      /// List all of the PDG values of neutrinos
      void ListNuPDGs() const;

      /// Helper function which opens the file that contains the cross section information
      /// This should normally be left to its default value, and typically is not a function called by the user
      void SetXSecFile(std::string override = "You_really_should_not_override_if_possible");

      /// Evaluates s(x), returning s(x) if s(x) > 0 and 0 otherwise, as cross sections are always positive
      double XSecEval(TSpline3* s, double x);

    private:
      /// Helper function for computing the cross section using a compound
      /// When splitting a compound string into its constituent atoms,
      /// this function adds the element and coefficient to a map
      /// The map links the coefficient to the base atom
      void AddElementToCompoundMap(std::string tar, int number,
                                   std::map<std::string, int>* map);

      /// Generates a cross section from a chemical compound, like CH2
      TGraph* GetGraphCompound(std::string compound, int pdg, std::string type);

      /// Helper functions which generate a string to be used as a histogram title
      std::string MakeXSecTitle(     int pdg,  std::string tar,  std::string type);
      std::string MakeXSecRatioTitle(int pdg1, std::string tar1, std::string type1,
                                     int pdg2, std::string tar2, std::string type2);

      /// Make adjustment to the cross section process for scattering off electrons if necessary
      void NuElectronCheck(std::string& type, int pdg, bool inSetXSecGenStr);

      /// Read the input file and set up the valid user inputs
      void SetupValidInputs();

      /// Helper function which generates the string which will match one of the directories in the cross section file
      void        SetXSecGenStr(int pdg, std::string tar, std::string& type);

      /// Helper functions that convert a string to a format appropriate for a histogram title
      std::string TitleFlavor (const int& pdg) const;
      std::string TitleProcess(const std::string& type) const;
      std::string TitleTarget (const std::string& target) const;

      std::set<std::string> fIntType; ///< Valid interaction types

      TDirectory* fLocalDir; ///< Directory before opening fXSecFile used for the histogram scope

      std::map<std::string, double> fMolarMass; ///< Targets (not necessarily valid) and their molar masses 

      std::map<int, std::string> fNuPDG; ///< Valid neutrino PDG inputs and associated title strings

      std::map<std::string, std::string> fTarget; ///< Valid targets and their molar masses

      TFile* fXSecFile; ///< File containing cross section information
      std::string fXSecGenStr; ///< Label pointing to the correct cross section in the cross section file
  };
}
