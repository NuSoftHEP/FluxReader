#pragma once

// C/C++ Includes
#include <vector>

// Root Includes
#include "TGraph.h"
#include "TH1.h"
#include "TSpline.h"

namespace flxrd
{
  /// Interface to cross-section information
  class XSec
  {
    public: 
      XSec();
      ~XSec();

      /// Helper function which sets the string for the file which contains the cross section information
      /// This should normally be left to its default value, and typically is not a function called by the user
      void        SetXSecFileName(std::string override = "You_really_should_not_override_if_possible");

      /// Get the file name containing the cross section information
      std::string GetXSecFileName() const { return fXSecFileName; }

      /// Get the string pointing to a directory containing specific cross section information
      std::string GetXSecGenStr() const { return fXSecGenStr; }

      /// List all the cross section interaction types; for example, "tot_cc".
      void ListIntTypes();

      /// Pulls the cross section for a given target, neutrino, interaction combinaton as a TGraph
      TGraph* GetGraph(int pdg, std::string tar, std::string type);

      /// Generates a cross section ratio as a TGraph
      /// Inputs with "1" are numerator; "2" are denominator
      TGraph* GetGraphRatio(int pdg1, std::string tar1, std::string type1,
                            int pdg2, std::string tar2, std::string type2);

      /// Generates a cross section plot as a TSpline3*
      TSpline3* GetXSec(int pdg, std::string tar, std::string type,
                        const char* opt = "", double begin_val = 0);

      /// Generates a cross section ratio as a TSpline3*
      /// Inputs with "1" are numerator; "2" are denominator
      TSpline3* GetXSecRatio(int pdg1, std::string tar1, std::string type1,
                             int pdg2, std::string tar2, std::string type2,
                             const char* opt = "", double begin_val = 0);

      /// Generates a cross section plot as a TH1* with equally spaced bins
      TH1* GetHist(TSpline3* s, int nbins, double min, double max);
      TH1* GetHist(TSpline3* s, int nbins, const double* edges);

      /// Evaluates s(x), returning s(x) if s(x) > 0 and 0 otherwise, as cross sections are always positive
      double XSecEval(TSpline3* s, double x);

      const std::vector<int>         kNuPDG   = {12,-12,14,-14,16-16}; ///< Valid neutrino PDG inputs

      const std::vector<std::string> kTarget  = {"H","C","N","O","S","Cl","Ti","Fe","CH2"}; ///< Valid target strings

      const std::vector<std::string> kIntType = {"qel_nc_p",
        "res_cc_p_1232P33","res_cc_p_1620S31","res_cc_p_1700D33","res_cc_p_1910P31","res_cc_p_1920P33","res_cc_p_1905F35",
        "res_cc_p_1950F37","res_nc_p_1232P33","res_nc_p_1535S11","res_nc_p_1520D13","res_nc_p_1650S11","res_nc_p_1700D13",
        "res_nc_p_1675D15","res_nc_p_1620S31","res_nc_p_1700D33","res_nc_p_1440P11","res_nc_p_1720P13","res_nc_p_1680F15",
        "res_nc_p_1910P31","res_nc_p_1920P33","res_nc_p_1905F35","res_nc_p_1950F37","res_nc_p_1710P11",
        "dis_cc_p_ubarsea","dis_cc_p_dval","dis_cc_p_dsea","dis_cc_p_ssea","dis_nc_p_sbarsea",
        "dis_nc_p_ubarsea","dis_nc_p_dbarsea","dis_nc_p_dval","dis_nc_p_dsea","dis_nc_p_uval",
        "dis_nc_p_usea","dis_nc_p_ssea","dis_cc_p_dval_charm","dis_cc_p_dsea_charm","dis_cc_p_ssea_charm",
        "qel_cc_p_charm4222","imd_cc","ve_nc","res_cc_p","res_cc_n","res_nc_p","res_nc_n",
        "dis_cc_p","dis_cc_n","dis_nc_p","dis_nc_n","dis_cc_p_charm","dis_cc_n_charm","dis_nc_p_charm","dis_nc_n_charm",
        "mec_cc","mec_nc","tot_cc","tot_cc_p","tot_cc_n","tot_nc","tot_nc_p","tot_nc_n"}; ///< Valid interaction types

    protected:
      /// Helper function which merges two TGraphs via c1*g1 + c2*g2.
      TGraph* GetGraphMath(TGraph* g1, double c1, TGraph* g2, double c2);

      /// Helper functions which generate a string to be used as a histogram title
      std::string MakeXSecTitle(     int pdg,  std::string tar,  std::string type);
      std::string MakeXSecRatioTitle(int pdg1, std::string tar1, std::string type1,
                                     int pdg2, std::string tar2, std::string type2);

      /// Helper function which generates the string which will match one of the directories in the cross section file
      void        SetXSecGenStr(int pdg, std::string tar, std::string type = "tot_cc");

      std::string fXSecFileName; ///< Filename containing the cross section information to pull
      std::string fXSecGenStr; ///< Label pointing to the correct cross section in the file pointed to by fXSecFileName
  };
}
