#include "XSec.h"

// C/C++ Includes
#include <cassert>
#include <iostream>
#include <stdlib.h>

// Root Includes
#include "TAxis.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TKey.h"
#include "TSpline.h"

// Package Includes
#include "Utilities.h"

namespace flxrd
{
  //----------------------------------------------------------------------
  XSec::XSec()
  {
    fLocalDir = gDirectory;
    fXSecFile = new TFile();

    SetXSecFile(); // Initialize the file

    // Source: http://www.chemeddl.org/resources/ptl/index.php
    fMolarMass["H"]  = 1.008;
    fMolarMass["C"]  = 12.011;
    fMolarMass["N"]  = 14.007;
    fMolarMass["O"]  = 15.999;
    fMolarMass["S"]  = 32.065;
    fMolarMass["Cl"] = 35.453;
    fMolarMass["Ar"] = 39.948;
    fMolarMass["Ti"] = 47.867;
    fMolarMass["Fe"] = 55.845;
  }

  //----------------------------------------------------------------------
  XSec::~XSec()
  {
    if(fXSecFile->IsOpen()) {
      fXSecFile->Close();
    }
  }

  //----------------------------------------------------------------------
  TGraph* XSec::GetGraph(int pdg, std::string tar, std::string type, bool eventRate)
  {
    TGraph* gRet = new TGraph();

    // Copy this string into a local variable that can be passed by reference
    std::string typeCopy(type);

    // Check for necessary recursion combinations
    // If target is a compound, need to add relevant graphs
    if(fTarget.find(tar) == fTarget.end()) {
      gRet = GetGraphCompound(tar, pdg, type);
    }
    else { // Base TGraph pulling scheme
      // Create the string pointing to the correct cross section directory
      SetXSecGenStr(pdg, tar, typeCopy);

      // Pull the cross section
      gRet = (TGraph*)fXSecFile->Get(GetXSecGenStr().c_str());
    }

    // If this cross section is for an event rate, scale the y values appropriately
    if(eventRate) {
      // Avogadro's Number * 10^-38 cm^2 * 10^9 g/kton / Molar Mass
      double scale = 0.0000060221413/fMolarMass[tar]; // Units: cm^2/kton

      // Scale each y value
      double x = 0., y = 0.;
      for(int i = 0, n = gRet->GetN(); i < n; ++i) {
        gRet->GetPoint(i, x, y);
        gRet->SetPoint(i, x, y*scale);
      }
    }

    // Set the graph title and axis labels
    std::string histTitle = MakeXSecTitle(pdg, tar, typeCopy) + ";Energy (GeV);";
    histTitle += (eventRate ? "cm^{2}/kton" : "10^{-38} cm^{2}"); // Add the appropriate y label
    gRet->SetTitle(histTitle.c_str());

    return gRet;
  }

  //----------------------------------------------------------------------
  TGraph* XSec::GetGraphRatio(int pdg1, std::string tar1, std::string type1,
                              int pdg2, std::string tar2, std::string type2,
                              bool eventRate)
  {
    TGraph* g1 = GetGraph(pdg1, tar1, type1); // Pull numerator graph with no scaling
    TGraph* g2 = GetGraph(pdg2, tar2, type2); // Pull denominator graph with no scaling

    // Make sure both graphs have the same number of points, and store that number
    // We could (should?) also check that x(g1) == x(g2), but this is currently not done
    assert(g1->GetN() == g2->GetN());
    const int n = g1->GetN();

    // Arrays to store the TGraph points
    double x[n], y1[n], y2[n], yRet[n];

    // For ratios, the scale for event rates is just the inverse ratio of the two molar masses
    double eventRateScale = 1.;
    if(eventRate) {
      eventRateScale = fMolarMass[tar2]/fMolarMass[tar1];
    }

    for(int i = 0; i < n; ++i) {
      // Get each of the points
      g1->GetPoint(i, x[i], y1[i]);
      g2->GetPoint(i, x[i], y2[i]);

      if(y2[i] != 0) {
        yRet[i] = eventRateScale*y1[i]/y2[i];
      }
      else {
        yRet[i] = 0.;
      }

      if( !(yRet[i] > 0. && yRet[i] < 1000000000.) ) { // Hack to make sure no NaN or Inf appear
        yRet[i] = 0.;
      }
    }

    TGraph* gRet = new TGraph(n, x, yRet); // Create the new graph

    // The interaction type passed to MakeXSecRatioTitle must be modified separately
    // But don't print the warning message a second time (SetXSecGenStr already will have)
    std::string typeCopy1(type1);
    NuElectronCheck(typeCopy1, pdg1, false);
    std::string typeCopy2(type2);
    NuElectronCheck(typeCopy2, pdg2, false);

    std::string histTitle = MakeXSecRatioTitle(pdg1, tar1, typeCopy1, pdg2, tar2, typeCopy2);
    histTitle += ";Energy (GeV);";
    gRet->SetTitle(histTitle.c_str());

    return gRet;
  }

  //----------------------------------------------------------------------
  TSpline3* XSec::GetXSec(int pdg, std::string tar, std::string type,
                          bool eventRate, const char* opt, double begin_val)
  {
    TGraph* g = GetGraph(pdg, tar, type, eventRate); // Get the graph to generate the spline

    TSpline3* s = new TSpline3("", g, opt, begin_val); // Generate spline

    s->SetTitle(g->GetTitle()); // This should pick up the actual title and axes labels

    return s;
  }

  //----------------------------------------------------------------------
  TSpline3* XSec::GetXSecRatio(int pdg1, std::string tar1, std::string type1,
                               int pdg2, std::string tar2, std::string type2,
                               bool eventRate, const char* opt, double begin_val)
  {
    TGraph* g = GetGraphRatio(pdg1, tar1, type1,
                              pdg2, tar2, type2,
                              eventRate); // Get the graph to generate the spline

    TSpline3* s = new TSpline3("", g, opt, begin_val); // Generate spline

    s->SetTitle(g->GetTitle()); // This should pick up the actual title and axes labels

    return s;
  }

  //----------------------------------------------------------------------
  TH1* XSec::GetHist(TSpline3* s, int nbins, double xmin, double xmax)
  {
    TH1* ret = new TH1D("", "", nbins, xmin, xmax);

    TAxis* ax = ret->GetXaxis();

    double bin_val = 0.; // This value will be assigned to a histogram bin 
    double delta_x = (xmax - xmin)/(10.*(double)nbins); // Width of (most) discrete integral trapezoids
    if(delta_x > 0.1) { // Max width of discrete integral trapezoids
      delta_x = 0.1;
    }

    for(int i = 1; i <= nbins; ++i) {
      // Calculate bin average from discrete integral over bin width
      // Calculate discrete integral of spline through current bin using trapezoids
      // Integral I(f(x), a, b) = delta_x * { [f(0) + f(N)]/2 + f(1) + f(2) + ... f(N-1) }, x(0) = a, x(N) = b
      if(delta_x < 0.1) { // Calculate discrete integral breaking up bin into 10 equal parts
        for(int i_x = 1; i_x < 10; i_x++) {
          bin_val += XSecEval(s, ax->GetBinLowEdge(i) + i_x*delta_x); // Adding f(1) through f(N-1)
        }

        bin_val += 0.5*( XSecEval(s, ax->GetBinLowEdge(i)) + XSecEval(s, ax->GetBinUpEdge(i)) ); // Add [f(0) and f(N)]/2
        bin_val *= delta_x; // Multiply by delta_x
        bin_val /= (ax->GetBinUpEdge(i) - ax->GetBinLowEdge(i)); // Compute average
      }
      else { // Calculate using trapezoids of width 0.1 GeV as is possible, use smaller width for last trapezoid
        int num_eval_pts = (ax->GetBinUpEdge(i) - ax->GetBinLowEdge(i))*10;
        double max_loop_x = ax->GetBinLowEdge(i) + num_eval_pts*delta_x;

        for(int i_x = 1; i_x < num_eval_pts; i_x++) {
          bin_val += s->Eval(ax->GetBinLowEdge(i) + i_x*delta_x);
        }

        bin_val += 0.5*(XSecEval(s, ax->GetBinLowEdge(i)) + XSecEval(s, max_loop_x));
        bin_val *= delta_x;

        bin_val += 0.5*(ax->GetBinUpEdge(i) - max_loop_x)
                      *(XSecEval(s, max_loop_x) + XSecEval(s, ax->GetBinUpEdge(i)));

        bin_val /= (ax->GetBinUpEdge(i) - ax->GetBinLowEdge(i));
      }

      ret->SetBinContent(i, bin_val); // Set bin content to histogram bin
      bin_val = 0; // Reset bin value calculation variable
    } // Loop over histogram bins

    ret->SetTitle(s->GetTitle()); // This should pick up the actual title and y axis label (if applicable)
    return ret;
  }

  //----------------------------------------------------------------------
  TH1* XSec::GetHist(TSpline3* s, int nbins, const double* edges)
  {
    // This function is similar to the equally spaced bins version above
    // See its comments for more details

    TH1* ret = new TH1D("", "", nbins, edges);

    TAxis* ax = ret->GetXaxis();

    double bin_val = 0.;
    double delta_x = 0.1;

    for(int i = 1; i <= nbins; ++i) {
      delta_x = (ax->GetBinUpEdge(i) - ax->GetBinLowEdge(i))/10.;

      if(delta_x < 0.1) {
        for(int i_x = 1; i_x < 10; i_x++) {
          bin_val += XSecEval(s, (ax->GetBinLowEdge(i) + i_x*delta_x) );
        }

        bin_val += 0.5*(XSecEval(s, ax->GetBinLowEdge(i)) + XSecEval(s, ax->GetBinUpEdge(i)));
        bin_val *= delta_x;
        bin_val /= ax->GetBinUpEdge(i) - ax->GetBinLowEdge(i);
      }
      else {
        delta_x = 0.1;

        int num_eval_pts = (ax->GetBinUpEdge(i) - ax->GetBinLowEdge(i))*10;
        double max_loop_x = ax->GetBinLowEdge(i) + num_eval_pts*delta_x;

        for(int i_x = 1; i_x < num_eval_pts; i_x++) {
          bin_val += XSecEval(s, ax->GetBinLowEdge(i) + i_x*delta_x);
        }

        bin_val += 0.5*(XSecEval(s, ax->GetBinLowEdge(i)) + XSecEval(s, max_loop_x));
        bin_val *= delta_x;

        bin_val += 0.5*(ax->GetBinUpEdge(i) - max_loop_x)
                      *(XSecEval(s, max_loop_x) + XSecEval(s, ax->GetBinUpEdge(i)));

        bin_val /= ax->GetBinUpEdge(i) - ax->GetBinLowEdge(i);
      }

      ret->SetBinContent(i, bin_val);
      bin_val = 0.;
    } // Loop over histogram bins

    ret->SetTitle(s->GetTitle()); // This should pick up the actual title and y axis label (if applicable)
    return ret;
  }

  //----------------------------------------------------------------------
  bool XSec::IsValidProcess(std::string type) const
  {
    for(const auto& allowedType : fIntType) {
      if(!type.compare(allowedType)) { return true; }
    }

    return false;
  }

  //----------------------------------------------------------------------
  void XSec::ListBaseTargets() const
  {
    for(const auto& allowedTar : fTarget) {
      std::cout << allowedTar.first << std::endl;
    }

    return;
  }

  //----------------------------------------------------------------------
  void XSec::ListIntTypes() const
  {
    for(const auto& allowedType : fIntType) {
      std::cout << allowedType << std::endl;
    }

    return;
  }

  //----------------------------------------------------------------------
  void XSec::ListMolarMasses() const
  {
    for(const auto& molarMass : fMolarMass) {
      std::cout << molarMass.first << ": " << molarMass.second << std::endl;
    }

    return;
  }

  //----------------------------------------------------------------------
  void XSec::ListNuPDGs() const
  {
    for(const auto& allowedNu : fNuPDG) {
      std::cout << allowedNu.first << std::endl;
    }

    return;
  }

  //----------------------------------------------------------------------
  void XSec::SetXSecFile(std::string override)
  {
    if(fXSecFile->IsOpen()) {
      fXSecFile->Close();
    }

    std::string xsecFileName;

    if(!override.compare("You_really_should_not_override_if_possible")) { // This is the default
      // This environment variable points to the folder with the most current version of the Genie cross sections
      xsecFileName = std::getenv("GENIEXSECPATH");

      xsecFileName += "/xsec_graphs_*_*.root"; // General format of the cross section filename
      std::vector<std::string> files = Wildcard(xsecFileName); // Get a list of all files matching the wildcard

      if(files.size() >= 1) {
        xsecFileName = files[0]; // Set the string to the first file in the list

        if(files.size() > 1) { // Show user all files if more than one appears
          std::cout << "More than one file was found." << std::endl; 
          for(const auto& fileName : files) {
            std::cout << fileName << std::endl;
          }
        }
      }
      else {
        std::cout << "No file could be found matching " << xsecFileName << std::endl;
      }

      assert(files.size() >= 1); // Make sure a file was found
    }
    else { // Warn the user about what he/she is doing.
      std::cout << "Warning: you are attempting to override the default file name." << std::endl;
      std::cout << "Please make sure you have a good reason for doing so!" << std::endl;
      std::cout << "The only error check that will be performed on your input" << std::endl;
      std::cout << "is whether it has a .root extension. Otherwise, you're on your own." << std::endl;

      std::string rootextension = override; // Check that the override string has .root at its end
      rootextension.erase(0, rootextension.size()-5);
      
      if(rootextension != ".root") { // Abort if .root is not at the end of the string
        std::cout << std::endl;
        std::cout << "This override file does not have a .root extension." << std::endl;
      }
      assert(rootextension == ".root");
      
      xsecFileName = override;
    }

    fXSecFile = new TFile(xsecFileName.c_str(), "READ");
    if(fXSecFile->IsOpen()) {
      SetupValidInputs();
    }
    else {
      std::cout << "Error: could not open file " << xsecFileName << "." << std::endl;
    }

    fLocalDir->cd();
    return;
  }

  //----------------------------------------------------------------------
  double XSec::XSecEval(TSpline3* s, double x)
  {
    double xsec = (s->Eval(x) > 0.) ? s->Eval(x) : 0.;
    return xsec;
  }

  //----------------------------------------------------------------------
  void XSec::AddElementToCompoundMap(std::string tar, int number,
                                     std::map<std::string, int>* map)
  {
    // First check that the base target is valid
    if(fTarget.find(tar) == fTarget.end()) {
      std::cout << "Warning: " << tar << " is not a valid target." << std::endl
                << "The compound must contain all valid targets." << std::endl
                << "Doing nothing with this target." << std::endl;
      return;
    }
    else { // The current target is valid...
      // Check if the target already exists in the map
      // Add the coefficient to the previous one if so
      if(map->find(tar) != map->end()) {
        map->at(tar) += number;
        std::cout << "Warning: " << tar << " was found in the compound multiple times."
                  << std::endl << "The total number of " << tar << " atoms is now "
                  << map->at(tar) << "." << std::endl;
      }
      else {
        // Create a new entry
        map->insert(std::pair<std::string, int>(tar, number));
      }
    }

    return;
  }

  //----------------------------------------------------------------------
  TGraph* XSec::GetGraphCompound(std::string compound, int pdg, std::string type)
  {
    // Character sets for string tokenization
    std::string upperCase = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    std::string lowerCase = "abcdefghijklmnopqrstuvwxyz";
    std::string number    = "0123456789";

    // Make sure the first character is correct so the following while loop works
    int i_compound = compound.find_first_of(upperCase);
    if(i_compound != 0) {
      std::cout << "Enter a valid compound. Each atom should have correct capitalization." << std::endl
                << "This means that the first letter should be capitalized..." << std::endl;
      abort();
    }

    // Map to contain the number of each atom
    std::map<std::string, int> targetsByNumber;

    // Loop over the compound capital letters
    while(i_compound != std::string::npos) {
      if(i_compound + 1 == compound.size() ||
         compound.find_first_of(upperCase, i_compound + 1) == i_compound + 1) {
        // Single letter atomic symbol, one atom

        std::string tar = compound.substr(i_compound, 1); // This is the atomic symbol

        AddElementToCompoundMap(tar, 1, &targetsByNumber); // Add this to the map
      } // End of conditional that the atomic symbol has one letter and one atom
      else if(compound.find_first_of(number, i_compound + 1) == i_compound + 1) {
        // Single letter atomic symbol, multiple atoms

        std::string tar = compound.substr(i_compound, 1);

        // The position in the compound string of the first character for the coefficient
        int numberBegin = compound.find_first_of(number, i_compound + 1);
        // The number of characters in the coefficient
        int numberLength = compound.find_first_not_of(number, numberBegin + 1) - numberBegin;
        // Create a substring that is just the coefficient; convert it to an int
        int number = std::atoi(compound.substr(numberBegin, numberLength).c_str());

        AddElementToCompoundMap(tar, number, &targetsByNumber);
      }
      else if(compound.find_first_of(lowerCase, i_compound + 1) == i_compound + 1) {
        // Double letter atomic symbol

        // If the symbol has more than 1 letter, it only has 2
        std::string tar = compound.substr(i_compound, 2);

        if(i_compound + 2 == compound.size() ||
           compound.find_first_of(upperCase, i_compound + 2) == i_compound + 2) {
          // Double letter atomic symbol, one atom
          AddElementToCompoundMap(tar, 1, &targetsByNumber);
        }
        else if(compound.find_first_of(number, i_compound + 2) == i_compound + 2) {
          // Double letter atomic symbol, multiple atoms
          int numberBegin = compound.find_first_of(number, i_compound + 1);
          int numberLength = compound.find_first_not_of(number, numberBegin) - numberBegin;
          int number = std::atoi(compound.substr(numberBegin, numberLength).c_str());

          AddElementToCompoundMap(tar, number, &targetsByNumber);
        }
        else {
          // Double letter atomic symbol followed by neither a capital letter nor number
          std::cout << "Warning: Unknown character following element " << tar << ". " << std::endl
                    << "Compounds should only contain letters and numbers," << std::endl
                    << "and each new element should begin with a capital letter." << std::endl
                    << "Skipping this element." << std::endl;
        }
      }
      else {
        // Single letter atomic symbol followed by neither a letter nor number
        std::string tar = compound.substr(i_compound, 1); // Get the atomic symbol
        std::cout << "Warning: Unknown character following element " << tar << ". " << std::endl
                  << "Compounds should only contain letters and numbers." << std::endl
                  << "Skipping this element." << std::endl;
      }

      // Move to the next atomic symbol
      i_compound = compound.find_first_of(upperCase, i_compound + 1);
    } // end of while loop over atomic symbols (capital letters)

    if(targetsByNumber.size() == 0) {
      std::cout << "Error: the compound contained no valid base atoms." << std::endl;
      abort();
    }

    // Get the number of points in each graph
    const int n = GetGraph(pdg, targetsByNumber.begin()->first, type)->GetN();
    double x[n]; // Array of x values
    double y = 0.; // Temp variables for the y value of each graph point
    double yTot[n]; // The final y values array

    double molarMass = 0.; // Molar mass of the compound

    // Initialize arrays
    for(int i = 0; i < n; ++i) {
      x[i] = 0.;
      yTot[i] = 0.;
    }

    // Loop over each target/number of atoms pair
    for(const auto& tarPair : targetsByNumber) {
      std::string tar = tarPair.first;

      // Get the relevant base graph
      TGraph* g = GetGraph(pdg, tar, type);

      // Cast the coefficient of the current target as a double
      double coeff = (double)tarPair.second;

      // Loop over each point in the current graph
      for(int i = 0; i < n; ++i) {
        // Get the point, add it to the running total with the proper coefficient
        g->GetPoint(i, x[i], y);
        yTot[i] += coeff*y;
      } // end of loop over points in the current graph

      molarMass += coeff*fMolarMass[tar];
    } // end of for loop over each target pair

    // Otherwise, add the molar mass to the map
    fMolarMass[compound] = molarMass;

    // Create the graph to be returned
    TGraph* gRet = new TGraph(n, x, yTot);
    return gRet;
  }

  //----------------------------------------------------------------------
  std::string XSec::MakeXSecTitle(int pdg, std::string tar, std::string type)
  {
    // All titles have the form:
    // "Cross Section: #nu_{<flavor>} <current> Scattering from <nucleus>" 
    std::string title = "Cross Section: " +
                        TitleFlavor(pdg) +
                        TitleProcess(type) +
                        " Scattering from " +
                        TitleTarget(tar);

    return title;
  }

  //----------------------------------------------------------------------
  std::string XSec::MakeXSecRatioTitle(int pdg1, std::string tar1, std::string type1,
                                       int pdg2, std::string tar2, std::string type2)
  {
    // All titles have the form:
    // "Cross Section Ratio: #nu_{<flavor>} <current> Scattering from <nucleus> to
    //                       #nu_{<flavor>} <current> Scattering from <nucleus>",
    // unless enough parameters are similar that this can be simplified

    std::string title = "Cross Section Ratio: ";

    // First check if the title can be simplified a bit.
    if((tar1 == tar2) && (type1 == type2)) {
      title += TitleFlavor(pdg1) + "to " + TitleFlavor(pdg2) +
               TitleProcess(type1) +
               " Scattering from " +
               TitleTarget(tar1);
    }
    else if(pdg1 == pdg2) {
      title += TitleFlavor(pdg1) +
               TitleProcess(type1) +
               " Scattering from " +
               TitleTarget(tar1) +
               " to " +
               TitleProcess(type2) +
               " Scattering from " +
               TitleTarget(tar2);
    }
    else {
      title += TitleFlavor(pdg1) +
               TitleProcess(type1) +
               " Scattering from " +
               TitleTarget(tar1) +
               " to " +
               TitleFlavor(pdg2) +
               TitleProcess(type2) +
               " Scattering from " +
               TitleTarget(tar2);
    }

    // Easter egg.
    if((pdg1 == pdg2) && (tar1 == tar2) && (type1 == type2)) {
      if(pdg1 == 12 || pdg1 == -12) {
        title = "http://www.youtube.com/watch?v=WM8bTdBs-cw";
      }
      if(pdg1 == 14 || pdg1 == -14) {
        title = "http://www.youtube.com/watch?v=ftjEcrrf7r0";
      }
      if(pdg1 == 16 || pdg1 == -16) {
        title = "http://www.youtube.com/watch?v=UiKcd7yPLdU";
      }
    }

    return title;
  }

  //----------------------------------------------------------------------
  void XSec::NuElectronCheck(std::string& type, int pdg, bool inSetXSecGenStr)
  {
    if     (!type.compare("ve_nc") && abs(pdg) == 12) {
      if(inSetXSecGenStr) {
        std::cout << type << " is not available for electron or anti-electron neutrinos." << std::endl;
        std::cout << "Changing type to ve_ccncmix." << std::endl;
      }
      type = "ve_ccncmix";
    }
    else if(!type.compare("ve_ccncmix") && abs(pdg) != 12) {
      if(inSetXSecGenStr) {
        std::cout << type << " is only available for electron and anti-electron neutrinos." << std::endl;
        std::cout << "Changing type to ve_nc." << std::endl;
      }
      type = "ve_nc";
    }

    return;
  }

  //----------------------------------------------------------------------
  void XSec::SetupValidInputs()
  {
    fIntType.clear();
    fNuPDG.clear();
    fTarget.clear();

    TDirectory* tmp = gDirectory;
    fXSecFile->cd();

    // Loop over all keys in the top directory
    TIter nextkeyTop(gDirectory->GetListOfKeys());
    TKey* key;
    TKey* oldkey = 0;

    while((key = (TKey*)nextkeyTop())) {
      // Keep only the highest cycle number for each key
      if(oldkey && !strcmp(oldkey->GetName(), key->GetName())) { continue; }

      std::string name(key->GetName()); // Get the name of the object (should be a folder)

      // Find the length of the flavor name, which has the form nu_<flav>_(bar_)
      int endOfNuName = ((name.find("bar_") != std::string::npos) ?
                         name.find("bar_") + 4 :
                         name.find('_', 3) + 1);

      // Copy the flavor name into a local variable, erase it from the full name
      std::string nuName = name.substr(0, endOfNuName);
      name.erase(0, nuName.size());

      // Calculate the PDG of the neutrino
      int pdg;
      if     (nuName.find("_e_")   != std::string::npos) { pdg = 12; }
      else if(nuName.find("_mu_")  != std::string::npos) { pdg = 14; }
      else if(nuName.find("_tau_") != std::string::npos) { pdg = 16; }
      else { std::cout << "Error: Unknown neutrino flavor." << std::endl; }

      if(nuName.find("bar") != std::string::npos) { pdg *= -1; }

      // Add the PDG and name to the valid flavor map
      fNuPDG.insert(std::pair<int, std::string>(pdg, nuName));

      // What is left should be Ab##, but Ab may only be one character
      // Get just the code in a local variable, but keep the full string
      std::string atomicCode = name.substr(0, name.find_first_of("0123456789"));

      // Add the target and its full name to the valid target map
      std::pair<std::string, std::string> targetPair(atomicCode, name);
      fTarget.insert(targetPair);
    } // end of while loop over TKeys

    // Next, get the interaction process types
    // One loop over a subdirectory should suffice,
    // with a small caveat discussed after the loop
    std::string subdir = fNuPDG.begin()->second + fTarget.begin()->second;
    fXSecFile->cd(subdir.c_str());

    TIter nextkeySub(gDirectory->GetListOfKeys());
    oldkey = 0;

    while((key = (TKey*)nextkeySub())) {
      // Keep only the highest cycle number for each key
      if(oldkey && !strcmp(oldkey->GetName(), key->GetName())) { continue; }

      std::string name(key->GetName()); // Get the name of the object (should be a TGraph)

      // Add the process to the list of valid processes
      fIntType.insert(name);
    } // end of while loop over TKeys

    // Scattering off of electrons differs based on the neutrino flavor
    // If one of the standard processes was found, add the other
    if     (fIntType.find("ve_ccncmix") != fIntType.end()) {
      fIntType.insert("ve_nc");
    }
    else if(fIntType.find("ve_nc") != fIntType.end()) {
      fIntType.insert("ve_ccncmix");
    }

    tmp->cd();

    return;
  }

  //----------------------------------------------------------------------
  void XSec::SetXSecGenStr(int pdg, std::string tar, std::string& type)
  {
    // String has form "pdg_(bar_)Ab##/InteractionType",
    // where Ab is the atom code and ## is the number of nucleons

    if(fNuPDG.find(pdg) != fNuPDG.end()) {
      fXSecGenStr = fNuPDG.find(pdg)->second;
    }
    else { // If invalid neutrino, show user list before quitting
      std::cout << "Invalid neutrino. Please input one of the following:" << std::endl;
      for(const auto& allowedPDG : fNuPDG) {
        std::cout << allowedPDG.first << std::endl;
      }

      abort();
    }

    if(fTarget.find(tar) != fTarget.end()) {
      fXSecGenStr += fTarget.find(tar)->second + "/";
    }
    else { // If invalid target, show user list before quitting
      std::cout << "Invalid target. Please input one of the following:" << std::endl;
      for(const auto& allowedTar : fTarget) {
        std::cout << allowedTar.first << std::endl;
      }

      abort();
    }

    // After directory, add the interaction type
    if(fIntType.find(type) != fIntType.end()) {
      // If the requested cross section is scattering off of electrons,
      // make sure the the appropriate process is used based on the neutrino flavor
      NuElectronCheck(type, pdg, true);

      fXSecGenStr += type;
    }
    else { // Check for valid interaction type. If invalid, give user advice before quitting
      std::cout << "Invalid interaction type." << std::endl;
      std::cout << "For the most general CC or NC, input tot_cc or tot_nc, respectively." << std::endl;
      std::cout << "For a full list of all the inputs, call the function ListIntTypes()." << std::endl;

      abort();
    }

    return;
  }

  //----------------------------------------------------------------------
  std::string XSec::TitleFlavor(const int& pdg) const
  {
    // The format for the neutrino flavor is (Anti-)#nu_{(#)flav}

    // Start with "#nu_flav_(bar_)"
    std::string titleHelper = "#" + fNuPDG.find(pdg)->second;

    // Find the first letter of the flavor
    int flavorBegin  = titleHelper.find("_") + 1;
    // Calculate the number of characters of the flavor
    int flavorLength = titleHelper.find("_", flavorBegin) - flavorBegin;

    // Create a helper string, start with just "{"
    std::string flavor = "{";
    // If the flavor is a Greek symbol, add a "#", to get "{#"
    if(titleHelper.find("nu_e_") == std::string::npos) { flavor += "#"; }
    // Add the flavor and final "}" to get "{(#)flav}"
    flavor += titleHelper.substr(flavorBegin, flavorLength) + "}";

    // Replace the flavor in the original string with the helper
    // "#nu_flav_(bar_)" becomes "#nu_{(#)flav}(bar_)
    titleHelper.replace(flavorBegin, flavorLength + 1, flavor);

    // Replace "bar_" at the end by "Anti-" at the beginning, if applicable
    if(titleHelper.find("bar_") != std::string::npos) {
      titleHelper.erase(titleHelper.find("bar_"), 4);
      titleHelper = "Anti-" + titleHelper;
    }

    titleHelper += " ";

    return titleHelper;
  }

  //----------------------------------------------------------------------
  std::string XSec::TitleProcess(const std::string& type) const
  {
    // The format for the process is exactly the process, unless it is CC or NC
    if(!type.compare("tot_cc")) { return "CC"; }
    if(!type.compare("tot_nc")) { return "NC"; }

    return type;
  }

  //----------------------------------------------------------------------
  std::string XSec::TitleTarget(const std::string& tar) const
  {
    // Character set of numbers
    std::string number = "0123456789";

    std::string titleHelper = tar;

    if(fTarget.find(tar) != fTarget.end()) { // Base atom
      // Start with "Ab##"
      titleHelper = fTarget.find(tar)->second;

      // Find where the coefficient begins
      int coeffBegin  = titleHelper.find_first_of(number);
      // Get the number of characters in the coefficient
      int coeffLength = titleHelper.find_last_of(number) - coeffBegin + 1;

      // Get the coefficient as a separate string, erase it from the helper
      std::string coeff = titleHelper.substr(coeffBegin, coeffLength);
      titleHelper.erase(coeffBegin, coeffLength);

      // Prepend the coefficient as a superscript to the left of the atomic symbol
      titleHelper = "^{" + coeff + "}" + titleHelper;
    }
    else { // Compound
      int coeffBegin = titleHelper.find_first_of(number); // Find the first coefficient

      while(coeffBegin != std::string::npos) {
        // Find the first non numeric character after the current one
        int coeffEnd = titleHelper.find_first_not_of(number, coeffBegin);
        // Calculate the number of characters in the coefficient
        int coeffLength = ((coeffEnd != std::string::npos) ?
                           coeffEnd - coeffBegin :
                           titleHelper.size() - coeffBegin);

        // Replace the coefficient with a subscript version
        std::string rep_str = "_{" + titleHelper.substr(coeffBegin, coeffLength) + "}";
        titleHelper.replace(coeffBegin, coeffLength, rep_str);

        // Move to the next coefficient, taking into account the additional string characters
        coeffBegin = titleHelper.find_first_of(number, 2 + coeffBegin + coeffLength);
      } // end of loop over compound coefficients
    }

    return titleHelper;
  }
}
