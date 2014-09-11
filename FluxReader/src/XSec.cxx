#include "XSec.h"

// C/C++ Includes
#include <cassert>
#include <iostream>
#include <stdlib.h>

// Root Includes
#include "TAxis.h"
#include "TFile.h"

// Package Includes
#include "Utilities.h"

namespace flxrd
{
  //----------------------------------------------------------------------
  XSec::XSec()
  {
    SetXSecFileName(); // Initialize the file name
  }

  //----------------------------------------------------------------------
  XSec::~XSec()
  {
  }

  //----------------------------------------------------------------------
  void XSec::SetXSecFileName(std::string override)
  {
    if(!override.compare("You_really_should_not_override_if_possible")) { // This is the default
      // This environment variable points to the folder with the most current version of the Genie cross sections
      fXSecFileName = std::getenv("GENIEXSECPATH");

      fXSecFileName += "/xsec_graphs_*_*.root"; // General format of the cross section filename
      std::vector<std::string> files = Wildcard(fXSecFileName); // Get a list of all files matching the wildcard

      if(files.size() >= 1) {
        fXSecFileName = files[0]; // Set the string to the first file in the list

        if(files.size() > 1) { // Show user all files if more than one appears
          std::cout << "More than one file was found." << std::endl; 
          for(unsigned int i = 0, n = files.size(); i < n; ++i) {
            std::cout << files[i] << std::endl;
          }
        }
      }
      else {
        std::cout << "No file could be found matching " << fXSecFileName << std::endl;
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
      
      fXSecFileName = override;
    }

    return;
  }

  //----------------------------------------------------------------------
  void XSec::ListIntTypes()
  {
    for(unsigned int i = 0, n = kIntType.size(); i < n; ++i) {
      std::cout << kIntType[i] << std::endl;
    }

    return;
  }

  //----------------------------------------------------------------------
  TGraph* XSec::GetGraph(int pdg, std::string tar, std::string type)
  {
    TGraph* gRet = new TGraph();

    // Check for necessary recursion combinations
    // If target is CH2, need to add graphs from C and 2 H's
    if(tar == "CH2") {
      gRet = GetGraphMath( GetGraph(pdg,"C",type), 1,
                           GetGraph(pdg,"H",type), 2 );
    }
    else { // Base TGraph pulling scheme
      // Create the string pointing to the correct cross section directory
      SetXSecGenStr(pdg, tar, type);

      // Open the cross section file
      TFile* f = new TFile( (GetXSecFileName()).c_str() );

      // Pull the cross section
      gRet = (TGraph*)f->Get( (GetXSecGenStr()).c_str() );

      f->Close(); // Close the file
    }

    gRet->SetTitle(MakeXSecTitle(pdg, tar, type).c_str()); // Set the graph title
    return gRet;
  }

  //----------------------------------------------------------------------
  TGraph* XSec::GetGraphRatio(int pdg1, std::string tar1, std::string type1,
                              int pdg2, std::string tar2, std::string type2)
  {
    TGraph* g1 = GetGraph(pdg1, tar1, type1); // Pull numerator graph
    TGraph* g2 = GetGraph(pdg2, tar2, type2); // Pull denominator graph

    // Make sure both graphs have the same number of points, and store that number
    // We could (should?) also check that x(g1) == x(g2), but this is currently not done
    assert(g1->GetN() == g2->GetN());
    const int n = g1->GetN();

    // Arrays to store the TGraph points
    double x[n], y1[n], y2[n], yRet[n];

    for(int i = 0; i < n; ++i) {
      // Get each of the points
      g1->GetPoint(i, x[i], y1[i]);
      g2->GetPoint(i, x[i], y2[i]);

      if(y2[i] != 0) {
        yRet[i] = y1[i]/y2[i];
      }
      else {
        yRet[i] = 0;
      }

      if( !(yRet[i] > 0 && yRet[i] < 1000000000) ) { // Hack to make sure no NaN or Inf appear
        yRet[i] = 0;
      }
    }

    TGraph* gRet = new TGraph(n, x, yRet); // Create the new graph
    gRet->SetTitle(MakeXSecRatioTitle(pdg1, tar1, type1, pdg2, tar2, type2).c_str());

    return gRet;
  }

  //----------------------------------------------------------------------
  double XSec::XSecEval(TSpline3* s, double x)
  {
    double xsec = (s->Eval(x) > 0) ? s->Eval(x) : 0;
    return xsec;
  }

  //----------------------------------------------------------------------
  TSpline3* XSec::GetXSec(int pdg, std::string tar, std::string type,
                          const char* opt, double begin_val)
  {
    TGraph* g = GetGraph(pdg, tar, type); // Get the graph to generate the spline

    TSpline3* s = new TSpline3("", g, opt, begin_val); // Generate spline
    s->SetTitle(g->GetTitle());

    return s;
  }

  //----------------------------------------------------------------------
  TSpline3* XSec::GetXSecRatio(int pdg1, std::string tar1, std::string type1,
                               int pdg2, std::string tar2, std::string type2,
                               const char* opt, double begin_val)
  {
    TGraph* g = GetGraphRatio(pdg1, tar1, type1,
                              pdg2, tar2, type2); // Get the graph to generate the spline

    TSpline3* s = new TSpline3("", g, opt, begin_val); // Generate spline
    s->SetTitle(g->GetTitle());

    return s;
  }

  //----------------------------------------------------------------------
  TH1* XSec::GetHist(TSpline3* s, int nbins, double xmin, double xmax)
  {
    TH1* ret =  new TH1D("",";Energy (GeV);", nbins, xmin, xmax);

    TAxis* ax = ret->GetXaxis();

    double bin_val = 0; // This value will be assigned to a histogram bin 
    double delta_x = (xmax-xmin)/(10*nbins); // Width of (most) discrete integral trapezoids
    if(delta_x > 0.1) { // Max width of discrete integral trapezoids
      delta_x = 0.1;
    }

    for(int i = 1; i <= nbins; ++i) {
      // Calculate bin average from discrete integral over bin width
      // Calculate discrete integral of spline through current bin using trapezoids
      // Integral I(f(x), a, b) = delta_x * { [f(0) + f(N)]/2 + f(1) + f(2) + ... f(N-1) }, x(0) = a, x(N) = b
      if(delta_x < 0.1) { // Calculate discrete integral breaking up bin into 10 equal parts
        for(int i_x = 1; i_x < 10; i_x++) {
          bin_val += XSecEval(s, (ax->GetBinLowEdge(i) + i_x*delta_x) ); // Adding f(1) through f(N-1)
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

        bin_val += 0.5*( XSecEval(s, ax->GetBinLowEdge(i)) + XSecEval(s, max_loop_x) );
        bin_val *= delta_x;

        bin_val += 0.5*(ax->GetBinUpEdge(i) - max_loop_x)
                      *( XSecEval(s, max_loop_x) + XSecEval(s, ax->GetBinUpEdge(i)) );

        bin_val /= (ax->GetBinUpEdge(i) - ax->GetBinLowEdge(i));
      }

      ret->SetBinContent(i, bin_val); // Set bin content to histogram bin
      bin_val = 0; // Reset bin value calculation variable
    } // Loop over histogram bins

    ret->SetTitle(s->GetTitle());
    return ret;
  }

  //----------------------------------------------------------------------
  TH1* XSec::GetHist(TSpline3* s, int nbins, const double* edges)
  {
    // This function is similar to the equally spaced bins version above
    // See its comments for more details

    TH1* ret =  new TH1D("",";Energy (GeV);", nbins, edges);

    TAxis* ax = ret->GetXaxis();

    double bin_val = 0;
    double delta_x = 0.1; 

    for(int i = 1; i <= nbins; ++i) {
      delta_x = (ax->GetBinUpEdge(i) - ax->GetBinLowEdge(i))/10;

      if(delta_x < 0.1) {
        for(int i_x = 1; i_x < 10; i_x++) {
          bin_val += XSecEval(s, (ax->GetBinLowEdge(i) + i_x*delta_x) );
        }

        bin_val += 0.5*( XSecEval(s, ax->GetBinLowEdge(i)) + XSecEval(s, ax->GetBinUpEdge(i)) );
        bin_val *= delta_x;
        bin_val /= ax->GetBinUpEdge(i) - ax->GetBinLowEdge(i);
      }
      else {
        delta_x = 0.1;

        int num_eval_pts = (ax->GetBinUpEdge(i) - ax->GetBinLowEdge(i))*10;
        double max_loop_x   = ax->GetBinLowEdge(i) + num_eval_pts*delta_x;

        for(int i_x = 1; i_x < num_eval_pts; i_x++) {
          bin_val += XSecEval(s, (ax->GetBinLowEdge(i) + i_x*delta_x) );
        }

        bin_val += 0.5*( XSecEval(s, ax->GetBinLowEdge(i)) + XSecEval(s, max_loop_x) );
        bin_val *= delta_x;

        bin_val += 0.5*(ax->GetBinUpEdge(i) - max_loop_x)
                      *( XSecEval(s, max_loop_x) + XSecEval(s, ax->GetBinUpEdge(i)) );

        bin_val /= ax->GetBinUpEdge(i) - ax->GetBinLowEdge(i);
      }

      ret->SetBinContent(i, bin_val);
      bin_val = 0;
    } // Loop over histogram bins

    ret->SetTitle(s->GetTitle());
    return ret;
  }

  //----------------------------------------------------------------------
  TGraph* XSec::GetGraphMath(TGraph* g1, double c1, TGraph* g2, double c2)
  {
    assert(g1->GetN() == g2->GetN());
    const int n = g1->GetN();
    double x[n], y1[n], y2[n], yRet[n];

    for(int i=0; i<n; ++i) {
      g1->GetPoint(i, x[i], y1[i]);
      g2->GetPoint(i, x[i], y2[i]);

      yRet[i] = c1*y1[i] + c2*y2[i];
    }

    TGraph* gRet = new TGraph(n, x, yRet);
    return gRet;
  }

  //----------------------------------------------------------------------
  std::string XSec::MakeXSecTitle(int pdg, std::string tar, std::string type)
  {
    // All titles have the form:
    // "Cross Section: #nu_{<flavor>} <current> Scattering from <nucleus>" 

    std::string title = "Cross Section: ";

    if(pdg < 0) { 
      title += "Anti-";
    }

    if(abs(pdg) == 12) {
      title += "#nu_{e} ";
    }
    if(abs(pdg) == 14) {
      title += "#nu_{#mu} ";
    }
    if(abs(pdg) == 16) {
      title += "#nu_{#tau} ";
    }

    if(     type == "tot_cc") {
      title += "CC";
    }
    else if(type == "tot_nc") {
      title += "NC";
    }
    else {
      title += type;
    }

    title += " Scattering from ";
    if(tar == "H")   {
      title += "^{1}H";
    }
    if(tar == "C")   {
      title += "^{12}C";
    }
    if(tar == "N")   {
      title += "^{14}N";
    }
    if(tar == "O")   {
      title += "^{16}O";
    }
    if(tar == "S")   {
      title += "^{32}S";
    }
    if(tar == "Cl")  {
      title += "^{35}Cl";
    }
    if(tar == "Ti")  {
      title += "^{48}Ti";
    }
    if(tar == "Fe")  {
      title += "^{56}Fe";
    }
    if(tar == "CH2") {
      title += "CH_{2}";
    }

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
      if(pdg1 < 0) {
        title += "Anti-";
      }
      if(abs(pdg1) == 12) {
        title += "#nu_{e} to ";
      }
      if(abs(pdg1) == 14) {
        title += "#nu_{#mu} to ";
      }
      if(abs(pdg1) == 16) {
        title += "#nu_{#tau} to ";
      }

      if(pdg2 < 0) {
        title += "Anti-";
      }
      if(abs(pdg2) == 12) {
        title += "#nu_{e} ";
      }
      if(abs(pdg2) == 14) {
        title += "#nu_{#mu} ";
      }
      if(abs(pdg2) == 16) {
        title += "#nu_{#tau} ";
      }

      if(     type1 == "tot_cc") {
        title += "CC";
      }
      else if(type1 == "tot_nc") {
        title += "NC";
      }
      else                       {
        title += type1;
      }

      title += " Scattering from ";
      if(tar1 == "H")   {
        title += "^{1}H";
      }
      if(tar1 == "C")   {
        title += "^{12}C";
      }
      if(tar1 == "N")   {
        title += "^{14}N";
      }
      if(tar1 == "O")   {
        title += "^{16}O";
      }
      if(tar1 == "S")   {
        title += "^{32}S";
      }
      if(tar1 == "Cl")  {
        title += "^{35}Cl";
      }
      if(tar1 == "Ti")  {
        title += "^{48}Ti";
      }
      if(tar1 == "Fe")  {
        title += "^{56}Fe";
      }
      if(tar1 == "CH2") {
        title += "CH_{2}";
      }
    }
    else if(pdg1 == pdg2) {
      if(pdg1 < 0) {
        title += "Anti-";
      }
      if(abs(pdg1) == 12) {
        title += "#nu_{e} ";
      }
      if(abs(pdg1) == 14) {
        title += "#nu_{#mu} ";
      }
      if(abs(pdg1) == 16) {
        title += "#nu_{#tau} ";
      }

      if(     type1 == "tot_cc") {
        title += "CC";
      }
      else if(type1 == "tot_nc") {
        title += "NC";
      }
      else                       {
        title += type1;
      }

      title += " Scattering from ";
      if(tar1 == "H")   {
        title += "^{1}H";
      }
      if(tar1 == "C")   {
        title += "^{12}C";
      }
      if(tar1 == "N")   {
        title += "^{14}N";
      }
      if(tar1 == "O")   {
        title += "^{16}O";
      }
      if(tar1 == "S")   {
        title += "^{32}S";
      }
      if(tar1 == "Cl")  {
        title += "^{35}Cl";
      }
      if(tar1 == "Ti")  {
        title += "^{48}Ti";
      }
      if(tar1 == "Fe")  {
        title += "^{56}Fe";
      }
      if(tar1 == "CH2") {
        title += "CH_{2}";
      }
      title += " to ";

      if(     type2 == "tot_cc") {
        title += "CC";
      }
      else if(type2 == "tot_nc") {
        title += "NC";
      }
      else                       {
        title += type2;
      }

      title += " Scattering from ";
      if(tar2 == "H")   {
        title += "^{1}H";
      }
      if(tar2 == "C")   {
        title += "^{12}C";
      }
      if(tar2 == "N")   {
        title += "^{14}N";
      }
      if(tar2 == "O")   {
        title += "^{16}O";
      }
      if(tar2 == "S")   {
        title += "^{32}S";
      }
      if(tar2 == "Cl")  {
        title += "^{35}Cl";
      }
      if(tar2 == "Ti")  {
        title += "^{48}Ti";
      }
      if(tar2 == "Fe")  {
        title += "^{56}Fe";
      }
      if(tar2 == "CH2") {
        title += "CH_{2}";
      }
    }
    else {
      if(pdg1 < 0) {
        title += "Anti-";
      }
      if(abs(pdg1) == 12) {
        title += "#nu_{e} ";
      }
      if(abs(pdg1) == 14) {
        title += "#nu_{#mu} ";
      }
      if(abs(pdg1) == 16) {
        title += "#nu_{#tau} ";
      }

      if(     type1 == "tot_cc") {
        title += "CC";
      }
      else if(type1 == "tot_nc") {
        title += "NC";
      }
      else                       {
        title += type1;
      }

      title += " Scattering from ";
      if(tar1 == "H")   {
        title += "^{1}H";
      }
      if(tar1 == "C")   {
        title += "^{12}C";
      }
      if(tar1 == "N")   {
        title += "^{14}N";
      }
      if(tar1 == "O")   {
        title += "^{16}O";
      }
      if(tar1 == "S")   {
        title += "^{32}S";
      }
      if(tar1 == "Cl")  {
        title += "^{35}Cl";
      }
      if(tar1 == "Ti")  {
        title += "^{48}Ti";
      }
      if(tar1 == "Fe")  {
        title += "^{56}Fe";
      }
      if(tar1 == "CH2") {
        title += "CH_{2}";
      }
      title += " to ";

      if(pdg2 < 0) {
        title += "Anti-";
      }
      if(abs(pdg2) == 12) {
        title += "#nu_{e} ";
      }
      if(abs(pdg2) == 14) {
        title += "#nu_{#mu} ";
      }
      if(abs(pdg2) == 16) {
        title += "#nu_{#tau} ";
      }

      if(     type2 == "tot_cc") {
        title += "CC";
      }
      else if(type2 == "tot_nc") {
        title += "NC";
      }
      else                       {
        title += type2;
      }

      title += " Scattering from ";
      if(tar2 == "H")   {
        title += "^{1}H";
      }
      if(tar2 == "C")   {
        title += "^{12}C";
      }
      if(tar2 == "N")   {
        title += "^{14}N";
      }
      if(tar2 == "O")   {
        title += "^{16}O";
      }
      if(tar2 == "S")   {
        title += "^{32}S";
      }
      if(tar2 == "Cl")  {
        title += "^{35}Cl";
      }
      if(tar2 == "Ti")  {
        title += "^{48}Ti";
      }
      if(tar2 == "Fe")  {
        title += "^{56}Fe";
      }
      if(tar2 == "CH2") {
        title += "CH_{2}";
      }
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
  void XSec::SetXSecGenStr(int pdg, std::string tar, std::string type)
  {
    // String has form "pdg_(bar_)Ab##/InteractionType",
    // where Ab is the atom code and ## is the number of nucleons
    bool valid_input = false;
    fXSecGenStr = "nu_"; // All directories in the cross section file begin with this string

    for(unsigned int i = 0, n = kNuPDG.size(); i < n; ++i) {
      if(pdg == kNuPDG[i]) { // Next component in the directory string is the neutrino type
        if(abs(pdg) == 12) {
          fXSecGenStr += "e_";
        }
        if(abs(pdg) == 14) {
          fXSecGenStr += "mu_";
        }
        if(abs(pdg) == 16) {
          fXSecGenStr += "tau_";
        }

        valid_input = true;
      }
    }
    if(!valid_input) { // Check for valid neutrino pdg. If invalid, show user list before quitting
      std::cout << "Invalid neutrino. Please input one of the following:" << std::endl;
      for(unsigned int i = 0, n = kNuPDG.size(); i < n; ++i) {
        std::cout << kNuPDG[i] << std::endl;
      }
    }
    assert(valid_input);
    valid_input = false; // Reset valid_input for next piece of label

    if(pdg < 0) {
      fXSecGenStr += "bar_"; // Add if necessary to directory string
    }

    for(unsigned int i = 0, n = kTarget.size(); i < n; ++i) {
      if(tar == kTarget[i]) { // Last component of directory string is target.
        if(tar == "H")  { 
          fXSecGenStr += "H1/";
        }
        if(tar == "C")  { 
          fXSecGenStr += "C12/";
        }
        if(tar == "N")  { 
          fXSecGenStr += "N14/";
        }
        if(tar == "O")  { 
          fXSecGenStr += "O16/";
        }
        if(tar == "S")  { 
          fXSecGenStr += "S32/";
        }
        if(tar == "Cl") {
          fXSecGenStr += "Cl35/";
        }
        if(tar == "Ti") {
          fXSecGenStr += "Ti48/";
        }
        if(tar == "Fe") {
          fXSecGenStr += "Fe56/";
        }

        valid_input = true;
      }
    }
    if(!valid_input) { // Check for valid target. If invalid, show user list before quitting
      std::cout << "Invalid target. Please input one of the following:" << std::endl;
      for(unsigned int i = 0, n = kTarget.size(); i < n; ++i) {
        std::cout << kTarget[i] << std::endl;
      }
    }
    assert(valid_input);
    valid_input = false;

    for(unsigned int i = 0, n = kIntType.size(); i < n; ++i) {
      if(type == kIntType[i]) { // After directory, add the interaction type
        fXSecGenStr += type;
        valid_input = true;
      }
    }
    if(!valid_input) { // Check for valid interaction type. If invalid, give user advice before quitting
      std::cout << "Invalid interaction type. Please input one of the following:" << std::endl;
      std::cout << "For the most general CC or NC, input tot_cc or tot_nc, respectively." << std::endl;
      std::cout << "For a full list of all the inputs, call the function ListIntTypes()." << std::endl;
    }
    assert(valid_input);

    return;
  }
}
