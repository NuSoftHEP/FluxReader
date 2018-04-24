#include "Utilities.h"

// C/C++ Includes
#include <cstring>
#include <fstream>
#include <iostream>
#include "sys/stat.h"
#include "wordexp.h"

// Other External Includes
#include "dk2nu.h"

namespace flxrd
{
  //---------------------------------------------------------------------------
  std::vector<double> Bins(int nbins, double min, double max)
  {
    std::vector<double> bins; // Vector of bin edges

    // Calculate the width of each bin
    const double bin_width = (max - min)/((double)nbins);

    // Fill the vector with each bin edge
    for(int i = 0; i <= nbins; ++i) {
      bins.push_back(min + i*bin_width);
    }

    return bins;
  }

  //---------------------------------------------------------------------------
  std::vector<double> LoadDetCoords(std::string detName)
  {
    std::vector<double> coords; // Return vector

    // Set the path to the detector locations file
    std::string locationFilePath = std::getenv("DK2NU");
    locationFilePath += "/etc/locations.txt";

    // Open the detector locations file
    std::ifstream locationFile;
    locationFile.open(locationFilePath.c_str());

    // Make sure the file opened properly
    if(!locationFile.good()) {
      std::cout << "Error opening location file." << std::endl;
      return coords;
    }

    // Variables for file reading
    const int maxChar = 256;    // Max number of characters to read from a single line from the input file
    char locationLine[maxChar]; // Array to hold a single line from the input file
    bool detFound = false;      // Flag whether the requested input detector has been found

    // Loop over lines from the input file until the detector is found or the end of the file is reached
    while(!(detFound || locationFile.eof())) {
      // Get the next line
      if(locationFile.getline(locationLine, maxChar).fail()) {
        std::cout << "Error accessing the next line. " << std::endl
                  << "It may just be the end of the file, "
                  << "and the EOF bit may not have been set as expected." << std::endl;
        break;
      }

      // Make a std::string from the char array
      // This must be done BEFORE tokenizing the character array,
      // as strtok "destorys" the array as it goes
      std::string stringLine(locationLine);

      // Get the first set of contiguous characters that skips leading and contains no spaces
      // If the first non-space character is a '#', the line is a comment and can be skipped
      char* token = strtok(locationLine, " ");
      if(token[0] == '#') { continue; }

      if(stringLine.find(detName) != std::string::npos) {
        // If the input detector name is in the current line, set the boolean flag to terminate the while loop
        detFound = true;

        // Remove the detector from the line
        // All that should remain is leading/tailing white space,
        // and three coordinates separated by (any number of) spaces
        // Then, put the result back into the character array
        stringLine.erase(stringLine.find(detName), detName.size());
        sprintf(locationLine, "%s", stringLine.c_str());
      } // end of condition that the detector name was in the current line from the location file
    } // end of loop until detector name is found or end of the location file is reached

    // If the detector was not found, warn the user and exit
    if(!detFound) {
      std::cout << "Could not find " << detName << "." << std::endl;
      return coords;
    }

    // Tokenize the remainder of the line
    // Get the first coordinate
    char* coordChar = strtok(locationLine, " ");
    for(int i = 0, n = 3; i < n; ++i) {
      // Warn the user if there were not 3 coordinates found
      if(!coordChar) {
        std::cout << "Could only find " << i << " position entries in the line." << std::endl;
        return coords;
      }

      // Take the current "token," convert it to a double, add it to the return vector
      double coord = atof(coordChar);
      coords.push_back(coord);

      // Get the next coordinate (or token)
      // This should probably return a null pointer on the third loop, but that is okay
      coordChar = strtok(NULL, " ");
    }

    return coords;
  }

  //---------------------------------------------------------------------------
  std::map<std::string, void*> OverrideAddresses(bsim::Dk2Nu* nu)
  {
    std::map<std::string, void*> ret; // Map to return

    // Make sure the Dk2Nu object has at least one index in each of its vectors
    if(!nu->nuray.size())
      nu->nuray.push_back(bsim::NuRay());
    if(!nu->ancestor.size())
      nu->ancestor.push_back(bsim::Ancestor());
    if(!nu->traj.size())
      nu->traj.push_back(bsim::Traj());

    // Set nuray branches (nuray is a vector in Dk2Nu)
    ret["nuray.px"]  = &nu->nuray[0].px;
    ret["nuray.py"]  = &nu->nuray[0].py;
    ret["nuray.pz"]  = &nu->nuray[0].pz;
    ret["nuray.E"]   = &nu->nuray[0].E;
    ret["nuray.wgt"] = &nu->nuray[0].wgt;

    // Set decay branches
    ret["decay.norig"]    = &nu->decay.norig;
    ret["decay.ndecay"]   = &nu->decay.ndecay;
    ret["decay.ntype"]    = &nu->decay.ntype;
    ret["decay.vx"]       = &nu->decay.vx;
    ret["decay.vy"]       = &nu->decay.vy;
    ret["decay.vz"]       = &nu->decay.vz;
    ret["decay.pdpx"]     = &nu->decay.pdpx;
    ret["decay.pdpy"]     = &nu->decay.pdpy;
    ret["decay.pdpz"]     = &nu->decay.pdpz;
    ret["decay.ppdxdz"]   = &nu->decay.ppdxdz;
    ret["decay.ppdydz"]   = &nu->decay.ppdydz;
    ret["decay.pppz"]     = &nu->decay.pppz;
    ret["decay.ppenergy"] = &nu->decay.ppenergy;
    ret["decay.ppmedium"] = &nu->decay.ppmedium;
    ret["decay.ptype"]    = &nu->decay.ptype;
    ret["decay.muparpx"]  = &nu->decay.muparpx;
    ret["decay.muparpy"]  = &nu->decay.muparpy;
    ret["decay.muparpz"]  = &nu->decay.muparpz;
    ret["decay.mupare"]   = &nu->decay.mupare;
    ret["decay.necm"]     = &nu->decay.necm;
    ret["decay.nimpwt"]   = &nu->decay.nimpwt;

    // Set ancestor branches (ancestor is a vector)
    ret["ancestor.pdg"]     = &nu->ancestor[0].pdg;
    ret["ancestor.startx"]  = &nu->ancestor[0].startx;
    ret["ancestor.starty"]  = &nu->ancestor[0].starty;
    ret["ancestor.startz"]  = &nu->ancestor[0].startz;
    ret["ancestor.startt"]  = &nu->ancestor[0].startt;
    ret["ancestor.startpx"] = &nu->ancestor[0].startpx;
    ret["ancestor.startpy"] = &nu->ancestor[0].startpy;
    ret["ancestor.startpz"] = &nu->ancestor[0].startpz;
    ret["ancestor.stoppx"]  = &nu->ancestor[0].stoppx;
    ret["ancestor.stoppy"]  = &nu->ancestor[0].stoppy;
    ret["ancestor.stoppz"]  = &nu->ancestor[0].stoppz;
    ret["ancestor.polx"]    = &nu->ancestor[0].polx;
    ret["ancestor.poly"]    = &nu->ancestor[0].poly;
    ret["ancestor.polz"]    = &nu->ancestor[0].polz;
    ret["ancestor.pprodpx"] = &nu->ancestor[0].pprodpx;
    ret["ancestor.pprodpy"] = &nu->ancestor[0].pprodpy;
    ret["ancestor.pprodpz"] = &nu->ancestor[0].pprodpz;
    ret["ancestor.nucleus"] = &nu->ancestor[0].nucleus;
    ret["ancestor.proc"]    = &nu->ancestor[0].proc;
    ret["ancestor.ivol"]    = &nu->ancestor[0].ivol;
    ret["ancestor.imat"]    = &nu->ancestor[0].imat;

    // Set tgtexit branches
    ret["tgtexit.tvx"]    = &nu->tgtexit.tvx;
    ret["tgtexit.tvy"]    = &nu->tgtexit.tvy;
    ret["tgtexit.tvz"]    = &nu->tgtexit.tvz;
    ret["tgtexit.tpx"]    = &nu->tgtexit.tpx;
    ret["tgtexit.tpy"]    = &nu->tgtexit.tpy;
    ret["tgtexit.tpz"]    = &nu->tgtexit.tpz;
    ret["tgtexit.tptype"] = &nu->tgtexit.tptype;
    ret["tgtexit.tgen"]   = &nu->tgtexit.tgen;

    // Set traj branches (traj is a vector)
    ret["traj.trkx"]  = &nu->traj[0].trkx;
    ret["traj.trky"]  = &nu->traj[0].trky;
    ret["traj.trkz"]  = &nu->traj[0].trkz;
    ret["traj.trkpx"] = &nu->traj[0].trkpx;
    ret["traj.trkpy"] = &nu->traj[0].trkpy;
    ret["traj.trkpz"] = &nu->traj[0].trkpz;

    // Set other Dk2Nu top level branches
    ret["job"]    = &nu->job;
    ret["potnum"] = &nu->potnum;

    return ret;
  }

  //---------------------------------------------------------------------------
  std::vector<std::string> Wildcard(std::string fileWildcard)
  {
    std::vector<std::string> filelist;

    if( fileWildcard.find("root://") != std::string::npos ){
      //This path is using XROOTD. Pray that the user knows what they are doing.
      filelist.push_back(fileWildcard);
    }
    else{
      // First, expand environment variables and wildcards like the shell would
      wordexp_t p;
      wordexp(fileWildcard.c_str(), &p, WRDE_SHOWERR);



      for(unsigned int i = 0; i < p.we_wordc; ++i){
	// Check the file exists before adding it
	struct stat sb;
	if(stat(p.we_wordv[i], &sb) == 0)
	  filelist.push_back(p.we_wordv[i]); // Add the file
      }

      wordfree(&p); // Clean up
    }
    return filelist;
  }
}
