#include "Utilities.h"

// C/C++ Includes
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
    // First, expand environment variables and wildcards like the shell would
    wordexp_t p;
    wordexp(fileWildcard.c_str(), &p, WRDE_SHOWERR);

    std::vector<std::string> filelist;

    for(unsigned int i = 0; i < p.we_wordc; ++i){
      // Check the file exists before adding it
      struct stat sb;
      if(stat(p.we_wordv[i], &sb) == 0)
        filelist.push_back(p.we_wordv[i]); // Add the file
    }

    wordfree(&p); // Clean up

    return filelist;
  }
}
