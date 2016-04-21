int load_flxrd()
{
  gSystem->Load("libCintex.so");
  Cintex::Enable();

  // Get all of the environment variable paths

  const char* root_path = gSystem->ExpandPathName("$(ROOTSYS)");
  if(!root_path) {
    std::cout << "ROOT needs to be setup." << std::endl;
    return -1;
  }

  const char* dk2nu_path = gSystem->ExpandPathName("$(DK2NU)");
  if(!dk2nu_path) {
    std::cout << "DK2NU needs to be setup." << std::endl;
    return -2;
  }

  const char* fluxreader_priv_path = gSystem->ExpandPathName("$(FLUXREADER_PRIV)");
  if(!fluxreader_priv_path) {
    std::cout << "FluxReader needs to be setup." << std::endl;
    return -3;
  }

  /*const char* fluxreader_pub_path = gSystem->ExpandPathName("$(FLUXREADER_PUB)");
  if(!fluxreader_pub_path) {
    std::cout << "FluxReader needs to be setup." << std::endl;
    return -4;
  }*/

  // Get the path to the load_dk2nu.C script within the Dk2Nu package
  TString load_dk2nu_path = dk2nu_path;
  load_dk2nu_path += "/scripts";
  load_dk2nu_path += "/load_dk2nu.C";

  // Create a command that can be run
  TString load_dk2nu_command = ".L " + load_dk2nu_path;

  // Process load_dk2nu.C and quit if there is an error
  if(gROOT->ProcessLine(load_dk2nu_command) != 0) {
    std::cout << "Error occurred in running " << load_dk2nu_path << "." << std::endl;
    return -5;
  }

  TString include_path = gSystem->GetIncludePath(); // Get the current list of include paths

  // Add on all of the other necessary include directory locations

  include_path += " -I";
  include_path += root_path;
  include_path += "/include";

  include_path += " -I";
  include_path += dk2nu_path;
  include_path += "/include/dk2nu/tree";

  include_path += " -I";
  include_path += fluxreader_priv_path;
  include_path += "/include";

  // Note that this MUST come AFTER the private path!
  /*include_path += " -I";
  include_path += fluxreader_pub_path;
  include_path += "/include";*/

  gSystem->SetIncludePath(include_path); // Set the include path as the new, appended list of include paths

  gSystem->Load("lib/libFluxReader.so"); // This library must be loaded before running a compiled script

  return 0;
}
