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

  /*const char* srt_priv_path = gSystem->ExpandPathName("$(SRT_PRIVATE_CONTEXT)");
  if(!srt_priv_path) {
    std::cout << "DK2NU needs to be setup." << std::endl;
    return -3;
  }*/

  /*const char* srt_pub_path = gSystem->ExpandPathName("$(SRT_PRIVATE_CONTEXT)");
  if(!srt_pub_path) {
    std::cout << "SRT_PUBLIC_CONTEXT needs to be setup." << std::endl;
    return -4;
  }*/

  TString include_path = gSystem->GetIncludePath(); // Get the current list of include paths

  include_path += " -I";
  include_path += root_path;
  include_path += "/include";

  include_path += " -I";
  include_path += dk2nu_path;
  include_path += "/tree";

  include_path += " -I";
  include_path += "/nova/app/users/gkafka/test_svn_development/FluxReader";
  include_path += "/include";

  // Only add a path to the SRT_PRIVATE_CONTEXT if it is setup
  /*if(srt_priv_path) {
    include_path += " -I";
    include_path += srt_priv_path;
    include_path += "/include";
  }*/

  // Add SRT_PUBLIC_CONTEXT path
  // Note that this MUST come AFTER SRT_PRIVATE_CONTEXT
  /*include_path += " -I";
  include_path += srt_pub_path;
  include_path += "/include";*/

  gSystem->SetIncludePath(include_path); // Set the include path as the new, appended list of include paths

  gSystem->Load("lib/libFluxReader.so"); // This library must be loaded before running a compiled script

  return 0;
}
