#pragma once

// Package Includes
#include "Detector.h"
#include "Utilities.h"

namespace flxrd
{
  // Null Detector
  static Detector kNullDetector("", "",
                                0., 0., 0., // position
                                0., 0., 0., // size
                                0);

  // Detector names, in the call to LoadDetCoords,
  // must match EXACTLY those found at ${DK2NU}/etc/locations.txt
  // The detector name in the first constructor argument can be anything

  // MicroBooNE
  static Detector kMicroBooNE("MicroBooNE", "Ar",
                              LoadDetCoords("MicroBooNE"),
                              {0., 0., 0.},
                              1);

  // Minerva
  static Detector kMinerva("Minerva", "CH2",
                           LoadDetCoords("Minerva"),
                           {0., 0., 0.},
                           1);

  // MiniBooNE
  static Detector kMiniBooNE("MiniBooNE", "CH2",
                             LoadDetCoords("MiniBooNE"),
                             {0., 0., 0.},
                             1);

  // MINOS Near Detector
  static Detector kMINOS_ND("MINOS-ND", "Fe",
                            LoadDetCoords("MINOS NearDet"),
                            {100., 100., 500.},
                            1);

  // MINOS Far Detector
  static Detector kMINOS_FD("MINOS-FD", "Fe",
                            LoadDetCoords("MINOS FarDet"),
                            {374., 374., 2800.},
                            1);

  // NOvA Near Detector
  static Detector kNOvA_ND("NOvA-ND", "CH2",
                           LoadDetCoords("NOvA NearDet"),
                           {262.14, 393.27, 1424.52698},
                           1);

  // NOvA Far Detector
  static Detector kNOvA_FD("NOvA-FD", "CH2",
                           LoadDetCoords("NOvA FarDet"),
                           {1560., 1560., 7800.},
                           1);

  // NOvA NDOS (IPND)
  static Detector kNOvA_IPND("NOvA-IPND", "CH2",
                             -29.,   9221.,  84176.,
                             262.14, 393.27, 1424.52698,
                             1);

  // SciBooNE
  static Detector kSciBooNE("SciBooNE", "CH2",
                            LoadDetCoords("SciBooNE"),
                            {0., 0., 0.},
                            1);
}
