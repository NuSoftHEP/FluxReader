#pragma once

// Package Includes
#include "Detector.h"

namespace flxrd
{
  // Null Detector
  static Detector kNullDetector("", "",
                                0., 0., 0.,
                                0., 0., 0.,
                                0);

  // Detector locations from https://cdcvs.fnal.gov/redmine/projects/dk2nu/repository/entry/trunk/dk2nu/etc/locations.txt

  // MicroBooNE
  static Detector kMicroBooNE("MicroBooNE", "CH2",
                              5300., 7600., 67900., // position
                              0., 0., 0.,           // size
                              1);

  // Minerva
  static Detector kMinerva("Minerva", "CH2",
                           -56.28, -53.29317, 103231.9, // position
                           0., 0., 0.,                  // size
                           1);

  // MiniBooNE
  static Detector kMiniBooNE("MiniBooNE", "CH2",
                             2604., 7864., 74487., // position
                             0., 0., 0.,           // size
                             1);

  // MINOS Near Detector
  static Detector kMINOS_ND("MINOS-ND", "CH2",
                            0.,   0.,   103648.8, // position
                            100., 100., 500.,     // size
                            1);

  // MINOS Far Detector
  static Detector kMINOS_FD("MINOS-FD", "CH2",
                            0.,   0.,   73534000., // position
                            374., 374., 2800.,     // size
                            1);

  // NOvA Near Detector
  static Detector kNOvA_ND("NOvA-ND", "CH2",
                           1150.172, -280.0719, 100099.2, // position
                           262.14, 393.27, 1424.52698,    // size
                           1);

  // NOvA Far Detector
  static Detector kNOvA_FD("NOvA-FD", "CH2",
                           1103746., -416264, 81042232., // position
                           1560.,    1560.,   7800.,     // size
                           1);

  // NOvA NDOS (IPND)
  static Detector kNOvA_IPND("NOvA-IPND", "CH2",
                             -29.,   9221.,  84176.,      // position
                             262.14, 393.27, 1424.52698,  // size
                             1);

  // SciBooNE
  static Detector kSciBooNE("SciBooNE", "CH2",
                            19760., 5340., 33940., // position
                            0., 0., 0.,            // size
                            1);
}
