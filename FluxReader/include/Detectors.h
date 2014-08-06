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

  // NOvA Near Detector
  static Detector knova_nd("nova_nd", "CH2",
                           1141.4, -345.6, 99466.5,     // position
                           262.14, 393.27, 1424.52698,  // size
                           1);

  // NOvA Far Detector
  static Detector knova_fd("nova_fd", "CH2",
                           1103746., -416264, 81042232., // position
                           1560.,    1560.,   7800.,     // size
                           1);

  // NOvA NDOS (IPND)
  static Detector knova_ipnd("nova_ipnd", "CH2",
                             -29.,   9221.,  84176.,      // position
                             262.14, 393.27, 1424.52698,  // size
                             1);

  // MINOS Near Detector
  static Detector kminos_nd("minos_nd", "CH2",
                            0.,   0.,   104000., // position
                            100., 100., 500.,    // size
                            1);

  // MINOS Far Detector
  static Detector kminos_fd("minos_fd", "CH2",
                            0.,   0.,   73534000., // position
                            374., 374., 2800.,     // size
                            1);
}
