#include "Weight.h"
// This include forces the code in "Weight.h" to get built during compilation

namespace flxrd
{
  // A constant weight of value c, which must be specified when the object is called
  const Weight kConstant(double c)
  {
    return Weight({},
                  [c](const double&, const bsim::Dk2Nu*, const int&, const TObject*)
                  { return c; });
  }
}
