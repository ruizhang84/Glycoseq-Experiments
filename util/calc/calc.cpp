#include "calc.h"

namespace util {
namespace calc {

const double Calc::DotProduct(
    const std::vector<double>& vect1, 
    const std::vector<double>& vect2) const
{
    return std::inner_product(vect1.begin(), vect1.end(),
            vect2.begin(), 0.0);
}

} // namespace calc
} // namespace util


