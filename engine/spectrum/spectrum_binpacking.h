#ifndef ENGINE_SPECTRUM_BINPACKING_H
#define ENGINE_SPECTRUM_BINPACKING_H
#include <vector>
#include <algorithm>  
#include "../../model/spectrum/spectrum.h"
#include "../../algorithm/base/binpacking.h"

namespace engine {
namespace spectrum {
 
class SpectrumBinPacking : 
    public algorithm::base::BinPacking<model::spectrum::Peak>
{
public:
    SpectrumBinPacking(double tol, double lower, double upper):
        BinPacking(tol, lower, upper){}
        
    std::vector<double> Packing
        (model::spectrum::Spectrum& spec);

protected:
    virtual double Merge
        (std::vector<model::spectrum::Peak>& peak);

    double Position(const 
        model::spectrum::Peak& pk) const override 
    { 
        return pk.MZ(); 
    }
};

} // namespace spectrum
} // namespace engine

#endif