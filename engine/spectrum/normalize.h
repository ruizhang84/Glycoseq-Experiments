#ifndef ENGINE_SPECTRUM_NORMALIZE_H
#define ENGINE_SPECTRUM_NORMALIZE_H

#include <vector>
#include <algorithm>
#include <numeric> 
#include "../../model/spectrum/spectrum.h"

namespace engine {
namespace spectrum {
 
class Normalizer
{
public:
    static void Transform(model::spectrum::Spectrum& spec)
    {   
        Transform(spec.Peaks());
    }   

    //normalization on total ion intensity sums 
    static void Transform(std::vector<model::spectrum::Peak>& peaks)
    {
        double sum = 0;
        for(auto& it : peaks)
        {
            sum += it.Intensity();
        }
        for(auto& it : peaks)
        {
            it.set_intensity(it.Intensity() / sum * 100.0);
        }
    }

protected:
    static bool IntensityCmp (const model::spectrum::Peak& i, const model::spectrum::Peak& j) 
        { return (i.Intensity() < j.Intensity()); }

};

} // namespace spectrum
} // namespace engine

#endif