#include "spectrum_binpacking.h"

namespace engine {
namespace spectrum {

using namespace model::spectrum;

std::vector<double> SpectrumBinPacking::Packing
    (model::spectrum::Spectrum& spec)
{
    
    std::vector<double> result;
    std::vector<Peak>& peaks = spec.Peaks();
    std::vector<std::vector<Peak>> peak_bins = BinPacking::Packing(peaks); 


    for(size_t i = 0; i < peak_bins.size(); i++)
    {
       
        result.push_back(Merge(peak_bins[i]));
    }
    return result;
}

double SpectrumBinPacking::Merge
    (std::vector<model::spectrum::Peak>& peak)
{
    if (peak.empty())
        return 0;

    return std::max_element(peak.begin(), peak.end())->Intensity(); 
}

} // namespace spectrum
} // namespace engine