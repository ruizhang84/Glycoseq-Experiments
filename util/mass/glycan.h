#ifndef UTIL_MASS_GLYCAN_H
#define UTIL_MASS_GLYCAN_H

#include "../../model/glycan/glycan.h"

namespace util {
namespace mass {

class GlycanMass
{
public:
    static double Compute(const model::glycan::Glycan& glycan) 
    {
        double ms = 0;
        for(const auto &it : glycan.CompositionConst())
        {
            switch (it.first)
            {
            case model::glycan::Monosaccharide::GlcNAc:
                ms += kHexNAc * it.second;
                break;
            case model::glycan::Monosaccharide::Gal:
                ms += kHex * it.second;
                break;     
            case model::glycan::Monosaccharide::Man:
                ms += kHex * it.second;
                break;     
            case model::glycan::Monosaccharide::Fuc:
                ms += kFuc * it.second;
                break;     
            case model::glycan::Monosaccharide::NeuAc:
                ms += kNeuAc * it.second;
                break;           
            case model::glycan::Monosaccharide::NeuGc:
                ms += kNeuGc * it.second;
                break;    
            default:
                break;
            }
        }
        return ms;
    }

    static double Compute
        (const std::map<model::glycan::Monosaccharide, int>& composite) 
    {
        double ms = 0;
        for(const auto &it : composite)
        {
            switch (it.first)
            {
            case model::glycan::Monosaccharide::GlcNAc:
                ms += kHexNAc * it.second;
                break;
            case model::glycan::Monosaccharide::Gal:
                ms += kHex * it.second;
                break;     
            case model::glycan::Monosaccharide::Man:
                ms += kHex * it.second;
                break;     
            case model::glycan::Monosaccharide::Fuc:
                ms += kFuc * it.second;
                break;     
            case model::glycan::Monosaccharide::NeuAc:
                ms += kNeuAc * it.second;
                break;           
            case model::glycan::Monosaccharide::NeuGc:
                ms += kNeuGc * it.second;
                break;    
            default:
                break;
            }
        }
        return ms;
    }

    static constexpr double kHexNAc = 203.0794;
    static constexpr double kHex = 162.0528;
    static constexpr double kFuc = 146.0579;
    static constexpr double kNeuAc = 291.0954;
    static constexpr double kNeuGc = 307.0903;
    static constexpr double kWater = 18.0105;
};



} // namespace mass
} // namespace util

#endif