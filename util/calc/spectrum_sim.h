#ifndef UTIL_CAL_SPECTRUM_SIM_H
#define UTIL_CAL_SPECTRUM_SIM_H

#include <cmath>
#include <algorithm>
#include "../../model/spectrum/spectrum.h"

namespace util {
namespace calc {

class SpectrumSim
{
public:
    SpectrumSim(): tolerance_(0.01){}
        double ComputeCosine(std::vector<model::spectrum::Peak> spec, std::vector<model::spectrum::Peak> other)
    {
        model::spectrum::Peak lower1 = *std::min_element(spec.begin(), spec.end());
        model::spectrum::Peak lower2 = *std::min_element(other.begin(), other.end());
        model::spectrum::Peak upper1 = *std::max_element(spec.begin(), spec.end());
        model::spectrum::Peak upper2 = *std::max_element(other.begin(), other.end());

        double lower = std::min(lower1.MZ(), lower2.MZ());
        double upper = std::max(upper1.MZ(), upper2.MZ());

        int bucket_size = (int) ((upper - lower) / tolerance_ + 1);

        std::vector<std::vector<model::spectrum::Peak>> q1, q2;
        q1.assign(bucket_size, std::vector<model::spectrum::Peak>());
        q2.assign(bucket_size, std::vector<model::spectrum::Peak>());

        for(auto& pk : spec)
        {
            int index = Index(pk, lower);
            q1[index].push_back(pk);
        }
        for(auto& pk : other)
        {
            int index = Index(pk, lower);
            q2[index].push_back(pk);
        }

        double numerator = 0;
        for (int i = 0; i < bucket_size; i++)
        {
            numerator += DotProduct(q1[i], q2[i]);
        }

        double denominator1 = 0;
        for (auto& pk : spec)
        {
            denominator1 += pk.Intensity() * pk.Intensity();
        }
        double denominator2 = 0;
        for (auto& pk : other)
        {
            denominator2 += pk.Intensity() * pk.Intensity();
        }
        double denominator = sqrt(denominator1) * sqrt(denominator2);
        return numerator / denominator;
    }

    double ComputeCosine(model::spectrum::Spectrum spec, model::spectrum::Spectrum other)
    {
        return ComputeCosine(spec.Peaks(), other.Peaks());
    }

    double Tolerance() { return tolerance_; }
    void set_tolerance(double tol) { tolerance_ = tol; }

protected:
    double DotProduct(std::vector<model::spectrum::Peak>& q1, std::vector<model::spectrum::Peak>& q2)
    {
        if (q1.empty() || q2.empty())
            return 0;
        
        std::vector<model::spectrum::Peak> p1, p2;
        p1.assign(q1.begin(), q1.end());
        p2.assign(q2.begin(), q2.end());

        // sort by intesntity, take min size, sort by mz
        if (p1.size() != p2.size())
        {
            std::sort(p1.begin(), p1.end(), IntensityGreater);
            std::sort(p2.begin(), p2.end(), IntensityGreater);
            if (p1.size() > p2.size())
                p1.erase(p1.begin() + p2.size());
            else
                p2.erase(p2.begin() + p1.size());

            std::sort(p1.begin(), p1.end(), MZComp);
            std::sort(p2.begin(), p2.end(), MZComp);
        }

        // pairwise multiple sum
        double numerator = 0;
        for (size_t i = 0; i < p1.size(); i++)
        {
            numerator += (p1[i].Intensity() * p2[i].Intensity());
        }
        return numerator;
    }
    static bool IntensityGreater (model::spectrum::Peak& i, model::spectrum::Peak& j) 
        { return (i.Intensity() > j.Intensity()); }
    static bool MZComp (model::spectrum::Peak& i, model::spectrum::Peak& j) 
        { return (i.MZ() < j.MZ()); }

    int Index(model::spectrum::Peak& pk, double lower) 
        { return (pk.MZ() - lower) / tolerance_; } 
    double tolerance_;

};


}  //  namespace calc
}  //  namespace util

#endif