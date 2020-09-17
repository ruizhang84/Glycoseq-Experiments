#ifndef ALGORITHM_BINARY_SEARCH_H
#define ALGORITHM_BINARY_SEARCH_H

#include <vector>
#include <cstdlib> 
#include <algorithm>
#include "search.h"
#include "../../util/mass/spectrum.h"

namespace algorithm {
namespace search {

// enhanced binary search
class BinarySearch
{
public:
    BinarySearch(double tol, ToleranceBy by):
        tolerance_(tol), by_(by), base_(-1), scale_(1) {};

    virtual void Init() 
    {
        if (!data_.empty())
            std::sort(data_.begin(), data_.end());
    }
    double Tolerance() const { return tolerance_; }
    std::vector<double>& Data() { return data_; }
    ToleranceBy ToleranceType() const { return by_; }
    double Base() const { return base_; }
    double Scale() const { return scale_; }
    void set_tolerance(double tol) { tolerance_ = tol; }
    void set_tolerance_by(ToleranceBy by) { by_ = by; }
    void set_data(std::vector<double> data) { data_ = data; }
    void set_base(double base) { base_ = base; }
    void set_scale(double scale) { scale_ = scale; }

    virtual bool Search(const double target)
    {
        if (data_.empty()) 
            return false;

        int start = 0, end = data_.size()-1;
        while (start <= end)
        {
            int mid = (end - start) / 2 + start;
            if (Match(data_[mid], target))
                return true;
            else if (data_[mid] < target)
                start = mid + 1;
            else
                end = mid - 1;
        }
        return false;
    }

protected:
    virtual bool Match(const double p, const double target)
    {
        switch (by_)
        {
        case ToleranceBy::PPM:
            if (base_ < 0)
                return util::mass::SpectrumMass::ComputePPM(p, target) < tolerance_;
            else
                return std::abs(p - target) / base_ * 1000000.0 < tolerance_; 
        case ToleranceBy::Dalton:
            return std::abs(p - target) < tolerance_ * scale_;
        default:
            break;
        }
        return false;
    }

    double tolerance_; 
    ToleranceBy by_;
    std::vector<double> data_;
    double base_;
    double scale_;
};

} // namespace algorithm
} // namespace search 

#endif
