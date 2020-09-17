#ifndef ALGORITHM_SEARCH_H
#define ALGORITHM_SEARCH_H

#include <vector>
#include <memory>
#include <cstdlib>
#include <algorithm> 
#include "point.h"
#include "../../util/mass/spectrum.h"

namespace algorithm {
namespace search {

enum class ToleranceBy { PPM, Dalton}; 

template <class T>
class BasicSearch
{
typedef std::vector<std::shared_ptr<Point<T>>> Points;
public:
    BasicSearch(double tol, ToleranceBy by):
        tolerance_(tol), by_(by), base_(-1), scale_{1} {};

    virtual void Init() 
    {
        if (!data_.empty())
            std::sort(data_.begin(), data_.end(), BasicComp);
    }

    double Tolerance() const { return tolerance_; }
    Points& Data() { return data_; }
    ToleranceBy ToleranceType() const { return by_; }
    double Base() const { return base_; }
    double Scale() const { return scale_; }
    void set_tolerance(double tol) { tolerance_ = tol; }
    void set_tolerance_by(ToleranceBy by) { by_ = by; }
    void set_data(std::vector<std::shared_ptr<Point<T>>> data)
        { data_ = data; }
    void set_base(double base) { base_ = base; }
    void set_scale(double scale) { scale_ = scale; }
    
    virtual std::vector<T> Query(const double target)
    {
        std::vector<T> result;
        if (data_.empty()) 
            return result;

        int start = 0, end = data_.size()-1;
        while (start <= end)
        {
            int mid = (end - start) / 2 + start;
            if (Match(data_[mid].get(), target))
            {
                for(int left = mid; left >= 0 && Match(data_[left].get(), target); left--)
                {
                    result.push_back(data_[left]->Content());
                }

                for (int right = mid+1; right < (int) data_.size() && Match(data_[right].get(), target); right++)
                {
                    result.push_back(data_[right]->Content());
                }
                break;
            }
            else if (data_[mid]->Value() < target)
                start = mid + 1;
            else
                end = mid - 1;
        }
        return result;
    }

    virtual bool Search(const double target)
    {
        if (data_.empty()) 
            return false;

        int start = 0, end = data_.size()-1;
        while (start <= end)
        {
            int mid = (end - start) / 2 + start;
            if (Match(data_[mid].get(), target))
                return true;
            else if (data_[mid]->Value() < target)
                start = mid + 1;
            else
                end = mid - 1;
        }
        return false;
    }

protected:
    virtual bool Match(const Point<T>* p, const double target)
    {
        double diff; 
        switch (by_)
        {
        case ToleranceBy::PPM:
            if (base_ < 0)
                diff = util::mass::SpectrumMass::ComputePPM(p->Value(), target);
            else
                diff = std::abs(p->Value() - target) / base_ * 1000000.0;
            return diff < tolerance_;
        case ToleranceBy::Dalton:
            diff = p->Value() - target;
            return std::abs(diff) < tolerance_ * scale_;
        default:
            break;
        }
        return false;
    }
    static bool BasicComp(const std::shared_ptr<Point<T>>& p1, const std::shared_ptr<Point<T>>& p2)
        { return p1->Value() < p2->Value(); }

    double tolerance_; 
    ToleranceBy by_;
    Points data_;
    double base_;
    double scale_;  // due to charge when compare mass
};

} // namespace algorithm
} // namespace search 

#endif
