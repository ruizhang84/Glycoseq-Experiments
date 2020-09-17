#ifndef ALGORITHM_BASE_BINPACKING_H
#define ALGORITHM_BASE_BINPACKING_H
#include <vector>
#include <cmath> 

namespace algorithm {
namespace base {

template <class T>
class BinPacking
{
public:
    BinPacking(double tol, double lower, double upper):
        tolerance_(tol), lower_(lower), upper_(upper) {};

    std::vector<std::vector<T>> Packing(
        const std::vector<T>& vect) const
    {
        std::vector<std::vector<T>> bins;
        bins.assign(Bucket(), std::vector<T>());
        for(auto& it : vect)
        {
            int index = Index(Position(it));
            if (index >= 0)
                bins[index].push_back(it);
        }
        return bins;
    }

    double Tolerance() { return tolerance_; }
    int BinSize() { return Bucket(); }
    void set_lower(double lower) { lower_ = lower; }
    void set_upper(double upper) { upper_ = upper; }
    void set_tolerance(double tol) { tolerance_ = tol; }

protected:
    virtual int Index(double pos) const
    { 
        if (pos < lower_ || pos > upper_)
            return -1;
        return (int) floor((pos - lower_) / tolerance_); 
    }
    virtual double Position(const T& elem) const = 0;
    virtual int Bucket() const
        { return (int) ceil((upper_ - lower_ + 1) / tolerance_); }

    double tolerance_;
    double lower_;
    double upper_;

};

} // namespace base
} // namespace algorithm

#endif