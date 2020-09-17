#ifndef MODEL_SPECTRUM_SPECTRUM_H_
#define MODEL_SPECTRUM_SPECTRUM_H_

#include <vector>
#include "peak.h"

namespace model {
namespace spectrum {

enum class SpectrumType
{ MS, EThcD, NONE };

class Spectrum
{
public:
    Spectrum() = default;
    Spectrum(const Spectrum& other){
        peaks_ = std::move(other.peaks_);
        scan_num_ = other.scan_num_;
        type_ = other.type_;
        precursor_mz_ = other.precursor_mz_;
        precursor_charge_ = other.precursor_charge_;
    }

    Spectrum& operator=(const Spectrum& other){
        peaks_ = std::move(other.peaks_);
        scan_num_ = other.scan_num_;
        type_ = other.type_;
        precursor_mz_ = other.precursor_mz_;
        precursor_charge_ = other.precursor_charge_;
        return *this;
    }

    int Scan() { return scan_num_; }
    void set_scan(int num) { scan_num_ = num; }
    SpectrumType Type() { return type_; }
    void set_type(SpectrumType type) { type_ = type; }

    std::vector<Peak>& Peaks() { return peaks_; }
    void set_peaks(std::vector<Peak>& peaks) 
        { peaks_ = std::move(peaks); }

    double PrecursorMZ() { return precursor_mz_; }
    double PrecursorCharge() { return precursor_charge_; }

    void set_parent_mz(double mz) { precursor_mz_ = mz;}
    void set_parent_charge(int charge) { precursor_charge_ = charge; }

protected:
    std::vector<Peak> peaks_;
    int scan_num_;
    SpectrumType type_;
    double precursor_mz_;
    int precursor_charge_;

};

}  //  namespace spectrum
}  //  namespace model

#endif