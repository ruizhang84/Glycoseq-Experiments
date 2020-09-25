#ifndef ENGINE_ANALYSIS_FDR_FILTER_H
#define ENGINE_ANALYSIS_FDR_FILTER_H

#include <vector>
#include <algorithm>
#include "../search/spectrum_search.h"


namespace engine {
namespace analysis {

class FDRFilter
{
public:
    FDRFilter(double fdr): fdr_(fdr), cutoff_(-1){}

    virtual void Init()
    {
         CutoffInit();
    }
    
    std::vector<engine::search::SearchResult>& Target() { return target_; }
    std::vector<engine::search::SearchResult>& Decoy() { return decoy_; }
    double Cutoff() const { return cutoff_; }
    void set_data(std::vector<engine::search::SearchResult>& targets, 
        std::vector<engine::search::SearchResult>& decoys) 
            { target_ = targets; decoy_ = decoys; }
    void set_cutoff(double cutoff) { cutoff_ = cutoff; }

    virtual std::vector<engine::search::SearchResult> Filter()
    {
        std::vector<engine::search::SearchResult> res;
        if (cutoff_ < 0) return target_;

        for(const auto& it : target_)
        {
            if (ComputeScore(it) >= cutoff_)
            {
                res.push_back(it);
            }
        }
        if (!res.empty())
            std::sort(res.begin(), res.end(), ScanOrder);
        return res;
    }

protected:
    static double ComputeScore(const engine::search::SearchResult& result)
        { return result.RawScore(); }

    static bool ScoreLess(const engine::search::SearchResult& r1, const engine::search::SearchResult& r2) 
        { return (ComputeScore(r1) < ComputeScore(r2)); }

    static bool ScanOrder(const engine::search::SearchResult& r1, const engine::search::SearchResult& r2)
        { return r1.Scan() < r2.Scan(); }

    virtual void CutoffInit()
    {
        // init
        cutoff_ = -1;
        if (decoy_.size() == 0 || target_.size() == 0 || 
            (decoy_.size() * 1.0 / (decoy_.size() + target_.size()) < fdr_))   //trivial case
        {
            return;
        }

        std::sort(target_.begin(), target_.end(), ScoreLess);
        std::sort(decoy_.begin(), decoy_.end(), ScoreLess);

        // compare and compute
        int i = 0, j = 0;
        int target_size = (int) target_.size();
        int decoy_size = (int) decoy_.size();
        while (i < target_size)
        {
            // decoy score is no less than targets
            while (j < decoy_size && 
                ComputeScore(decoy_[j]) <= ComputeScore(target_[i]))
            {
                j++;
            }
            // compute fdr rate
            double rate = (decoy_size - j ) * 1.0 / (target_size + decoy_size - i - j);
            if (rate <= fdr_)
            {
                cutoff_ = ComputeScore(target_[i]);
                return;
            }
            else
            {
                i++;
            }
        }

        // set max
        cutoff_ = INT64_MAX;
    }

    double fdr_;
    double cutoff_;
    std::vector<engine::search::SearchResult> target_;
    std::vector<engine::search::SearchResult> decoy_;

};

} // namespace analysis
} // namespace engine 

#endif