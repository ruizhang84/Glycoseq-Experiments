#ifndef ENGINE_SEARCH_SEARCH_SCORE_H
#define ENGINE_SEARCH_SEARCH_SCORE_H


#include <vector>
#include <map>
#include <cmath> 
#include <numeric>
#include "../../model/spectrum/spectrum.h"
#include "../../util/mass/glycan.h"
#include "../../util/mass/peptide.h"
#include "../../util/mass/spectrum.h"
#include <iostream>

namespace engine{
namespace search{

enum class SearchType { Core, Branch, Terminal, Oxonium, Peptide, Base };
enum class ScoreType { Precursor, Elution };

class SearchResult
{
public:
    SearchResult() = default;

    const int Scan() const { return scan_; }
    const int ModifySite() const { return pos_; }
    const std::string Sequence() const { return peptide_; }
    const std::string Glycan() const { return glycan_; }
    const double RawScore() const 
    { 
        if (score_.size() == 0) return 0.0;
        if (simple_)
            return std::accumulate(score_.begin(), score_.end(), 0.0);
        return std::sqrt(std::accumulate(score_.begin(), score_.end(), 0.0));
    }
    const double Value() const { return value_; }
    const std::vector<double> Score() const { return score_; }

    void set_simple(bool simple) { simple_ = simple; }

    const double ExtraScore(ScoreType type) const 
    {
        const auto& it = extra_.find(type);
        if (it == extra_.end())
            return 0;
        return it->second;
    }

    void set_scan(int scan) { scan_ = scan; }
    void set_site(int pos) { pos_ = pos; }
    void set_peptide(std::string seq) { peptide_ = seq; }
    void set_glycan(std::string glycan) { glycan_ = glycan; }
    void set_score(std::vector<double> score) { score_ = score; }
    void set_value(double value) { value_ = value; }
    void set_extra(double score, ScoreType type) { extra_[type] = score; }

    static double PeakValue(const std::vector<model::spectrum::Peak>& peaks, bool simple=true)
    { 
        double sum = 0;
        for(const auto& it : peaks)
        {
            if (simple)
                sum += it.Intensity();
            else
                sum += it.Intensity() * it.Intensity();
        }
        return sum;
    }
    
    static double PrecursorValue(const std::string peptide, const std::string composite,
        double precursor_mass, double isotopic)
    {
        double mass = util::mass::PeptideMass::Compute(peptide)
            + util::mass::GlycanMass::Compute(model::glycan::Glycan::Interpret(composite));
        
        double ppm = kPPM;
        for (int i = 0; i <= isotopic; i ++)
        {
            double isotopic_mass = util::mass::SpectrumMass::kIon * i + mass;
            double isotopic_ppm = util::mass::SpectrumMass::ComputePPM(isotopic_mass, precursor_mass);
            ppm = (ppm <= isotopic_ppm) ? ppm : isotopic_ppm;
        }
        return 1.0 - ppm / kPPM;
    }
    static constexpr double kPPM = 50.0;

protected:
    bool simple_ = false;
    int scan_;
    std::string peptide_;
    std::string glycan_;
    int pos_;
    std::vector<double> score_;
    double value_;
    std::map<ScoreType, double> extra_;
    
};


class ResultCollector{
public:
    ResultCollector() = default;

    void set_score_compute(bool simple){
        simple_ = simple;
    }
    
    std::vector<SearchResult> Result()
    {
        // remove if hits over 20
        if ((int) results_.size() > max_hits)
        {
            std::sort(results_.begin(), results_.end(), 
                [](const SearchResult& r1, const SearchResult& r2) -> bool { return r1.RawScore() > r1.RawScore(); });
            results_.erase(results_.begin() + max_hits, results_.end());
        }

        return results_;
    }
    std::vector<SearchResult> BestResult()
    {
        // max score
        std::vector<SearchResult> best_rest;
        double max_score = 0;
        for (const auto& it : results_)
        {
            double score = it.RawScore();
            if (score >= max_score){
                if (score > max_score){
                    best_rest.clear();
                }
                best_rest.push_back(it);
                max_score = score;
            }
        }
        // update extra
        for (auto& it : best_rest)
        {
            double score = SearchResult::PrecursorValue(
                it.Sequence(), it.Glycan(), precursor_mass_, isotopic_);
            it.set_extra(score, ScoreType::Precursor);
        }
        // pick tie by extra
        std::vector<SearchResult> res;
        max_score = 0;
        for (const auto& it : best_rest)
        {
            double score = it.ExtraScore(ScoreType::Precursor);
            if (score >= max_score){
                if (score > max_score){
                    res.clear();
                } 
                res.push_back(it);
                max_score = score;
            }
        }
        return res;
    }
    void Update(int scan, const std::string& sequence, const std::string& composite)
    {
        for(const auto& pos_it : peptide_)
        {
            // compute score
            std::vector<double> score_vec = ComputeScore(pos_it.second);
            // emplace results
            Emplace(scan, sequence, composite, pos_it.first, score_vec);
        }
    }

    void BestUpdate(int scan, const std::string& sequence, const std::string& composite)
    {
        for(const auto& pos_it : peptide_)
        {
            // compute score
            std::vector<double> score_vec = ComputeScore(pos_it.second);
            double score = std::accumulate(score_vec.begin(), score_vec.end(), 0.0);
            if (score >= best_)
            {
                if (score > best_)
                    results_.clear();
                best_ = score;
                // emplace results
                Emplace(scan, sequence, composite, pos_it.first, score_vec);
            }
        }
    }
    
    void SpectrumBase(const std::vector<model::spectrum::Peak>& spectrum_peaks)
    {
        spectrum_ = SearchResult::PeakValue(spectrum_peaks, simple_);
    }
    void InitCollect()
    {
        peptide_.clear(); 
        glycan_core_.clear();
        glycan_branch_.clear();
        glycan_terminal_.clear();
    }

    void OxoniumCollect(const std::vector<model::spectrum::Peak>& oxonium_peaks)
    {
        if (oxonium_peaks.empty()) return;
        oxonium_ = SearchResult::PeakValue(oxonium_peaks, simple_);
    }
    void PeptideCollect(const std::vector<model::spectrum::Peak>& peptide_peaks, int pos)
    {
        if (!peptide_peaks.empty())
        {  
            peptide_[pos] = SearchResult::PeakValue(peptide_peaks, simple_);
        }
    }
    void GlycanCollect(const std::vector<model::spectrum::Peak>& glycan_peaks, 
        std::string isomer, SearchType type)
    {
        if (!glycan_peaks.empty())
        {
            switch (type)
            {
                case SearchType::Core:
                    glycan_core_[isomer] = SearchResult::PeakValue(glycan_peaks, simple_);
                    break;
                case SearchType::Branch:
                    glycan_branch_[isomer] = SearchResult::PeakValue(glycan_peaks, simple_);
                case SearchType::Terminal:
                    glycan_terminal_[isomer] = SearchResult::PeakValue(glycan_peaks, simple_);
                default:
                    break;
            }
        }
    }
    void PrecursorCollect(double precursor_mass, int isotopic)
    {
        precursor_mass_ = precursor_mass;
        isotopic_ = isotopic;
    }
    bool OxoniumMiss() { return oxonium_ <= 0; }
    bool PeptideMiss()
    {
        return peptide_.empty();
    }
    bool GlycanMiss(const std::string& isomer)
    {
        return (glycan_core_.find(isomer) == glycan_core_.end());
    }
    bool GlycanMiss()
    {
        return glycan_core_.empty();
    }
    bool Empty() { return results_.empty(); }

protected:
    std::vector<double> ComputeScore(double peptide_score)
    {
        double score = 0;
        std::vector<double> score_vec(5, 0.0);
        for(const auto& isomer_it : glycan_core_)
        {
            std::string isomer = isomer_it.first;
            double glycan_score = glycan_core_[isomer] + glycan_branch_[isomer] + glycan_terminal_[isomer]; 
            if (glycan_score > score)
            {
                score = glycan_score;
                score_vec[0] = glycan_core_[isomer];
                score_vec[1] = glycan_branch_[isomer];
                score_vec[2] = glycan_terminal_[isomer];                
            }
        }
        score_vec[3] = oxonium_;
        score_vec[4] = peptide_score; 
        if (simple_)
        {
            for(int i = 0; i < (int) score_vec.size(); i++)
            {
                score_vec[i] /= 100.0;
            }
        }else
        {
            for(int i = 0; i < (int) score_vec.size(); i++)
            {
                score_vec[i] = score_vec[i] * 1.0 / spectrum_;
            }
        }
           
        return score_vec;
    }

    void Emplace(int scan, const std::string& sequence, 
        const std::string& composite, int site, const std::vector<double>& score_vec)
    {
        SearchResult res;
        res.set_scan(scan);
        res.set_peptide(sequence);
        res.set_glycan(composite);
        res.set_site(site);
        res.set_score(score_vec);
        results_.push_back(res);
    }

    const int max_hits = 20;
    double best_ = 0.0;
    double spectrum_ = 0.0;
    double oxonium_ = 0.0;
    bool simple_ = false;
    std::map<int, double> peptide_;
    std::map<std::string, double> glycan_core_, glycan_branch_, glycan_terminal_;
    double precursor_mass_; 
    int isotopic_;
    std::vector<SearchResult> results_;

};



} // namespace engine
} // namespace search

#endif