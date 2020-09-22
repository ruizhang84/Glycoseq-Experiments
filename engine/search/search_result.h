#ifndef ENGINE_SEARCH_SEARCH_SCORE_H
#define ENGINE_SEARCH_SEARCH_SCORE_H


#include <vector>
#include <map>
#include <unordered_map>
#include <cmath> 
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
    int Scan() const { return scan_; }
    int ModifySite() const { return pos_; }
    std::string Sequence() const { return peptide_; }
    std::string Glycan() const { return glycan_; }
    double Score() const { return score_; }
    double ExtraScore(ScoreType type) const 
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
    void set_score(double score) { score_ = score; }
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
    int scan_;
    std::string peptide_;
    std::string glycan_;
    int pos_;
    double score_;
    std::map<ScoreType, double> extra_;
    
};


class ResultCollector{
public:
    ResultCollector(): best_(0.0), oxonium_(0){}

    void set_score_compute(bool simple){
        simple_ = simple;
    }
    
    std::vector<SearchResult> Result()
    {
        // remove if hits over 20
        if ((int) results_.size() > max_hits)
        {
            std::sort(results_.begin(), results_.end(), 
                [](const SearchResult& r1, const SearchResult& r2) -> bool { return r1.Score() > r1.Score(); });
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
            double score = it.Score();
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
            double score = ComputeScore(pos_it.second);
            // emplace results
            Emplace(scan, sequence, composite, pos_it.first, score);
        }
    }

    void BestUpdate(int scan, const std::string& sequence, const std::string& composite)
    {
        for(const auto& pos_it : peptide_)
        {
            // compute score
            double score = ComputeScore(pos_it.second);
            if (score >= best_)
            {
                if (score > best_)
                    results_.clear();
                best_ = score;
                // emplace results
                Emplace(scan, sequence, composite, pos_it.first, score);
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
                    glycan_branch_[isomer] =  SearchResult::PeakValue(glycan_peaks, simple_);
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
    double ComputeScore(double peptide_score)
    {
        double score = 0;
        for(const auto& isomer_it : glycan_core_)
        {
            std::string isomer = isomer_it.first;
            double glycan_score = glycan_core_[isomer] + glycan_branch_[isomer] + glycan_terminal_[isomer]; 
            score = std::max(score, glycan_score);  
        }
        score += peptide_score + oxonium_;
        if (!simple_)
            score = std::sqrt(score) * 1.0 / std::sqrt(spectrum_);
        return score;
    }

    void Emplace(int scan, const std::string& sequence, 
        const std::string composite, int site, double score)
    {
        SearchResult res;
        res.set_scan(scan);
        res.set_peptide(sequence);
        res.set_glycan(composite);
        res.set_site(site);
        res.set_score(score);
        results_.push_back(res);
    }

    const int max_hits = 20;
    double best_;
    double spectrum_;
    double oxonium_;
    bool simple_;
    std::unordered_map<int, double> peptide_;
    std::unordered_map<std::string, double> glycan_core_, glycan_branch_, glycan_terminal_;
    double precursor_mass_; 
    int isotopic_;
    std::vector<SearchResult> results_;

};



} // namespace engine
} // namespace search

#endif