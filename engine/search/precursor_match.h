#ifndef ENGINE_SEARCH_PRECURSOR_MATCH_H
#define ENGINE_SEARCH_PRECURSOR_MATCH_H

#include <string>
#include <vector>
#include <unordered_set>
#include "../../algorithm/search/search.h"
#include "../../util/mass/peptide.h"
#include "../../model/glycan/glycan.h"
#include "../../util/mass/glycan.h"
#include "../../util/mass/spectrum.h"
#include "../../engine/glycan/glycan_builder.h"
#include <iostream>

namespace engine{
namespace search{

class MatchResultStore
{
public:
    std::unordered_map<std::string, 
        std::unordered_set<std::string>> Map() const { return map_; }
    bool Empty() const { return peptides_.size() == 0; }
    std::vector<std::string> Peptides() const { return peptides_; }
    std::vector<std::string> Glycans() const
    {
        std::vector<std::string> res;
        for (const auto& peptide : Peptides())
        {
            if (map_.find(peptide) != map_.end())
            {
                std::unordered_set<std::string> glycans = Map()[peptide];
                res.insert(res.end(), glycans.begin(), glycans.end());
            }
        }
        return res;
    }
    std::unordered_set<std::string> Glycans(const std::string& peptide) const
    {
        if (map_.find(peptide) != map_.end())
        {
            return Map()[peptide];
        }
        return std::unordered_set<std::string>();
    }
    void Add(const std::string& peptide, const std::string& glycan)
    {
        if (map_.find(peptide) == map_.end())
        {
            peptides_.push_back(peptide);
            map_[peptide] = std::unordered_set<std::string>();
        }
        map_[peptide].insert(glycan);
    }

protected:
    std::vector<std::string> peptides_;
    std::unordered_map<std::string, 
        std::unordered_set<std::string>> map_;
};

class PrecursorMatcher
{
public:
    PrecursorMatcher(double tol, algorithm::search::ToleranceBy by, 
        engine::glycan::GlycanStore isomer): tolerance_(tol), by_(by),
            searcher_(algorithm::search::BasicSearch<std::string>(tol, by)),
                isomer_(isomer){}

    void Init(const std::vector<std::string>& peptides, const std::vector<std::string>& glycans)
    {
        // set up glycans
        set_glycans(glycans);
        // set up search box
        set_peptides(peptides);
    }

    std::vector<std::string>& Glycans() { return glycans_; }
    std::vector<std::string>& Peptides() { return peptides_; }
    virtual void set_glycans(const std::vector<std::string>& glycans) { glycans_ = glycans; }
    virtual void set_peptides(const std::vector<std::string>& peptides)
    {
        std::vector<std::shared_ptr<algorithm::search::Point<std::string>>> points;
        for(const auto& peptide : peptides)
        {
            double mass = util::mass::PeptideMass::Compute(peptide);
            std::shared_ptr<algorithm::search::Point<std::string>> p = 
                std::make_shared<algorithm::search::Point<std::string>>(mass, peptide);
            points.push_back(std::move(p));
        }
        searcher_.set_data(std::move(points));
        searcher_.Init();
    }

    double Tolerance() const { return tolerance_; }
    algorithm::search::ToleranceBy ToleranceType() const { return by_; }
    void set_tolerance(double tol) 
        { tolerance_ = tol; searcher_.set_tolerance(tol); searcher_.Init(); }
    void set_tolerance_by(algorithm::search::ToleranceBy by) 
        { by_ = by; searcher_.set_tolerance_by(by); searcher_.Init(); }

    virtual MatchResultStore Match(const double target, int charge)
    {
        return Match(target, charge, 0);
    }

    virtual MatchResultStore Match(const double target, int charge, const int isotope)
    {
        MatchResultStore res;
        if (searcher_.ToleranceType() == algorithm::search::ToleranceBy::PPM)
            searcher_.set_base(target);
        else if (searcher_.ToleranceType() == algorithm::search::ToleranceBy::Dalton)
            searcher_.set_scale(charge);

        for(const auto& glycan : glycans_)
        {
            double delta = target - isomer_.QueryMass(glycan);
            if (delta <= 0 ) continue;

            for (int i = 0; i <= isotope; i++)
            {
                double q = delta - i * util::mass::SpectrumMass::kIon;
                std::vector<std::string> peptides = searcher_.Query(q);
                for(const auto& peptide : peptides)
                {
                    res.Add(peptide, glycan);
                }
            }
        }
        return res;
    }

protected:
    double tolerance_;
    algorithm::search::ToleranceBy by_;
    algorithm::search::BasicSearch<std::string> searcher_;
    engine::glycan::GlycanStore isomer_;
    std::vector<std::string> glycans_;
    std::vector<std::string> peptides_;

}; 

} // namespace engine
} // namespace search

#endif