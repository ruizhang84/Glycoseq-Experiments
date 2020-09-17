#ifndef ENGINE_SEARCH_SPECTRUM_MATCH_H
#define ENGINE_SEARCH_SPECTRUM_MATCH_H

#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <algorithm>
#include "precursor_match.h"
#include "search_result.h"

#include "../../algorithm/search/bucket_search.h"
#include "../../algorithm/search/binary_search.h"
#include "../../util/mass/peptide.h"
#include "../../model/glycan/glycan.h"
#include "../../model/spectrum/spectrum.h"
#include "../../util/mass/glycan.h"
#include "../../util/mass/spectrum.h"
#include "../../util/mass/ion.h"
#include "../../engine/glycan/glycan_builder.h"
#include "../../engine/protein/protein_ptm.h"

#include <iostream>

namespace engine{
namespace search{

class SpectrumSearcher
{
public:
    SpectrumSearcher(const double tol, const algorithm::search::ToleranceBy by, int isotope,
        engine::glycan::NGlycanBuilder* builder, bool decoy_search):
            tolerance_(tol), by_(by), isotopic_(isotope), builder_(builder), decoy_search_(decoy_search),
                searcher_(algorithm::search::BucketSearch<model::spectrum::Peak>(tol, by)),
                binary_(algorithm::search::BinarySearch(tol, by)){}

    void Init()
    {
        glycan_isomer_ = builder_->Isomer();
        glycan_core_ = builder_->Core();
        glycan_branch_ = builder_->Branch();
        glycan_terminal_ = builder_->Terminal();
    }

    model::spectrum::Spectrum& Spectrum() { return spectrum_; }
    MatchResultStore& Candidate() { return candidate_; }
    void set_spectrum(const model::spectrum::Spectrum& spectrum) { spectrum_ = spectrum; }
    void set_candidate(const MatchResultStore& candidate) { candidate_ = candidate; }

    double Tolerance() const { return tolerance_; }
    algorithm::search::ToleranceBy ToleranceType() const { return by_; }
    int Isoptoic() const { return isotopic_; }
    void set_tolerance(double tol) 
        { tolerance_ = tol; searcher_.set_tolerance(tol); searcher_.Init(); }
    void set_tolerance_by(algorithm::search::ToleranceBy by) 
        { by_ = by; searcher_.set_tolerance_by(by); searcher_.Init(); }
    void set_isotopic(int isotope)
        { isotopic_ = isotope; }

    std::vector<SearchResult> Search()
    {
        SearchInit();
        ResultCollector collector;

        collector.OxoniumCollect(SearchOxonium());
        if (collector.OxoniumMiss()) 
            return collector.Result();

        collector.SpectrumBase(spectrum_.Peaks());
        for(const auto& peptide : candidate_.Peptides())
        {
            for(const auto& composite: candidate_.Glycans(peptide))
            {
                collector.InitCollect();
                for (const auto& pos : engine::protein::ProteinPTM::FindNGlycanSite(peptide))
                {
                   collector.PeptideCollect(SearchPeptides(peptide, composite, pos), pos);
                }
                if (collector.PeptideMiss()) continue;
                        

                std::unordered_map<std::string, double> result_core, result_branch, result_terminal;
                for(const auto & isomer : glycan_isomer_.Query(composite))
                {
                    collector.GlycanCollect(SearchGlycans(peptide, isomer, glycan_core_), 
                        isomer, SearchType::Core);
                    if (collector.GlycanMiss(isomer)) continue;

                    collector.GlycanCollect(SearchGlycans(peptide, isomer, glycan_branch_), 
                        isomer, SearchType::Branch);
                    collector.GlycanCollect(SearchGlycans(peptide, isomer, glycan_terminal_), 
                        isomer, SearchType::Terminal);
                }
                if (collector.GlycanMiss()) continue;
                          
                if (decoy_search_)
                    collector.Update(spectrum_.Scan(), peptide, composite);
                else
                    collector.BestUpdate(spectrum_.Scan(), peptide, composite);
            }
        }
        if (collector.Empty())
            return collector.Result();
            
        // compute precursor differ
        double precursor_mass = 
            util::mass::SpectrumMass::Compute(spectrum_.PrecursorMZ(), spectrum_.PrecursorCharge());
        collector.PrecursorCollect(precursor_mass, isotopic_);
        
        // save 
        if (decoy_search_)
            return collector.Result();
        return collector.BestResult();   
    }

protected:
    void SearchInit()
    {
        std::vector<std::shared_ptr<algorithm::search::Point<model::spectrum::Peak>>> mz_points;
        for(const auto& it : spectrum_.Peaks())
        {
            std::shared_ptr<algorithm::search::Point<model::spectrum::Peak>> p = 
                std::make_shared<algorithm::search::Point<model::spectrum::Peak>>(it.MZ(), it);
            mz_points.push_back(std::move(p));
        }
                    
        searcher_.set_data(std::move(mz_points));
        searcher_.Init();
    }

    std::vector<model::spectrum::Peak> SearchOxonium()
    {
        std::vector<model::spectrum::Peak> res;
        for (const auto& mass : oxonium_)
        {
            for(int charge = 1; charge <= spectrum_.PrecursorCharge(); charge++)
            {
                double mz = util::mass::SpectrumMass::ComputeMZ(mass, charge);
                std::vector<model::spectrum::Peak> p = searcher_.Query(mz);
                if (! p.empty())
                {
                    res.push_back(*std::max_element(p.begin(), p.end(), IntensityCmp));
                }
            }
        }
        return res;
    }

    std::vector<model::spectrum::Peak> SearchPeptides
        (const std::string& seq, const std::string& composite, const int pos)
    {
        std::vector<model::spectrum::Peak> res;
        std::vector<double> peptides_mz;
       
        // speed up
        std::string key = seq + std::to_string(pos);
        if (peptides_ptm_mz_.find(key) == peptides_ptm_mz_.end())
        {
            peptides_ptm_mz_[key] = ComputePTMPeptideMass(seq, pos);
            std::sort(peptides_ptm_mz_[key].begin(), peptides_ptm_mz_[key].end());

            peptides_mz_[key] = ComputeNonePTMPeptideMass(seq, pos);
            std::sort(peptides_mz_[key].begin(), peptides_mz_[key].end());
        }

        // search ptm
        binary_.set_data(peptides_ptm_mz_[key]);
        double extra = util::mass::GlycanMass::Compute(model::glycan::Glycan::Interpret(composite));
        for(const auto& pk : spectrum_.Peaks())
        {
            for (int charge = 1; charge <= spectrum_.PrecursorCharge(); charge++)
            {
                double target = util::mass::SpectrumMass::Compute(pk.MZ(), charge);
                if (binary_.ToleranceType() == algorithm::search::ToleranceBy::PPM)
                    binary_.set_base(target);
                else if (binary_.ToleranceType() == algorithm::search::ToleranceBy::Dalton)
                    binary_.set_scale(charge);
                if (target > extra && binary_.Search(target-extra))
                {
                    res.push_back(pk);
                    break;
                }
            }
        }

        // search peptides
        binary_.set_data(peptides_mz_[key]);
        for(const auto& pk : spectrum_.Peaks())
        {
            for (int charge = 1; charge <= spectrum_.PrecursorCharge(); charge++)
            {
                double target = util::mass::SpectrumMass::Compute(pk.MZ(), charge);
                if (binary_.ToleranceType() == algorithm::search::ToleranceBy::PPM)
                    binary_.set_base(target);
                else if (binary_.ToleranceType() == algorithm::search::ToleranceBy::Dalton)
                    binary_.set_scale(charge);
                if (binary_.Search(target))
                {
                    res.push_back(pk);
                    break;
                }
            }
        }
        return res;
    }

    std::vector<model::spectrum::Peak> SearchGlycans
        (const std::string& seq, const std::string& id, 
        engine::glycan::GlycanMassStore& glycan_mass_)
    {
        std::vector<model::spectrum::Peak> res;
        std::unordered_set<double> subset = glycan_mass_.Query(id);
        std::vector<double> subset_mass;
        subset_mass.insert(subset_mass.end(), subset.begin(), subset.end());
        binary_.set_data(subset_mass);
        binary_.Init();

        double extra = util::mass::PeptideMass::Compute(seq);
        for(const auto& pk : spectrum_.Peaks())
        {
            for(int charge = 1; charge <= spectrum_.PrecursorCharge(); charge++)
            {
                double mass = util::mass::SpectrumMass::Compute(pk.MZ(), charge);
                if (binary_.ToleranceType() == algorithm::search::ToleranceBy::PPM)
                    binary_.set_base(mass);
                else if (binary_.ToleranceType() == algorithm::search::ToleranceBy::Dalton)
                    binary_.set_scale(charge);

                if (mass > extra && binary_.Search(mass-extra))
                {        
                    res.push_back(pk);
                    break;
                }
            }
        }
        return res;
    }

    // for computing the peptide ions
    static std::vector<double> ComputePTMPeptideMass(const std::string& seq, const int pos)
    {
        std::vector<double> mass_list;
        for (int i = pos; i < (int) seq.length() - 1; i++) // seldom at n
        {
            double mass = util::mass::IonMass::Compute(seq.substr(0, i+1), util::mass::IonType::b);
            mass_list.push_back(mass);
            mass = util::mass::IonMass::Compute(seq.substr(0, i+1), util::mass::IonType::c);
            mass_list.push_back(mass);
        }
        for (int i = 1; i <= pos; i++)
        {
            double mass = util::mass::IonMass::Compute(seq.substr(i, seq.length()-i), util::mass::IonType::y);
            mass_list.push_back(mass);
            mass = util::mass::IonMass::Compute(seq.substr(i, seq.length()-i), util::mass::IonType::z);
            mass_list.push_back(mass);
        }
        return mass_list;
    }

    static std::vector<double> ComputeNonePTMPeptideMass(const std::string& seq, const int pos)
    {
        std::vector<double> mass_list;
        for (int i = 0; i < pos; i++) // seldom at n
        {
            double mass = util::mass::IonMass::Compute(seq.substr(0, i+1), util::mass::IonType::b);
            mass_list.push_back(mass);
            mass = util::mass::IonMass::Compute(seq.substr(0, i+1), util::mass::IonType::c);
            mass_list.push_back(mass);
        }
        for (int i = pos + 1; i < (int) seq.length(); i++)
        {
            double mass = util::mass::IonMass::Compute(seq.substr(i, seq.length()-i), util::mass::IonType::y);
            mass_list.push_back(mass);
            mass = util::mass::IonMass::Compute(seq.substr(i, seq.length()-i), util::mass::IonType::z);
            mass_list.push_back(mass);
        }
        return mass_list;
    }


    static bool IntensityCmp(const model::spectrum::Peak& p1, const model::spectrum::Peak& p2)
        { return (p1.Intensity() < p2.Intensity()); }

    double tolerance_;
    algorithm::search::ToleranceBy by_;
    int isotopic_; // up to isotopic
    engine::glycan::NGlycanBuilder* builder_;
    bool decoy_search_;
    algorithm::search::BucketSearch<model::spectrum::Peak> searcher_;
    algorithm::search::BinarySearch binary_;
    MatchResultStore candidate_;
    model::spectrum::Spectrum spectrum_;
    std::unordered_map<std::string, std::vector<double>> peptides_ptm_mz_;
    std::unordered_map<std::string, std::vector<double>> peptides_mz_; 

    engine::glycan::GlycanStore glycan_isomer_;
    engine::glycan::GlycanMassStore glycan_core_, glycan_branch_, glycan_terminal_;

    const std::vector<double> oxonium_ 
    {
        util::mass::GlycanMass::kHexNAc,
        util::mass::GlycanMass::kHexNAc - util::mass::GlycanMass::kWater,
        util::mass::GlycanMass::kHexNAc - util::mass::GlycanMass::kWater * 2,
        util::mass::GlycanMass::kHexNAc + util::mass::GlycanMass::kHex
    };
}; 

} // namespace engine
} // namespace search

#endif