#ifndef APP_SEARCH_WORK_DISTRIBUTOR_H
#define APP_SEARCH_WORK_DISTRIBUTOR_H

#include <deque>
#include <thread>  
#include <mutex> 

#include "search_parameter.h"
#include "../../engine/spectrum/normalize.h"
#include "../../engine/search/spectrum_search.h"

class SearchQueue
{
public:
    SearchQueue(const std::vector<model::spectrum::Spectrum>& spectra)
        { GenerateQueue(spectra); }

    SearchQueue(const SearchQueue& other)
    {
        queue_ = other.queue_;
    }

    virtual void GenerateQueue(
        std::vector<model::spectrum::Spectrum> spectra)
    {
        for(const auto& it : spectra)
        {
            queue_.push_back(it);
        }
    }

    virtual model::spectrum::Spectrum TryGetSpectrum()
    {
        model::spectrum::Spectrum spec;
        mutex_.lock();
            if (! queue_.empty())
            {
                spec = queue_.front();
                queue_.pop_front();
            }
            else
            {
                spec.set_scan(-1);
            }
        mutex_.unlock();
        return spec;
    }
    


protected:
    std::deque<model::spectrum::Spectrum> queue_;
    std::mutex mutex_; 
};



class SearchDispatcher
{
public:
    SearchDispatcher(const std::vector<model::spectrum::Spectrum>& spectra, 
        engine::glycan::NGlycanBuilder* builder, const std::vector<std::string>& peptides, 
            SearchParameter parameter): queue_(SearchQueue(spectra)), builder_(builder), 
                peptides_(peptides), parameter_(parameter){}

    engine::glycan::NGlycanBuilder* Builder() { return builder_; }
    std::vector<std::string> Peptides() { return peptides_; }
    SearchParameter Parameter() { return parameter_; }
    void set_builder(engine::glycan::NGlycanBuilder* builder)
        { builder_ = builder; }
    void set_peptides(std::vector<std::string> peptides) 
        { peptides_ = peptides; }
    void set_parameter(SearchParameter parameter) 
        { parameter_ = parameter; }

    void set_score_compute(bool simple){
        simple_ = simple;
    }    

    std::vector<engine::search::SearchResult> Dispatch()
    {
        std::vector<engine::search::SearchResult> results;
        std::vector< std::thread> thread_pool;
        for (int i = 0; i < parameter_.n_thread; i ++)
        {
            std::thread worker(&SearchDispatcher::SearchingWorker, this, std::ref(results), false);
            thread_pool.push_back(std::move(worker));
        }
        for (auto& worker : thread_pool)
        {
            worker.join();
        }
        return results;
    }

    std::vector<engine::search::SearchResult> DecoyDispatch()
    {
        std::vector<engine::search::SearchResult> results;
        std::vector< std::thread> thread_pool;
        for (int i = 0; i < parameter_.n_thread; i ++)
        {
            std::thread worker(&SearchDispatcher::SearchingWorker, this, std::ref(results), true);
            thread_pool.push_back(std::move(worker));
        }
        for (auto& worker : thread_pool)
        {
            worker.join();
        }
        return results;
    }

protected:
    void SearchingWorker(
        std::vector<engine::search::SearchResult>& results, bool decoy_search)
    {
        engine::search::PrecursorMatcher precursor_runner
            (parameter_.ms1_tol, parameter_.ms1_by, builder_->Isomer());
        engine::search::SpectrumSearcher spectrum_runner
            (parameter_.ms2_tol, parameter_.ms2_by, parameter_.isotopic_count, builder_, decoy_search, parameter_.weight);
        std::vector<std::string> glycans_str = builder_->Isomer().Collection();
        precursor_runner.Init(peptides_, glycans_str);
        spectrum_runner.Init();
        spectrum_runner.set_score_compute(simple_);

        std::vector<engine::search::SearchResult> temp_result;
        
        while (true)
        {
            model::spectrum::Spectrum spec = queue_.TryGetSpectrum();
            if (spec.Scan() < 0) break;
            
            // precusor
            double target = 
                util::mass::SpectrumMass::Compute(spec.PrecursorMZ(), spec.PrecursorCharge());
            engine::search::MatchResultStore r = 
                precursor_runner.Match(target, spec.PrecursorCharge(), parameter_.isotopic_count);
            if (r.Empty()) continue;

            // process spectrum by normalization
            engine::spectrum::Normalizer::Transform(spec);

            // msms
            spectrum_runner.set_spectrum(spec);
            spectrum_runner.set_candidate(r);
            std::vector<engine::search::SearchResult> res = spectrum_runner.Search();
            if (res.empty()) continue;

            temp_result.insert(temp_result.end(), res.begin(), res.end());
        }
        
        mutex_.lock();
            results.insert(results.end(), temp_result.begin(), temp_result.end());
        mutex_.unlock();
    }

    std::mutex mutex_; 
    SearchQueue queue_;
    engine::glycan::NGlycanBuilder* builder_;
    std::vector<std::string> peptides_;
    SearchParameter parameter_;
    bool simple_ = false;

};

#endif